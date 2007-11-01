`poLCA` <-
function(formula,data,nclass=2,maxiter=1000,graphs=FALSE,tol=1e-10,na.rm=TRUE,probs.start=NULL,nrep=1,verbose=TRUE) {
    starttime <- Sys.time()
    if (!na.rm) {
        mframe <- model.frame(formula,data,na.action=NULL)
        y <- model.response(mframe)
        x <- model.matrix(formula,mframe)
        y[is.na(y)] <- 0
        data <- data.frame(y,x) 
    }
    mframe <- model.frame(formula,data)
    y <- model.response(mframe)
    x <- model.matrix(formula,mframe)
    N <- nrow(y)
    J <- ncol(y)
    K.j <- t(matrix(apply(y,2,max)))
    R <- nclass
    S <- ncol(x)
    eflag <- FALSE
    probs.start.ok <- TRUE
    ret <- list()
    if (R==1) {
        ret$probs <- list()
        for (j in 1:J) { ret$probs[[j]] <- matrix(table(y[,j])/sum(table(y[,j])),nrow=1) }
        ret$probs.start <- ret$probs
        ret$P <- 1
        ret$posterior <- ret$predclass <- prior <- matrix(1,nrow=N,ncol=1)
        ret$llik <- sum(log(poLCA.ylik.C(poLCA.vectorize(ret$probs),y)))
        se <- poLCA.se(y,x,ret$probs,prior,ret$posterior)
        ret$probs.se <- se$probs           # standard errors of class-conditional response probabilities
        ret$P.se <- se$P                   # standard errors of class population shares
        ret$numiter <- 1
        ret$probs.start.ok <- TRUE
        ret$coeff <- NULL
        ret$coeff.se <- NULL
        ret$coeff.V <- NULL
        ret$eflag <- FALSE
        if (S>1) {
            cat("\n ALERT: covariates not allowed when nclass=1; will be ignored. \n \n")
            S <- 1
        }
    } else {
        if (graphs) ifelse(max(K.j)==2,
                            layout(matrix(c(1,2),2,1),heights=c(8,1)),
                            layout(matrix(seq(1,(R+1)),R+1,1),heights=c(rep(5,R),1)))
        if (!is.null(probs.start)) { # error checking on user-inputted probs.start
            if ((length(probs.start) != J) | (!is.list(probs.start))) {
                probs.start.ok <- FALSE
            } else {
                if (sum(sapply(probs.start,dim)[1,]==R) != J) probs.start.ok <- FALSE
                if (sum(sapply(probs.start,dim)[2,]==K.j) != J) probs.start.ok <- FALSE
                if (sum(round(sapply(probs.start,rowSums),4)==1) != (R*J)) probs.start.ok <- FALSE
            }
        }
        ret$llik <- -Inf
        ret$attempts <- NULL
        for (repl in 1:nrep) { # automatically reestimate the model multiple times to locate the global max llik
            error <- TRUE; firstrun <- TRUE
            probs <- probs.init <- probs.start
            while (error) { # error trap
                error <- FALSE
                b <- rep(0,S*(R-1))
                prior <- poLCA.updatePrior(b,x,R)
                if ((!probs.start.ok) | (is.null(probs.start)) | (!firstrun) | (repl>1)) { # only use the specified probs.start in the first nrep
                    probs <- list()
                    for (j in 1:J) { 
                        probs[[j]] <- matrix(runif(R*K.j[j]),nrow=R,ncol=K.j[j])
                        probs[[j]] <- probs[[j]]/rowSums(probs[[j]]) 
                    }
                    probs.init <- probs
                }
                vp <- poLCA.vectorize(probs)
                iter <- 1
                llik <- matrix(NA,nrow=maxiter,ncol=1)
                llik[iter] <- -Inf
                dll <- Inf
                while ((iter <= maxiter) & (dll > tol) & (!error)) {
                    iter <- iter+1
                    rgivy <- poLCA.postClass.C(prior,vp,y)      # calculate posterior
                    vp$vecprobs <- poLCA.probHat.C(rgivy,y,vp)  # update probs
                    if (S>1) {
                        dd <- poLCA.dLL2dBeta.C(rgivy,prior,x)
                        b <- b + ginv(-dd$hess) %*% dd$grad     # update betas
                        prior <- poLCA.updatePrior(b,x,R)       # update prior
                    } else {
                        prior <- matrix(colMeans(rgivy),nrow=N,ncol=R,byrow=TRUE)
                    }
                    llik[iter] <- sum(log(rowSums(prior*poLCA.ylik.C(vp,y))))
                    dll <- llik[iter]-llik[iter-1]
                    if (is.na(dll)) {
                        error <- TRUE
                    } else if ((S>1) & (dll < -1e-7)) {
                        error <- TRUE
                    }
                    if (graphs) {
                        if (max(K.j)==2) {
                            poLCA.makeplot.dich(poLCA.unvectorize(vp),colMeans(rgivy),y,NULL)
                        } else {
                            for (r in 1:R) { poLCA.makeplot.poly(poLCA.unvectorize(vp),r,y,K.j,paste("Class",r,": p=",round(colMeans(rgivy)[r],3))) }
                        }
                        par(mar=c(0,0,0,0))
                        plot(0,main=paste("\n iteration",iter,": log-lik =",llik[iter]),cex.main=1.5,
                             col="white",col.axis="white",col.lab="white",xaxt="n",yaxt="n",bty="n")
                    }
                }
                if (!error) { 
                    se <- poLCA.se(y,x,poLCA.unvectorize(vp),prior,rgivy)
                } else {
                    eflag <- TRUE
                }
                firstrun <- FALSE
            } # finish estimating model without triggering error
            ret$attempts <- c(ret$attempts,llik[iter])
            if (llik[iter] > ret$llik) {
                ret$llik <- llik[iter]             # maximum value of the log-likelihood
                ret$probs.start <- probs.init      # starting values of class-conditional response probabilities
                ret$probs <- poLCA.unvectorize(vp) # estimated class-conditional response probabilities
                ret$probs.se <- se$probs           # standard errors of class-conditional response probabilities
                ret$P.se <- se$P                   # standard errors of class population shares
                ret$posterior <- rgivy             # NxR matrix of posterior class membership probabilities
                ret$predclass <- apply(ret$posterior,1,which.max)   # Nx1 vector of predicted class memberships, by modal assignment
                ret$P <- colMeans(ret$posterior)   # estimated class population shares
                ret$numiter <- iter-1              # number of iterations until reaching convergence
                ret$probs.start.ok <- probs.start.ok # if starting probs specified, logical indicating proper entry format
                if (S>1) {
                    b <- matrix(b,nrow=S)
                    rownames(b) <- colnames(x)
                    rownames(se$b) <- colnames(x)
                    ret$coeff <- b                 # coefficient estimates (when estimated)
                    ret$coeff.se <- se$b           # standard errors of coefficient estimates (when estimated)
                    ret$coeff.V <- se$var.b        # covariance matrix of coefficient estimates (when estimated)
                } else {
                    ret$coeff <- NULL
                    ret$coeff.se <- NULL
                    ret$coeff.V <- NULL
                }
                ret$eflag <- eflag                 # error flag, true if estimation algorithm ever needed to restart with new initial values
            }
            if (nrep>1) { cat("Model ",repl,": llik = ",llik[iter]," ... best llik = ",ret$llik,"\n",sep=""); flush.console() }
        } # end replication loop
    }
    names(ret$probs) <- names(ret$probs.se) <- colnames(y)
    ret$npar <- (R*sum(K.j-1)) + (R-1)                  # number of degrees of freedom used by the model (number of estimated parameters)
    if (S>1) { ret$npar <- ret$npar + (S*(R-1)) - (R-1) }
    ret$aic <- (-2 * ret$llik) + (2 * ret$npar)         # Akaike Information Criterion
    ret$bic <- (-2 * ret$llik) + (log(N) * ret$npar)    # Schwarz-Bayesian Information Criterion
    ret$Nobs <- sum(rowSums(y==0)==0)                   # number of fully observed cases (if na.rm=F)
    if (all(rowSums(y==0)>0)) { # if no rows are fully observed
        ret$Chisq <- NA
        ret$Gsq <- NA
        ret$predcell <- NULL
    } else {
        compy <- poLCA.compress(y[(rowSums(y==0)==0),])
        datacell <- compy$datamat
        rownames(datacell) <- NULL
        freq <- compy$freq
        if (!na.rm) {
            fit <- matrix(ret$Nobs * (poLCA.ylik.C(poLCA.vectorize(ret$probs),datacell) %*% ret$P))
            ret$Chisq <- sum((freq-fit)^2/fit) + (ret$Nobs-sum(fit)) # Pearson Chi-square goodness of fit statistic for fitted vs. observed multiway tables
        } else {
            fit <- matrix(N * (poLCA.ylik.C(poLCA.vectorize(ret$probs),datacell) %*% ret$P))
            ret$Chisq <- sum((freq-fit)^2/fit) + (N-sum(fit))
        }
        ret$predcell <- data.frame(datacell,observed=freq,expected=round(fit,3)) # Table that gives observed vs. predicted cell counts
        ret$Gsq <- 2 * sum(freq*log(freq/fit))  # Likelihood ratio/deviance statistic
    }
    y[y==0] <- NA
    for (j in 1:J) { dimnames(ret$probs[[j]]) <- list(paste("class ",1:R,": ",sep=""),
                                                  paste("Pr(",1:ncol(ret$probs[[j]]),")",sep="")) }
    ret$y <- data.frame(y)             # outcome variables
    ret$x <- data.frame(x)             # covariates, if specified
    ret$N <- N                         # number of observations
    ret$maxiter <- maxiter             # maximum number of iterations specified by user
    ret$resid.df <- min(ret$N,(prod(K.j)-1))-ret$npar # number of residual degrees of freedom
    class(ret) <- "poLCA"
    if (verbose) print.poLCA(ret)
    ret$time <- Sys.time()-starttime   # how long it took to run the model
    return(ret)
}


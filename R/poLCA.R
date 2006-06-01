"poLCA" <-
function(formula,data,nclass=2,maxiter=1000,graphs=FALSE,tol=1e-10,na.rm=TRUE) {
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
    Nobs <- sum(rowSums(y==0)==0)
    J <- ncol(y)
    K.j <- t(matrix(apply(y,2,max)))
    R <- nclass
    S <- ncol(x)
    probs <- list()
    eflag <- FALSE
    if (R==1) {
        for (j in 1:ncol(y)) {
            probs[[j]] <- matrix(NA,nrow=1,ncol=K.j[j])
            for (k in 1:K.j[j]) { probs[[j]][k] <- sum(y[,j]==k)/sum(table(y[,j])) }
        }
        vp <- poLCA.vectorize(probs)
        P <- 1
        rgivy <- matrix(1,nrow=N,ncol=1)
        ml <- sum(log(poLCA.ylik.C(vp,y)))
        iter <- 1
        if (S>1) {
            cat("\n ALERT: covariates not allowed when nclass=1; will be ignored. \n \n")
            S <- 0
        }
    } else {
        if (graphs) ifelse(max(K.j)==2,layout(matrix(c(1,2),2,1),heights=c(3.2,1)),layout(matrix(seq(1,(R+1)),R+1,1),heights=c(rep(2,R),1)))
        error <- T
        while (error) { # error trap
            error <- F
            b <- rep(0,S*(R-1))
            prior <- poLCA.updatePrior(b,x,R)
            for (j in 1:J) { 
                probs[[j]] <- matrix(runif(R*K.j[j]),nrow=R,ncol=K.j[j])
                probs[[j]] <- probs[[j]]/rowSums(probs[[j]]) 
            }
            vp <- poLCA.vectorize(probs)
            iter <- 1
            llik <- matrix(NA,nrow=maxiter,ncol=1)
            llik[iter] <- -Inf
            dll <- Inf
            while ((iter <= maxiter) & (dll > tol) & (!error)) {
                iter <- iter+1
                rgivy <- poLCA.postClass.C(prior,vp,y)
                vp$vecprobs <- poLCA.probHat.C(rgivy,y,vp)
                if (S>1) {
                    dd <- poLCA.dLL2dBeta.C(rgivy,prior,x)
                    if (!is.matrix(try(solve(-dd$hess),silent=TRUE))) {
                        error <- T
                    } else {
                        b <- b + solve(-dd$hess) %*% dd$grad
                    }
                    prior <- poLCA.updatePrior(b,x,R)
                } else {
                    prior <- matrix(colMeans(rgivy),nrow=N,ncol=R,byrow=TRUE)
                }
                llik[iter] <- sum(log(rowSums(prior*poLCA.ylik.C(vp,y))))
                dll <- llik[iter]-llik[iter-1]
                if (is.na(dll)) {
                    error <- T
                } else if ((S>1) & (dll < -1e-7)) {
                    error <- T
                }
                if (graphs) {
                    if (max(K.j)==2) {
                        poLCA.makeplot.dich(poLCA.unvectorize(vp),colMeans(rgivy),y,NULL)
                    } else {
                        for (r in 1:R) { poLCA.makeplot.poly(poLCA.unvectorize(vp),r,y,K.j,paste("Class",r,": p=",round(colMeans(rgivy)[r],3))) }
                    }
                    par(mar=c(2.8,2.5,1.5,1),mgp=c(1.5,0.5,0))
                    plot(llik,type="l",xlab="iteration",ylab="log-likelihood",main=paste("iteration",iter,": log-lik =",llik[iter]),
                        lwd=2,cex.axis=0.9,cex.lab=0.9)
                }
            }
            if ((!error) & (S>1)){
                if (sum(diag(solve(-dd$hess))<0)>0) error <- T
            }
            P <- colMeans(rgivy)
            ml <- llik[iter]
            if (error) eflag <- T
        }
    }
    cl.pred <- apply(rgivy,1,which.max)
    npar <- (R*sum(K.j-1)) + (R-1)
    if (S>1) npar <- npar + (S*(R-1)) - (R-1)
    aic <- (-2 * ml) + (2 * npar)
    bic <- (-2 * ml) + (log(N) * npar)
    compy <- poLCA.compress(y[(rowSums(y==0)==0),])
    datacell <- compy$datamat
    rownames(datacell) <- NULL
    freq <- compy$freq
    if (!na.rm) {
        fit <- matrix(Nobs * (poLCA.ylik.C(vp,datacell) %*% P))
        Chisq <- sum((freq-fit)^2/fit) + (Nobs-sum(fit))
    } else {
        fit <- matrix(N * (poLCA.ylik.C(vp,datacell) %*% P))
        Chisq <- sum((freq-fit)^2/fit) + (N-sum(fit))
    }
    predcell <- data.frame(datacell,observed=freq,expected=round(fit,3))
    Gsq <- 2 * sum(freq*log(freq/fit))
    resid.df <- min(N,(prod(K.j)-1))-npar
    cat("Conditional item response (column) probabilities,\n by outcome variable, for each class (row) \n \n")
    probs <- poLCA.unvectorize(vp)
    names(probs) <- colnames(y)
    print(probs)
    cat("Estimated class population shares \n", P, "\n \n")
    cat("Predicted class memberships (by modal posterior prob.) \n",table(cl.pred)/N, "\n \n")
    cat("=========================================== \n")
    cat("Fit for", R, "latent classes: \n")
    cat("=========================================== \n")
    if (S>1) {
        b <- matrix(b,nrow=S)
        #se <- matrix(sqrt(diag(solve(-dd$hess))),nrow=S)
        se <- NULL
        for (r in 2:R) {
            cat(r,"/ 1 \n")
            #disp <- cbind(b[,(r-1)]),se[,(r-1)])
            disp <- matrix(b[,(r-1)],ncol=1)
            rownames(disp) <- colnames(x)
            colnames(disp) <- c("coefficient")#,"standard error")
            print(disp)
            cat("=========================================== \n")
        }
    }
    cat("number of observations:", N, "\n")
    if(!na.rm) cat("number of fully observed cases:", Nobs, "\n")
    cat("number of estimated parameters:", npar, "\n")
    cat("residual degrees of freedom:", resid.df, "\n")
    cat("maximum log-likelihood:", ml, "\n \n")
    cat("AIC(",R,"):",aic,"\n")
    cat("BIC(",R,"):",bic,"\n")
    if (S==1) cat("G^2(",R,"):",Gsq,"(Likelihood ratio/deviance statistic) \n")
    cat("X^2(",R,"):",Chisq,"(Chi-square goodness of fit) \n \n")
    if ((iter-1)==maxiter) cat("iterations finished, MAXIMUM LIKELIHOOD NOT FOUND \n \n")
    if (npar>N) cat("ALERT: number of parameters estimated (",npar,") exceeds number of observations (",N,") \n \n")
    if (resid.df<0) cat("ALERT: negative degrees of freedom; respecify model \n \n")
    if (eflag) cat("ALERT: estimation algorithm automatically restarted with new initial values \n \n")
    ret <- list()
    ret$y <- data.frame(y)             # outcome variables
    ret$x <- data.frame(x)             # covariates, if specified
    ret$N <- N                         # number of observations
    ret$Nobs <- Nobs                   # number of fully observed cases (if na.rm=F)
    ret$probs <- probs                 # estimated class-conditional "response" probabilities
    ret$P <- P                         # estimated class population shares
    ret$posterior <- rgivy             # NxR matrix of posterior class membership probabilities
    ret$predclass <- cl.pred           # Nx1 vector of predicted class memberships, by modal assignment
    ret$predcell <- predcell           # Table that gives observed vs. predicted cell counts
    ret$llik <- ml                     # maximum value of the log-likelihood
    ret$numiter <- iter-1              # number of iterations until reaching convergence
    if (S>1) {
        rownames(b) <- colnames(x)
        #rownames(se) <- colnames(x)
        ret$coeff <- b                 # coefficient estimates (when estimated)
        ret$coeff.se <- se             # standard errors of coefficient estimates (when estimated)
    } else {
        ret$coeff <- NULL
        ret$coeff.se <- NULL
    }
    ret$aic <- aic                     # Akaike Information Criterion
    ret$bic <- bic                     # Schwarz-Bayesian Information Criterion
    ret$Gsq <- Gsq                     # Likelihood ratio/deviance statistic
    ret$Chisq <- Chisq                 # Pearson Chi-square goodness of fit statistic for fitted vs. observed multiway tables
    ret$time <- Sys.time()-starttime   # how long it took to run the model
    ret$npar <- npar                   # number of degrees of freedom used by the model (number of estimated parameters)
    ret$resid.df <- resid.df           # number of residual degrees of freedom
    ret$eflag <- eflag                 # error flag, true if estimation algorithm ever needed to restart with new initial values
    return(ret)
}


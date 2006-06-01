"poLCA.simdata" <-
function(N=5000,probs=NULL,nclass=2,ndv=4,nresp=NULL,niv=0,b=NULL,classdist=NULL,missval=FALSE,pctmiss=NULL) {
    if (is.null(probs)) {
        if (is.null(nresp)) nresp <- ceiling(runif(ndv,min=1,max=5))
        if (!is.null(classdist)) nclass <- length(classdist)
        ndv <- length(nresp)
        probs <- list()
        for(i in 1:ndv) {
            probs[[i]] <- matrix(runif(nclass*nresp[i]),nrow=nclass,ncol=nresp[i])
            probs[[i]] <- probs[[i]]/rowSums(probs[[i]])
        }
    } else {
        ndv <- length(probs)
        nclass <- nrow(probs[[1]])
        nresp <- ncol(probs[[1]])
        for (i in 2:ndv) {
          nresp <- c(nresp,ncol(probs[[i]]))
        }
    }
    if (nclass==1) niv <- 0
    if (niv > 0) {
        x <- matrix(rnorm(N*niv),nrow=N,ncol=niv)
        colnames(x) <- paste("X",c(1:niv),sep="")
        if (is.null(b)) b <- matrix(round(runif(((nclass-1)*(niv+1)),min=-2,max=2)),nrow=(niv+1))
        prior <- poLCA.updatePrior(b,cbind(1,x),nclass)
    } else {
        x <- NULL
        if (nrow(probs[[1]])!=length(classdist)) {
            classdist <- runif(nclass)
            classdist <- classdist/sum(classdist)
        }
        prior <- matrix(classdist,byrow=TRUE,nrow=N,ncol=nclass)
    }
    ifelse(ncol(prior)>1,group <- rmulti(prior),group <- matrix(1,nrow=N,ncol=1))
    y <- rmulti(probs[[1]][group,])
    for (j in 2:ndv) {
        y <- cbind(y,rmulti(probs[[j]][group,]))
    }
    colnames(y) <- paste("Y",c(1:ndv),sep="")
    if (niv > 0) classdist <- colMeans(poLCA.postClass.C(prior,poLCA.vectorize(probs),y))
    if (missval) {
        if (is.null(pctmiss)) pctmiss <- runif(1,min=0.05,max=0.4)
        make.na <- cbind(ceiling(runif(round(pctmiss*N*ndv),min=0,max=N)),ceiling(runif(round(pctmiss*N*ndv),min=0,max=ndv)))
        y[make.na] <- NA
    }
    ret <- list()
    if (is.null(x)) {
        ret$dat <- data.frame(y)
    } else {
        ret$dat <- data.frame(y,x)
    }
    ret$trueclass <- group
    ret$probs <- probs
    ret$nresp <- nresp
    ret$b <- b
    ret$classdist <- classdist
    ret$pctmiss <- pctmiss
    return(ret)
}


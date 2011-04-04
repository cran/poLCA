poLCA.table <-
function(formula,condition,lc) {
    # tabulate formula given condition, based on lc model
    y <- lc$y
    mf <- model.frame(formula,y)
    ret <- NULL
    trap <- FALSE
    if (any(condition<=0) | any(condition>apply(y[names(condition)],2,max,na.rm=T))) {
        cat("Error: Some 'condition' values are not observed in data set. \n")
        trap <- TRUE
    } else if (any(table(c(names(condition),names(mf)))>1)) {
        cat("Error: Variables can only be specified once in 'formula' or 'condition'. \n")
        trap <- TRUE
    } else if (ncol(mf)>2) {
        cat("Error: 'formula' must be of form 'y~1' or 'y~x'. \n")
        trap <- TRUE
    }
    if (!trap) {
        grp <- F
        sel <- list()
        for (j in 1:ncol(y)) {
            if (names(y)[j] %in% names(mf)) {
                sel[[j]] <- c(1:max(mf[,which(names(mf) == names(y)[j])]))
            } else {
                if (sum(names(condition) == names(y)[j])==0) {
                    sel[[j]] <- c(1:max(y[,j],na.rm=T))
                    grp <- TRUE
                } else {
                    sel[[j]] <- condition[[which(names(condition) == names(y)[j])]]
                }
            }
        }
        names(sel) <- names(y)
        yc <- expand.grid(sel)
        predcell <- lc$N * (poLCA.ylik.C(poLCA.vectorize(lc$probs),yc) %*% lc$P)
        if (ncol(mf)>1) {
            ord <- 1+(which(names(mf)[1]==names(y)) > which(names(mf)[2]==names(y)))
            if (grp) {
                pc.col <- NULL
                for (i1 in 1:max(mf[,3-ord])) {
                    for (i2 in 1:max(mf[,ord])) {
                        pc.col <- c(pc.col,sum(predcell[yc[,which(names(y) %in% names(mf)[3-ord])]==i1 &
                                                        yc[,which(names(y) %in% names(mf)[ord])]==i2]))
                    }
                }
                predcell <- pc.col
            }
            nr <- apply(mf,2,max)[ord]
            nc <- apply(mf,2,max)[3-ord]
            ret <- matrix(predcell,nrow=nr,ncol=nc)
            rownames(ret) <- paste(names(mf)[ord],c(1:max(mf[,ord])))
            colnames(ret) <- paste(names(mf)[3-ord],c(1:max(mf[,3-ord])))
            if (ord==2) { ret <- t(ret) }
        } else {
            if (grp) {
                pc.col <- NULL
                for (i in 1:max(mf)) {
                    pc.col <- c(pc.col,sum(predcell[yc[,which(names(y) %in% names(mf))]==i]))
                }
                predcell <- pc.col
            }
            ret <- matrix(predcell,nrow=1,ncol=max(mf))
            rownames(ret) <- ""
            colnames(ret) <- paste(names(mf)[1],c(1:max(mf)))
        }
        return(ret)
    } else {
        invisible(NULL)
    }
}


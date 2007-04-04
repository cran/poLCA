`poLCA.ylik.C` <-
function(vp,y) {
    res <- .C("ylik",
                as.double(vp$vecprobs),
                as.integer(t(y)),
                as.integer(dim(y)[1]),
                as.integer(length(vp$numChoices)),
                as.integer(vp$numChoices),
                as.integer(vp$classes),
                lik = double(dim(y)[1]*vp$classes)
            )
    res$lik <- matrix(res$lik,ncol=vp$classes,byrow=TRUE)
    return(res$lik)
}


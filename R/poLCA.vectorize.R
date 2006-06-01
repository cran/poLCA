"poLCA.vectorize" <-
function(probs) {
    vecprobs <- as.numeric(t(probs[[1]]))
    numChoices <- dim(probs[[1]])[2]
    classes <- dim(probs[[1]])[1]             
    for (i in 2:length(probs)) {
        vecprobs <- append(vecprobs,as.numeric(t(probs[[i]])))
        numChoices <- append(numChoices,dim(probs[[i]])[2])         
    } 
    return(list(vecprobs=vecprobs,numChoices=numChoices,classes=classes))
}


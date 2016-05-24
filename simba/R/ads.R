"ads" <-
function(nspec, nplots, avSR=NULL, anc=NULL, grad.v=NULL, reord=TRUE, cf=0.2, puq=0.01) {
    if (!is.null(anc)){
        nspec <- ncol(anc)
        nplots <- nrow(anc)
        if (is.null(avSR)){
            anc <- ifelse(anc > 0, 1, 0)
            avSR <- mean(rowSums(anc))
        }
    }
    prob.v <- 1/(c(2:(nspec+1))^cf)
    prob.v <- prob.v - min(prob.v) + 1/nplots
    prob.v <- prob.v/max(prob.v)
    prob.mat <- sapply(prob.v, function(x) c(rep(1, round(nplots*x, 0)), rep(0, (nplots- round(nplots*x, 0)))))
        
    if (!is.null(grad.v)){
        if (!is.null(puq)){
            nuq <- round(nspec*puq, 0)
            muq <- replicate(nuq, rep(1, nplots))
        }
        else {
            muq <- NULL
            nuq <- 0
        }
    if(length(grad.v)!=nplots){
        stop("gradient vector must have a length of nplots")
    }    
    grad.max <- max(grad.v)
    grad.mat1 <- sapply(grad.v, function(x) sample(c(rep(1, x), rep(0, grad.max-x))))
    grad.mat2 <- sapply(abs(grad.v-max(grad.v)), function(x) sample(c(rep(1, x), rep(0, grad.max-x))))
    divid <- runif(nspec-nuq)
    divid1 <- c(1:sum(divid <= 0.5))
    divid2 <- c(1:(length(divid) - length(divid1)))
    mgrad1 <- sapply(divid1, function(x) apply(grad.mat1, 2, "sample", 1))
    mgrad2 <- sapply(divid2, function(x) apply(grad.mat2, 2, "sample", 1))
    mat <- cbind(muq, mgrad1, mgrad2)
    if(reord){
        mat <- mat[,c(sample(c(1:ncol(mat))))]
        }
    nspec <- ncol(mat)
    dis.v <- colSums(prob.mat)
    corr <- colSums(mat) < colSums(prob.mat)
    for (i in 1:sum(corr)){
        mat[,corr][,i] <- sample(c(rep(1, dis.v[corr][i]), rep(0,     nplots-dis.v[corr][i])))
        }
    for (i in 1:nspec){
        mat[mat[,i]==1,i] <- sample(c(rep(1, dis.v[i]), rep(0, length(mat[mat[,i]==1,i])-dis.v[i])))
        }
    }
    else {
    mat <- apply(prob.mat, 2, sample)
    }
    nspec <- ncol(mat)
    mat <- data.frame(mat)
    rownames(mat) <- c(1:nrow(mat))
    names(mat) <- c(1:ncol(mat))
    attr(mat, "lengths") <- nspec
    return(mat)
}

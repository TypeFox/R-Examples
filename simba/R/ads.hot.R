"ads.hot" <-
function(nspec, nplots, avSR=NULL, anc=NULL, grad.v=NULL, frac=0.5, reord=TRUE, cf=0.2, puq=0.01) {
    if (!is.null(anc)){
        nspec <- ncol(anc)
        nplots <- nrow(anc)
        if (is.null(avSR)){
            anc <- ifelse(anc > 0, 1, 0)
            avSR <- mean(rowSums(anc))
        }
    }
    if(length(grad.v)!=nplots){
        stop("gradient vector must have a length of nplots")
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
    grad.max <- max(grad.v)
    grad.mat <- sapply(grad.v, function(x) sample(c(rep(1, x), rep(0, grad.max-x))))
    divid1 <- c(1:(frac*(nspec-nuq)))
    divid2 <- c(1:(nspec-nuq-length(divid1)))
    mgrad <- sapply(divid1, function(x) apply(grad.mat, 2, "sample", 1))
    mrand <- sapply(divid2, function(x) ifelse(runif(nplots) > 0.5, 1, 0))
    mat <- cbind(muq, mgrad, mrand)
    if(reord){
        mat <- mat[,c(sample(c(1:ncol(mat))))]
        }
#    applying natural species occurrences onto the places were species could occur regarding to random distribution and hotspot gradient
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

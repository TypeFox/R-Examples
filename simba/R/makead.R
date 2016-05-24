"makead" <-
function(nspec, nplots, avSR=NULL, anc=NULL, grad.v=NULL, cf=0.2, puq=0.01) {
#    grad.v is a vector describing the gradient setting, if it is set a gradient is applied.
#    puq gives the proportion of ubiquituous species which are allowed to grow everywhere
#    if ancestor is given, nspec, nplots and avSR is generated from this matrix, but avSR can still be changed
    if (!is.null(anc)){
        nspec <- ncol(anc)
        nplots <- nrow(anc)
        if (is.null(avSR)){
            anc <- ifelse(anc > 0, 1, 0)
            avSR <- mean(rowSums(anc))
        }
    }
    probv <- 1/c(2:(nspec+1))^cf
    probv <- probv - min(probv) + 1/nplots
    probv <- probv/max(probv)
    mat <- matrix(NA, nplots, nspec)
    for(i in 1:nspec) {
        mat[,i] <- as.vector(rbind(matrix(1, round(nplots*probv[i],0), 1), matrix(0, (nplots-round(nplots*probv[i],0)), 1)))
    }
    mat <- apply(mat, 2, "sample")
    if (!is.null(grad.v)){
        if (!is.null(puq)){
            nuq <- round(nspec*puq, 0)
            vuq <- rep(1, nuq)
            vuq <- sample(c(vuq,vuq-1))
            vuq2 <- rep(0, ncol(mat)-length(vuq))
            vuq <- c(vuq, vuq2)
            vuqs <- vuq == 1
            vrest <- vuq == 0
            muq <- mat[,vuqs]
            muq <- muq + matrix(rnorm(ncol(muq)*nrow(muq), mean=0.5, sd=0.2), nrow(muq), ncol(muq))
            mgrad <- mat[,vrest]
        }
        else {
            mgrad <- mat
            muq <- NULL
        }
    divid <- runif(ncol(mgrad))
    divid1 <- ifelse(divid <= 0.5, TRUE, FALSE)
    divid2 <- ifelse(divid <= 0.5, FALSE, TRUE)
    grad.v <- as.numeric(grad.v)
    vgrad1 <- grad.v/max(grad.v)
    vgrad2 <- 1-vgrad1
    mgrad1 <- mgrad[,divid1]*vgrad1
    mgrad2 <- mgrad[,divid2]*vgrad2
    mrand <- matrix(rnorm(ncol(mgrad)*nrow(mgrad), mean=0.5, sd=0.2), nrow(mgrad), ncol(mgrad))
    mrand <- ifelse(mrand<0, 0, mrand)
    mgrad1 <- mgrad1 + mrand[,divid1]
    mgrad2 <- mgrad2 + mrand[,divid2]
    mat <- cbind(muq, mgrad1, mgrad2)
    contr <- nspec*nplots
    x.f <- 1                                                  
    for (i in 1:100){
      mat.tmp <- ifelse(mat>=(i/100), 1, 0)
      x.f <- ifelse(abs(sum(mat.tmp)-avSR*nplots)<=contr,i/100, x.f)  # keep optimum value
      contr <- ifelse(abs(sum(mat.tmp)-avSR*nplots)<=contr, abs(sum(mat.tmp)-avSR*nplots), contr) 
      }
    mat <- ifelse(mat>=x.f, 1, 0)
    ##save rare species
    ##mat[1,colSums(mat)==0] <- 1
    nspec <- c(ncol(muq), ncol(mgrad1), ncol(mgrad2))
    }
    mat <- mat[,c(sample(c(1:ncol(mat))))]
    mat <- data.frame(mat, row.names=as.character(c(1:nrow(mat))))
    names(mat) <- as.character(c(1:ncol(mat)))
    attr(mat, "lengths") <- nspec
    return(mat)
}
###############################################################################
## Computation of distances for RFLP data - alternative approach
###############################################################################

## compute number of bands
nrBands <- function(x){
    sort(unique(sapply(split(x, x$Sample), nrow)))
}

## x: data.frame with RFLP data
## distfun: function to compute distance (e.g. ?dist)
## nrBands: number of bands 
## nrMissing: number of bands which may be nrMissing
## diag: see ?dist
## upper: see ?dist
## compares samples with number of bands in: nrBands, nrBands + 1, ..., nrBands + nrMissing
RFLPdist2 <- function(x, distfun = dist, nrBands, nrMissing, LOD = 0, diag = FALSE, upper = FALSE){
    stopifnot(is.data.frame(x))
    stopifnot(is.function(distfun))
    if(missing(nrMissing))
        stop("'nrMissing' is not specified!")
    if(nrMissing == 0)
        stop("'nrMissing == 0', please use function 'RFLPdist'!")
    if(nrMissing <= 0)
        stop("'nrMissing' has to be a positive interger!")    
    if(missing(nrBands))
        stop("'nrBands' is not specified!")
    if(nrBands <= 0)
        stop("'nrBands' has to be a positive interger!")
    if(LOD < 0)
        stop("'LOD' has to be non-negative!")

    if(LOD == 0){
        x1 <- split(x, x$Sample)
        x1.bands <- sapply(x1, nrow)

        temp <- do.call("rbind", x1[x1.bands %in% c(nrBands:(nrBands+nrMissing))])
        temp1 <- split(temp[,"MW"], factor(temp[,"Sample"]))
        N <- length(temp1)
        temp2 <- matrix(NA, ncol = nrBands + nrMissing, nrow = N)
        rownames(temp2) <- names(temp1)
        for(i in 1:N){
            temp2[i,1:length(temp1[[i]])] <- temp1[[i]]
        }
        d <- matrix(NA, nrow = N, ncol = N)
        dfun1 <- function(x, y){
            y <- y[!is.na(y)]
            m <- sum(!is.na(x))
            min(as.matrix(distfun(rbind(x[1:m], t(combn(y, m)))))[-1,1])
        }
        for(i in 1:N){
            for(j in 1:i){
                if(sum(!is.na(temp2[i,])) == sum(!is.na(temp2[j,]))){
                    m <- sum(!is.na(temp2[i,]))
                    d[i,j] <- as.vector(distfun(rbind(temp2[i,1:m], temp2[j,1:m])))
                }
                if(sum(!is.na(temp2[i,])) > sum(!is.na(temp2[j,])))
                    d[i,j] <- dfun1(temp2[j,], temp2[i,])
                if(sum(!is.na(temp2[i,])) < sum(!is.na(temp2[j,])))
                    d[i,j] <- dfun1(temp2[i,], temp2[j,])
            }
        }
    }else{
        if(nrMissing == 0)
            stop("'nrMissing' has to be at least 1")
        ## consider bands >= LOD
        x1 <- x[x$MW >= LOD,]
        x2 <- split(x1, x1$Sample)
        x2.bands <- sapply(x2, nrow)
        
        ## remove samples with no band >= LOD
        x3 <- split(x, x$Sample)
        x4 <- x3[names(x3) %in% names(x2)]
        
        ## consider only samples where number of bands with MW >= LOD
        ## is equal to nrBands (missing bands only for MW < LOD)
        x5 <- do.call("rbind", x4[x2.bands == nrBands])
        x6 <- split(x5, x5$Sample)
        x6.bands <- sapply(x6, nrow)

        ## extract samples where up to nrMissing number of bands occur
        ## for MW < LOD
        temp <- do.call("rbind", x6[x6.bands %in% c(nrBands:(nrBands+nrMissing))])
        temp1 <- split(temp[,"MW"], factor(temp[,"Sample"]))
        N <- length(temp1)
        temp2 <- matrix(NA, ncol = nrBands + nrMissing, nrow = N)
        rownames(temp2) <- names(temp1)
        for(i in 1:N){
            temp2[i,1:length(temp1[[i]])] <- temp1[[i]]
        }
        dfun2 <- function(x, y, LOD){
            y <- y[!is.na(y)]
            m <- sum(!is.na(x))
            m.lod <- sum(!is.na(x[x < LOD]))
            y1 <- t(combn(y[y < LOD], m.lod))
            y2 <- matrix(rep(y[y >= LOD], nrow(y1)), nrow = nrow(y1), byrow = TRUE)
            y3 <- cbind(y2, y1)
            min(as.matrix(distfun(rbind(x[1:m], y3)))[-1,1])
        }
        d <- matrix(NA, nrow = N, ncol = N)
        for(i in 1:N){
            for(j in 1:i){
                if(sum(!is.na(temp2[i,])) == sum(!is.na(temp2[j,]))){
                    m <- sum(!is.na(temp2[i,]))
                    d[i,j] <- as.vector(distfun(rbind(temp2[i,1:m], temp2[j,1:m])))
                }
                if(sum(!is.na(temp2[i,])) > sum(!is.na(temp2[j,])))
                    d[i,j] <- dfun2(temp2[j,], temp2[i,], LOD = LOD)
                if(sum(!is.na(temp2[i,])) < sum(!is.na(temp2[j,])))
                    d[i,j] <- dfun2(temp2[i,], temp2[j,], LOD = LOD)
            }
        }        
    }
    d <- d[lower.tri(d)]
    attr(d, "Size") <- N
    attr(d, "Labels") <- dimnames(temp2)[[1]]
    attr(d, "Diag") <- diag
    attr(d, "Upper") <- upper
    attr(d, "method") <- distfun
    attr(d, "call") <- match.call()
    class(d) <- "dist"
    return(d)
}

`rankabundance` <-
function(x,y="",factor="",level,digits=1,t=qt(0.975,df=n-1)) {
    if(class(y) == "data.frame" && factor != "") {
        subs <- y[,factor]==level
        for (q in 1:length(subs)) {
            if(is.na(subs[q])) {subs[q]<-F}
        }
        x <- x[subs,,drop=F]
        freq <- apply(x,2,sum)
        subs <- freq>0
        x <- x[,subs,drop=F]
    }
    if(dim(as.matrix(x))[1]==0) {
        result <- array(NA,dim=c(1,8))
        colnames(result) <- c("rank","abundance","proportion","plower","pupper","accumfreq","logabun","rankfreq")
        rownames(result) <- "none"
        return(result)
    }    
    total <- apply(x,1,sum)    
    p <- ncol(x)
    n <- nrow(x)
    mu <- sum(total)/n
    result <- array(dim=c(p,8))
    colnames(result) <- c("rank","abundance","proportion","plower","pupper","accumfreq","logabun","rankfreq")
    rownames(result) <- colnames(x)
    for (j in 1:p) {
        spec <- x[,j]
        pi <- spec/total
        p <- sum(spec)/sum(total)
        sigma2 <- 0
        for (i in 1:n) {
            sigma2 <- sigma2 + (total[i]^2 * (pi[i]-p)^2)
        }
        sigma2 <- sigma2 / (n * (n-1) * mu * mu)
        sigma <- sigma2^0.5
        result[j,2] <- sum(spec)
        result[j,3] <- p*100
        result[j,4] <- (p - t*sigma)*100
        result[j,5] <- (p + t*sigma)*100
    }
    p <- ncol(x)
    result2 <- result
    seq <- rev(order(result[,2],-order(rownames(result))))
    result[1:p,] <- result2[seq,]
    rownames(result)[1:p] <- rownames(result2)[seq]
    result[,1] <- c(1:ncol(x))
    result[,6] <- cumsum(result[,3])
    result[,7] <- log(result[,2],base=10)
    result[,8] <- result[,1]/ncol(x)*100
    result <- round(result,digits=digits)
    return(result)
}


`accumresult` <-
function(x,y="",factor="",level,scale="",method="exact",permutations=100,conditioned=T,gamma="boot",...){
    op <- options()
    options(warn=-1)
    subs <- c(1:nrow(x))
    if(class(y) == "data.frame" && factor != "") {
        subs <- y[,factor]==level
        for (q in 1:length(subs)) {
            if(is.na(subs[q])) {subs[q]<-F}
        }
        x <- x[subs,,drop=F]
        freq <- apply(x,2,sum)
        subs2 <- freq>0
        x <- x[,subs2,drop=F]
    }
    if(dim(as.matrix(x))[1]==0) {
        result <- list(call = match.call(), method = method, sites = 0, richness = NA, sd = NA, perm = NA)
        return(result)
    }
    result <- specaccum(x,method=method,permutations=permutations,conditioned=conditioned,gamma=gamma,...)
    if (scale != "") {    
        y <- y[subs,,drop=F]
        tot <- mean(y[,scale])
        result$sites <- result$sites * tot
    }
    options(op)
    return(result)
}


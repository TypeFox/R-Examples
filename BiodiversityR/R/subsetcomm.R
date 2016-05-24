`subsetcomm` <-
function(x, y, factor, level, returncomm=TRUE) {
    subs <- y[,factor]==level
    for (q in 1:length(subs)) {
        if(is.na(subs[q])) {subs[q]<-F}
    }
    if (returncomm==T) {
        x <- x[subs,,drop=F]
        freq <- apply(x,2,sum)
        subs <- freq>0
        x <- x[,subs,drop=F]
        return(x)
    }else{
        y <- y[subs,,drop=F]
        for (i in 1:ncol(y)) {
            if (is.factor(y[,i])) {y[,i] <- factor(y[,i][drop=T])}
        }
        return(y)
    }
}

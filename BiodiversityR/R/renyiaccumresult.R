`renyiaccumresult` <-
function(x,y="",factor,level,scales=c(0,0.25,0.5,1,2,4,8,Inf),permutations=100,...) {
    if(class(y) == "data.frame") {
        subs <- y[,factor]==level
        for (q in 1:length(subs)) {
            if(is.na(subs[q])) {subs[q]<-F}
        }
        x <- x[subs,,drop=F]
        freq <- apply(x,2,sum)
        subs <- freq>0
        x <- x[,subs,drop=F]
    }
    result <- renyiaccum(x,scales=scales,permutations=permutations,...)
    return(result)
}



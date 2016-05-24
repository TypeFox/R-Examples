`radfitresult` <-
function(x,y="",factor,level,plotit=T){
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
    x <- as.matrix(apply(x,2,sum))
    result1 <- radfit(x)
    result2 <- fisherfit(x)
    result3 <- prestonfit(x)
    if(plotit==T) {graphics::plot(result1)}
    result <- list(radfit=result1, fisherfit=result2, prestonfit=result3)
    return(result)
}


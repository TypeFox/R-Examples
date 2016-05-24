divchain.deldir <- function (x,...) {
#
    z <- x$summary$z 
    if(!is.factor(z)) {
        xc <- deparse(substitute(x))
        whinge <- paste("The class deldir object",xc,"was created without\n",
                        "a factor-valued \"weights\" argument \"z\" being supplied.\n")
        stop(whinge)
    }
    ddd  <- x$dirsgs
    ddd  <- ddd[z[ddd$ind1] != z[ddd$ind2],]
    id1  <- as.matrix(ddd[,c("ind1","ind2","thirdv1")])
    id2  <- as.matrix(ddd[,c("ind1","ind2","thirdv2")])
    id1  <- t(apply(id1,1,function(x){if(x[3] > 0) sort(x) else c(sort(x[1:2]),x[3])}))
    id2  <- t(apply(id2,1,function(x){if(x[3] > 0) sort(x) else c(sort(x[1:2]),x[3])}))
    rslt <- cbind(ddd[,1:4],id1,id2)
    names(rslt) <- c("x0","y0","x1","y1","v01","v02","v03","v11","v12","v13")
    class(rslt) <- c("divchain","data.frame")
    attr(rslt,"rw") <- x$rw
    rslt
}

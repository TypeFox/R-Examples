"liste" <-
function(obj, x="NBX", y="NBY", entry=NULL, factorize=TRUE, splist=FALSE){
    if (class(obj)=="dist"){
        obj.v <- as.vector(obj)
        if(is.null(attr(obj, "Labels"))) attr(obj, "Labels")  <- as.character(c(1:attr(obj, "Size")))
        name.mat <- outer(attr(obj, "Labels"), attr(obj, "Labels"), "paste")
        name.list <- name.mat[row(name.mat) > col(name.mat)]
    }
    else {
        obj.v <- as.vector(as.matrix(obj))
        obj <- data.frame(obj)
        name.mat <- outer(rownames(obj), names(obj),  "paste")
        name.list <- as.vector(name.mat)
    }
    list.split <- strsplit(name.list, " ")
    NBY <- sapply(list.split, function(x) x[1])
    NBX <- sapply(list.split, function(x) x[2])
    dat.lst <- data.frame(NBX, NBY, obj.v)
    if (factorize){
    	dat.lst[,1] <- as.factor(dat.lst[,1])
    	dat.lst[,2] <- as.factor(dat.lst[,2])
    	}
    if (is.null(entry)){
        entry <- "we"
    }
    names(dat.lst) <- c(x, y, entry)
    if(splist) {
        dat.lst <- dat.lst[(dat.lst[,3]!=0),c(2,1,3)]
        dat.lst <- dat.lst[order(dat.lst[,1]),]
        rownames(dat.lst) <- c(1:nrow(dat.lst))
        names(dat.lst) <- c("plot", "spec", "occ")
    }
    out <- dat.lst
    return(out)
}
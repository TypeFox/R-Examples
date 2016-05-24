`maxpear` <-
function(psm, cls.draw=NULL, method=c("avg","comp","draws","all"), max.k=NULL){
    
    if(any(psm !=t(psm)) | any(psm >1) | any(psm < 0) | sum(diag(psm)) != nrow(psm) ){
     stop("psm must be a symmetric matrix with entries between 0 and 1 and 1's on the diagonals")}
    
    method <- match.arg(method, choices=method)
    if(method %in% c("draws","all") & is.null(cls.draw)) stop("cls.draw must be provided if method=''draws''")
    
    if(method == "avg" | method == "all"){
        if(is.null(max.k)) max.k <- ceiling(dim(psm)[1]/8)
        hclust.avg <- hclust(as.dist(1-psm), method="average")
        cls.avg <-  t(apply(matrix(1:max.k),1,function(x) cutree(hclust.avg,k=x)))
        pears.avg <- pear(cls.avg,psm)
        val.avg <- max(pears.avg)
        cl.avg <- cls.avg[which.max(pears.avg),] 
        if(method== "avg") return(list(cl=cl.avg, value=val.avg, method="avg"))
    }
    
    if(method == "comp" | method == "all"){
        if(is.null(max.k)) max.k <- ceiling(dim(psm)[1]/8)
        hclust.comp <- hclust(as.dist(1-psm), method="complete")
        cls.comp <-  t(apply(matrix(1:max.k),1,function(x) cutree(hclust.comp,k=x)))
        pears.comp <- pear(cls.comp,psm)
        val.comp <- max(pears.comp)
        cl.comp <- cls.comp[which.max(pears.comp),] 
        if(method== "comp") return(list(cl=cl.comp, value=val.comp, method="comp"))
    }
    
    if(method == "draws" | method == "all"){
        compIpi <- function(cl){
            mat <- cltoSim(cl)*psm
            sum(mat[lower.tri(mat)])
            }
        n <- ncol(psm)
        no2 <- choose(n,2)
        sumpij <- sum(psm[lower.tri(psm)])
        sumIij <- apply(cls.draw,1, function(x) sum(choose(table(x),2)))
        sumIpi <- apply(cls.draw,1, compIpi)
        correc <- (sumIij*sumpij)/no2
        pears.draws <- (sumIpi - correc)/ (0.5*(sumpij+sumIij)-correc)
        val.draws <- max(pears.draws)
        cl.draw <- cls.draw[which.max(pears.draws),] 
        names(cl.draw) <- NULL
        if(method== "draws") return(list(cl=cl.draw, value=val.draws, method="draws"))
    }
    
    vals <- c(val.avg, val.comp, val.draws)
    cls <- rbind(cl.avg,cl.comp,cl.draw)
    cls <- rbind(cls[which.max(vals),], cls)
    vals <- c(max(vals), vals)
    rownames(cls) <- names(vals) <- c("best","avg","comp","draws")
    colnames(cls) <- NULL
    list(cl=cls, value=vals, method="all")
}

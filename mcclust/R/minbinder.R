`minbinder` <-
function(psm, cls.draw=NULL, method=c("avg","comp","draws","laugreen","all"), max.k=NULL, include.lg=FALSE, start.cl=NULL, tol=0.001){
    
    if(any(psm !=t(psm)) | any(psm >1) | any(psm < 0) | sum(diag(psm)) != nrow(psm) ){
     stop("psm must be a symmetric matrix with entries between 0 and 1 and 1's on the diagonals")}
    
    method <- match.arg(method, choices=method)
    if(method %in% c("draws","all") & is.null(cls.draw)) stop("cls.draw must be provided if method=''draws''")
    
    if(method == "avg" | method == "all"){
        if(is.null(max.k)) max.k <- ceiling(dim(psm)[1]/4)
        hclust.avg <- hclust(as.dist(1-psm), method="average")
        cls.avg <-  t(apply(matrix(1:max.k),1,function(x) cutree(hclust.avg,k=x)))
        binder.avg <- binder(cls.avg,psm)
        val.avg <- min(binder.avg)
        cl.avg <- cls.avg[which.min(binder.avg),] 
        if(method== "avg") return(list(cl=cl.avg, value=val.avg, method="avg"))
    }
    
    if(method == "comp" | method == "all"){
        if(is.null(max.k)) max.k <- ceiling(dim(psm)[1]/4)
        hclust.comp <- hclust(as.dist(1-psm), method="complete")
        cls.comp <-  t(apply(matrix(1:max.k),1,function(x) cutree(hclust.comp,k=x)))
        binder.comp <- binder(cls.comp,psm)
        val.comp <- min(binder.comp)
        cl.comp <- cls.comp[which.min(binder.comp),] 
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
        binder.draws <- sumIij - 2*sumIpi + sumpij 
        val.draws <- min(binder.draws)
        cl.draw <-  cls.draw[which.min(binder.draws),] 
        names(cl.draw) <- NULL
        if(method== "draws") return(list(cl=cl.draw, value=val.draws, method="draws"))
    }
    
    if(method == "laugreen" | (method == "all" & include.lg)){
        res.lg <- laugreen(psm, start.cl=start.cl, tol=tol)
        if(method=="laugreen") return(res.lg)  
    }
    
    vals <- c(val.avg, val.comp, val.draws)
    cls <- rbind(cl.avg,cl.comp,cl.draw)
    if(include.lg){
            vals <- c(vals,res.lg$value)
            cls <- rbind(cls,res.lg$cl)
    }
    cls <- rbind(cls[which.min(vals),], cls)
    vals <- c(min(vals), vals)
    if(include.lg){ rownames(cls) <- names(vals) <- c("best","avg","comp","draws","laugreen")
    } else rownames(cls) <- names(vals) <- c("best","avg","comp","draws")
        colnames(cls) <- NULL    
        res <- list(cl=cls, value=vals)
    if(include.lg) return(c(res,list(iter.lg=res.lg$iter.lg)))
    res
}

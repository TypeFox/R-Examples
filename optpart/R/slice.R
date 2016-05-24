slice <- function (clust, k=NULL) 
{
    if (is.null(k)) {
        tmp <- locator(n=1)
        abline(tmp$y,0,col=2)
        tmp <- list(clustering=cutree(clust,h=tmp$y))
        cat(paste("Number of clusters = ",max(tmp$clustering),"\n"))
    }
    else {
        tmp <- list(clustering=cutree(clust,k=k))
    }
    attr(tmp,"class") <- "clustering"
    invisible(tmp)
}


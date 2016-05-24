mergeclust <- function(clustering,from,to) 
{
    clustering <- as.integer(clustify(clustering))

    clustering[clustering==from] <- to
    out <- list()
    out$clustering <- as.numeric(factor(clustering))
    out
}

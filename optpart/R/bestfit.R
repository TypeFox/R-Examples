bestfit <- function(x,cluster)
{
    if (inherits(x,'partana')) {
       if (cluster < 1 || cluster > ncol(x$ptc))
            stop('Cluster must be > 1 and <= number of clusters')
       clustering <- as.numeric(clustify(x))
       count <- sum(clustering==cluster)
       if (count > 1) {
           rows <- x$ptc[clustering==cluster,]
           vals <- rows[,cluster]
           names <- x$names[clustering==cluster]
           fit <- data.frame(names[rev(order(vals))],sort(vals,decreasing=TRUE))
           names(fit) <- c('ID','fit')
           return(fit)
       } else {
           print(paste("Cluster ",cluster," has too few members"))
       }
    } else if (inherits(x,'silhouette')) {
       if (cluster < 1 || cluster > max(x[,1]))
            stop('Cluster must be > 1 and <= number of clusters')
       rows <- x[,1] == cluster
       names <- seq(1,nrow(x))
       count <- length(rows)
       if (count > 1) {
           vals <- as.numeric(x[rows,3])
           names <- names[rows]
           fit <- data.frame(names[rev(order(vals))],sort(vals,decreasing=TRUE))
           names(fit) <- c('ID','fit')
           return(fit)
       } else {
           print(paste("Cluster ",cluster," has too few members"))
       }
    } else {
        stop('You must pass an object of class partana or silhouette')
    }
}


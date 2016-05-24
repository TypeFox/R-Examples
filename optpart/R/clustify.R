clustify <- function(clustering)
{
    if (inherits(clustering,c("partana", "partition", "clustering"))) {
        clustering <- factor(clustering$clustering)
    } else if (is.character(clustering)) {
        clustering <- factor(clustering)
    } else if (is.numeric(clustering)) {
        clustering <- factor(clustering)
    } else if (is.logical(clustering)) {
        clustering <- factor(clustering)
    } else if (!is.factor(clustering)) 
        stop('Cannot understand passed clustering')
    clustering
}

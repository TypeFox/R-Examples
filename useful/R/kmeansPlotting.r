# k-means plotting

#' fortify.kmeans
#' 
#' Fortify a kmeans model with its data
#' 
#' Prepares a kmeans object to be plotted using \code{\link{cmdscale}} to compute the projected x/y coordinates.  If \code{data} is not provided, then just the center points are calculated.
#' 
#' @aliases fortify.kmeans
#' @export
#' @export fortify.kmeans
#' @importFrom stats cmdscale dist
#' @author Jared P. Lander
#' @seealso kmeans fortify ggplot plot.kmeans
#' @param model \code{\link{kmeans}} model
#' @param data Data used to fit the model
#' @param \dots Not Used
#' @return The original data with extra columns:
#'      \item{.x}{The projected x position.}
#'      \item{.y}{The projected y position.}
#'      \item{.Cluster}{The cluster that point belongs to.}
#' @examples
#' 
#' k1 <- kmeans(x=iris[, 1:4], centers=3)
#' hold <- fortify(k1, data=iris)
#' head(hold)
#' hold2 <- fortify(k1)
#' head(hold2)
#' 
fortify.kmeans <- function(model, data=NULL, ...)
{
    # get the names of columns used
    usedCols <- colnames(model$centers)
    
    if(is.null(data))
    {
        # get 2 dimensional scaling of the centers
        centerPoints <- data.frame(cmdscale(d=dist(model$centers), k=2))
        names(centerPoints) <- c(".x", ".y")
        centerPoints$.Cluster <- as.factor(rownames(centerPoints))
        
        return(centerPoints)
    }
    
    # make a 2 dimensional scaling of the data
    points <- data.frame(cmdscale(d=dist(data[, usedCols]), k=2))
    names(points) <- c(".x", ".y")
    
    # tack centers onto the points
    points$.Cluster <- as.factor(model$cluster)
    
    data <- cbind(data, points)
    
    return(data)
}

#' plot.kmeans
#' 
#' Plot the results from a k-means object
#' 
#' Plots the results of k-means with color-coding for the cluster membership.  If \code{data} is not provided, then just the center points are calculated.
#' 
#' @aliases plot.kmeans
#' @export
#' @export plot.kmeans
#' @author Jared P. Lander
#' @seealso kmeans fortify ggplot plot.kmeans
#' @param x A \code{\link{kmeans}} object.
#' @param data The data used to kit the \code{\link{kmeans}} object.
#' @param class Character name of the "true" classes of the data.
#' @param legend.position Character indicating where the legend should be placed.
#' @param title Title for the plot.
#' @param xlab Label for the x-axis.
#' @param ylab Label for the y-axis.
#' @param \dots Not Used.
#' @return A ggplot object
#' @examples
#' 
#' k1 <- kmeans(x=iris[, 1:4], centers=3)
#' plot(k1)
#' plot(k1, data=iris)
#' 
plot.kmeans <- function(x, data=NULL, class=NULL, legend.position=c("right", "bottom", "left", "top", "none"), 
                        title="K-Means Results",
                        xlab="Principal Component 1", ylab="Principal Component 2", ...)
{
    # fortify the model and data so it is convenient to plot in ggplot
    toPlot <- fortify(model=x, data=data)
    
    # get the legend position
    legend.position <- match.arg(legend.position)
    
    # convert class to factor just in case it is not already
    if(!is.null(class)) toPlot[, class] <- factor(toPlot[, class])
    
    ggplot(toPlot, aes_string(x=".x", y=".y", colour=".Cluster")) + 
        geom_point(aes_string(shape=class)) + 
        scale_color_discrete("Cluster") +
        theme(legend.position=legend.position) +
        labs(title=title, x=xlab, y=ylab)
}

## Functions for determining ideal number of clusters for kmeans

## Plots the results from the Hartigan's Rule run
## prints and returns the ggplot object
## @hartigan (data.frame) the results from fitting Hartigan's Rule
## @title (character) the title of the plot
## @linecolor (numeric) the color of the line indicating 10
## @linetype (numeric) the style of the line indicating 10
## @linesize (numeric) the size of the line indicating 10
## @minor (logical) whether the minor grid lines should be displayed



#' Plot a series of Hartigan's Numbers
#' 
#' After fitting a series of Hartigan's Numbers (see \code{\link{FitKMeans}}) this will plot the results so it is easy to visualize
#' 
#' Displays a graphical representation of the results of \code{\link{FitKMeans}}
#' 
#' @param hartigan The results from
#' \code{\link{FitKMeans}}
#' @param title Title to be used in the plot
#' @param smooth logical; if true a smoothed line will be fit to the points, otherwise it will be a piecewise line
#' @param linecolor Color of the horizontal line denoting 10
#' @param linetype Type of the horizontal line denoting 10
#' @param linesize Size of the horizontal line denoting 10
#' @param minor logical; if true minor grid
#' lines will be plotted
#' @return a ggplot object
#' @author Jared P. Lander
#' www.jaredlander.com
#' @import ggplot2
#' @export PlotHartigan
#' @seealso \code{\link{kmeans}} \code{\link{FitKMeans}}
#' @references #' http://www.stat.columbia.edu/~madigan/DM08/descriptive.ppt.pdf
#' @keywords cluster kmeans hartigan clustering
#' @examples
#' 
#' data(iris)
#' hartiganResults <- FitKMeans(iris[, -ncol(iris)])
#' PlotHartigan(hartiganResults)
#' 
PlotHartigan <- function(hartigan, title="Hartigan's Rule", smooth=FALSE, linecolor="grey", linetype=2L, linesize=1L, minor=TRUE)
{
    thePlot <- ggplot(data=hartigan, aes_string(x="Clusters", y="Hartigan")) + 
        geom_hline(aes(yintercept=10), linetype=linetype, colour=linecolor, size=linesize) + 
        #geom_line() +
        geom_point(aes_string(colour="AddCluster")) +
        scale_colour_discrete(name="Add Cluster") +
        ggtitle(label=title) + if(minor) scale_x_continuous(minor_breaks=(1:(max(hartigan$Clusters)+1)))
    
    if(smooth)
    {
        thePlot <- thePlot + geom_smooth(method="lm", formula=y~log(x))
    }else
    {
        thePlot <- thePlot + geom_line()
    }
    
    return(thePlot)
}


## Compute Hartigan's Rule given a kmeans cluster WSS and a k+1means cluster WSS and the number of rows in the data
## returns the number
#' Compute Hartigan's Number
#'
#' Runs the computation found in http://www.stat.columbia.edu/~madigan/DM08/descriptive.ppt.pdf
#'
#' Not exported, only used by \code{\link{FitKMeans}}
#'
#' @aliases ComputeHartigan
#' @param FitActualWSS the WSS from a kmeans fit
#' @param FitPlus1WSS the WSS from a kmeans fit
#' @param nrow the number of rows in the original dataset
#' @return The computed Hartigan Number
#' @references http://www.stat.columbia.edu/~madigan/DM08/descriptive.ppt.pdf
#' @author Jared P. Lander
#' www.jaredlander.com
#' @seealso \code{\link{kmeans}} \code{\link{FitKMeans}}
#' @keywords cluster kmeans hartigan clustering
#' @examples
#' data(iris)
#' hartiganResults <- FitKMeans(iris[, -ncol(iris)])
#' PlotHartigan(hartiganResults)
#'
ComputeHartigan <- function(FitActualWSS, FitPlus1WSS, nrow)
{
    return((sum(FitActualWSS) / sum(FitPlus1WSS) - 1) * (nrow - length(FitActualWSS) - 1))
}


## this function fits a series of kmeans and returns a data.frame listing the number of clusters and the result of applying Hartigan's Rule
## returns the data.frame of Hartigan results
## @x (data.frame or matrix) the data to fit kmeans on
## @max.clusters (numeric) the number of clusters to try
## @spectral (logical) whether it is fitting using spectral methods
## @nstart (numeric) the number of random starts for kmeans to use
## @iter.max (numeric) the maximum number of iterations for kmeans before giving up on convergence
## @seed (numeric) the random seed to be set



#' Fit a series of kmeans clusterings and compute Hartigan's Number
#' 
#' Given a numeric dataset this function fits a series of kmeans clusterings with increasing number of centers.  k-means is compared to k+1-means using Hartigan's Number to determine if the k+1st cluster should be added.
#' 
#' A consecutive series of kmeans is computed with increasing k (number of centers).  Each result for k and k+1 are compared using Hartigan's Number.  If the number is greater than 10, it is noted that having k+1 clusters is of value.
#' 
#' @param x The data, numeric, either a matrix or data.frame
#' @param max.clusters The maximum number of clusters that should be tried
#' @param spectral logical; If the data being fit are eigenvectors for spectral clustering
#' @param nstart The number of random starts for the kmeans algorithm to use
#' @param iter.max Maximum number of tries before the kmeans algorithm gives up on conversion
#' @param algorithm The desired algorithm to be used for kmeans.  Options are c("Hartigan-Wong", "Lloyd", "Forgy", "MacQueen").  See \code{\link{kmeans}}
#' @param seed If not null, the random seed will be reset before each application of the kmeans algorithm
#' @return A data.frame consisting of columns, for the number of clusters, the Hartigan Number and whether that cluster should be added, based on Hartigan's Number.
#' @author Jared P. Lander
#' www.jaredlander.com
#' @export FitKMeans
#' @importFrom stats kmeans
#' @seealso \code{\link{kmeans}} \code{\link{PlotHartigan}}
#' @references http://www.stat.columbia.edu/~madigan/DM08/descriptive.ppt.pdf
#' @keywords cluster kmeans hartigan clustering
#' @examples
#' 
#' data(iris)
#' hartiganResults <- FitKMeans(iris[, -ncol(iris)])
#' PlotHartigan(hartiganResults)
#' 
FitKMeans <- function(x, max.clusters=12L, spectral=FALSE, nstart=1L, iter.max=10L, algorithm=c("Hartigan-Wong", "Lloyd", "Forgy", "MacQueen"), seed=NULL)
{
    # get algorithm choice
    algorithm <- match.arg(algorithm)
    
	# data.frame for keeping track of Hartigan number
	hartigan <- data.frame(Clusters=2:(max.clusters), Hartigan=NA, AddCluster=NA)
 
 	# compute the number of rows and columns just once
 	nRowX <- nrow(x)
 	nColX <- ncol(x)
 	
    ## new algorithm
    # in each loop build one partition
    # compare to old partition
    # make new partition into old partition
    
    # first compute partition for 1 cluster
    if(!is.null(seed))
    {
        set.seed(seed=seed)
    }
    FitActual <- kmeans(x[, 1:(nColX - (nColX-(2-1))*spectral)], centers=2-1, nstart=nstart, iter.max=iter.max, algorithm=algorithm)
    
    ## now build loop
    for(i in 2:(max.clusters))
    {
        # calculate FitPlus1, which in this case will just be i
        if(!is.null(seed))
        {
            set.seed(seed=seed)
        }
        FitPlus1 <- kmeans(x[, 1:(nColX - (nColX-(i+0))*spectral)], centers=i, nstart=nstart, iter.max=iter.max, algorithm=algorithm)
        
        # calculate Hartigan and record in table
        hartigan[i-1, "Hartigan"] <- ComputeHartigan(FitActualWSS=FitActual$withinss, FitPlus1WSS=FitPlus1$withinss, nrow=nRowX)
        
        # now turn FitPlus1 into FitActual for use in the next iteration
        FitActual <- FitPlus1
        rm(FitPlus1); gc()          # housekeeping
    }
    
    ## could be made more efficient by fitting kmeans for each value of i one time, then comepare the consecutive pair wise results, would take mroe memory though
 	# compute kmeans repeatedly
#     for(i in 2:(max.clusters))
#     {
#         # for k
#         if(!is.null(seed))
#         {
#             set.seed(seed=seed)
#         }
#         FitActual <- kmeans(x[, 1:(nColX - (nColX-(i-1))*spectral)], centers=i-1, nstart=nstart, iter.max=iter.max, algorithm=algorithm)
#         
#         # for k+1
#         if(!is.null(seed))
#         {
#             set.seed(seed=seed)
#         }
#         FitPlus1 <- kmeans(x[, 1:(nColX - (nColX-(i+0))*spectral)], centers=i, nstart=nstart, iter.max=iter.max, algorithm=algorithm)
#         
#         # calculate Hartigan
#         hartigan[i-1, "Hartigan"] <- ComputeHartigan(FitActualWSS=FitActual$withinss, FitPlus1WSS=FitPlus1$withinss, nrow=nRowX)
#     }
 
    # if Hartigan is greater than 10 then the cluster should be added
    hartigan$AddCluster <- ifelse(hartigan$Hartigan > 10, TRUE, FALSE)
 
    return(hartigan)
}

# make compiled versions
# saving for a future version when compiler is more common
# ComputeHartigan <- cmpfun(ComputeHartigan)
# PlotHartigan <- cmpfun(PlotHartigan)
# FitKMeans <- cmpfun(FitKMeans)

#' Cluster detected events
#' 
#' This function groups the detected events into clusters. 
#' 
#' @details The clustering is based on statistical characteristics of
#' event. Each extracted event is first described using a feature vector, and then the events are clustered according
#' to the Euclidean distances among the feature vectors. Note that before clustering, we apply principal 
#' component analysis (PCA) to the feature vectors to reduce the correlation as well as the dimension of the feature space.

#' 
#' @param events an object of class `events'.
#' @param k0 the number of clusters.
#' @seealso \code{\link{measures}}
#' @references Xiaozhe Wang, Kate Smith-Miles and  Rob Hyndman (2005). Characteristic-Based Clustering for Time Series Data.
#' \emph{Data Mining and Knowledge Discovery}. \bold{13}(3), 335-364. \url{http://dx.doi.org//10.1007/s10618-005-0039-x}

#' @references Yanfei Kang, Danijel Belusic, Kate Smith-Miles (2014). Detecting and Classifying Events in Noisy Time Series. \emph{J. Atmos. Sci.}, \bold{71}, 1090-1104.
#' \url{http://dx.doi.org/10.1175/JAS-D-13-0182.1}.

#' @references Gregory S. Poulos, William Blumen, David C. Fritts, Julie K. Lundquist, Jielun Sun, Sean P. Burns, 
#' Carmen Nappo, Robert Banta, Rob Newsom, Joan Cuxart, Enric Terradellas, Ben Balsley, and Michael Jensen. 
#' CASES-99: A comprehensive investigation of the stable nocturnal boundary layer (2002). \emph{Bulletin of the American 
#' Meteorological Society}, \bold{83}(4):555-581. 

#' @return a list consisting of:
#' 
#' \item{cl}{a vector indicating which cluster each event belongs to.}
#' 
#' \item{center}{a matrix which gives cluster centroids.} 
#' 
#' \item{pca}{PCA results for characteristics of the detected events.}
#' 
#' @export
#' @examples
#' ##################################
#' #   An artificial example
#' ##################################
#' set.seed(123)
#' n <- 128
#' types <- c('box','rc','cr','sine')
#' shapes <- matrix(NA,20,n)
#' for (i in 1:20){
#'   shapes[i,] <- cbfs(type=types[sample(1:4,1)])
#' }
#' whitenoise <- ts2mat(rnorm(128*20),128)
#' # generate x which randomly combine the four types of events with each two of them 
#' # separated by noise
#' x <- c(rnorm(128),t(cbind(shapes,whitenoise)))
#' \dontrun{
#' plot(x,ty='l')
#' }
#' # specify a sliding window size
#' w <- 128
#' # specify a significant level
#' alpha <- 0.05
#' # event detection
#' \dontrun{
#' events <- eventDetection(x,w,'white',parallel=FALSE,alpha, 'art')
#' # clustering
#' cc <- eventCluster(events,4)
#' myclkm <- cc$cl
#' }
#' ##################################
#' #   CASES-99 dataset (9.5m)
#' ##################################
#' # a sliding window length chosen by the user
#' w <- 120 
#' # specify a significant level
#' alpha <- 0.05
#' data(CASES99)
#' \dontrun{
#' CASESevents <- eventDetection(CASES99,w,'red',parallel=FALSE,0.05,'real')
#' cc <- eventCluster(CASESevents,3)
#' cc$center
#' myclkm <- cc$cl
#' # visualise the clustering in 2-dimension PCA space
#' pc.cr <- cc$pca
#' pca.dim1 <- pc.cr$scores[,1]
#' pca.dim2 <- pc.cr$scores[,2]
#' plot(pca.dim1,pca.dim2,col=myclkm+1,main='PCA plots for k-means clustering',pch=16)
#' }


eventCluster <- function(events, k0) {
    x = events$x
    a = events$start
    b = events$end
    eventmeasures = foreach(i = 1:length(a), .combine = rbind) %dopar% {
        event = x[a[i]:b[i]]
        measures(event)
    }
    # transform measures to get them prepared for clustering
    tranmeasures = matrix(0, length(a), dim(eventmeasures)[2])
    for (i in 1:dim(eventmeasures)[2]) {
        x = eventmeasures[, i]
        if (min(x, na.rm = TRUE) <= 0) {
            x = 1 + eventmeasures[, i] - min(eventmeasures[, i], na.rm = TRUE)
        } else {
            x = x
        }
        tranmeasures[, i] = (x)^boxcoxfit(x)$lambda
    }
    # do pca to reduce correlation
    pc.cr <- princomp(scale(tranmeasures))
    ratio_kmeans = rep(NA, 100)
    for (i in 1:100) {
        assign(paste("cc", i, sep = ""), kmeans(pc.cr$scores[, 1:5], k0, iter.max = 1000))
        ratio_kmeans[i] = get(paste("cc", i, sep = ""))$tot.withinss/get(paste("cc", i, sep = ""))$totss
    }
    index = which.min(ratio_kmeans)
    cc = get(paste("cc", index, sep = ""))
    myclkm = cc$cluster
    index = list()
    center = matrix(NA, k0, dim(eventmeasures)[2])
    for (i in 1:k0) {
        index[[i]] = which(myclkm == i)
        if (!is.array(eventmeasures[index[[i]], ])) {
            center[i, ] = (eventmeasures[index[[i]], ])
        } else {
            center[i, ] = colMeans(eventmeasures[index[[i]], ])
        }
    }
    colnames(center) = colnames(eventmeasures)
    return(list(cl = myclkm, center = center, pca = pc.cr))
} 

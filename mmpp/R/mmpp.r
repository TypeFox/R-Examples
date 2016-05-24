#' This packages is used to calculate various similarity and distance measures for sample sequences of both simple and marked temporal point processes. 
#'
#' A simple temporal point process (SPP) is an important class of time series, where the sample realization of the process is solely composed of the times at which events occur. Particular examples of point process data are neuronal spike patterns or spike trains, and a large number of distance and similarity metrics for those data have been proposed. A marked point process (MPP) is an extension of a simple temporal point process, in which a certain vector valued mark is associated with each of the temporal points in the SPP. Analysis of MPPs are of practical importance because instances of MPPs include recordings of natural disasters such as earthquakes and tornadoes.
#' This package implements a number of distance and similarity metrics for SPP, and also extends those metrics for dealing with MPP.
#' It provides a systematic and unified platform for calculating the similarities and distances between SPP, and support marked point process to offer a platform for performing metric-based analysis of earthquakes, tornados, epidemics, or stock exchange data.
#'
#' The package has functions \code{coocmetric}, \code{fmetric}, \code{ieimetric}, and \code{iipmetric} for calculating similarity or distance between two sample sequences. A sample dataset \code{Miyagi20030626} is included in the package. It offers utility two functions: \code{splitMPP}, which splits a sample sequence into a list of partial sequences by using a sliding window, and \code{k2d}, which transforms a similarity matrix to a distance matrix and vice versa.
#'
#' @name mmpp-package
#' @docType package
#' @title A package for Computing Similarity and Distance Metrics for Marked Point Process Data
#' @aliases mmpp
#' @author 
#' Author: Hideitsu Hino \email{hinohide@@cs.tsukuba.ac.jp}, Ken Takano, Yuki Yoshikawa, and Noboru Murata
#' @references R. Quian Quiroga, T. Kreuz, and P. Grassberger. Event synchronization: a simple and fast method to measure synchronicity and time delay patterns, Physical Review E, Vol. 66(4), 041904, 2002.
#' @references J. D. Hunter and G. Milton. Amplitude and frequency dependence of spike timing: implications for dynamic regulation, Journal of Neurophysiology, Vol. 90, pp. 387-94, 2003.
#' @references M. C. W. van Rossum. A Novel Spike Distance. Neural Computation, Vol. 13(4), pp. 751-763, 2001.
#' @references S. Schreiber, J.M. Fellous, P.H. Tiesinga, and T.J. Sejnowski. A new correlation-based measure of spike timing reliability, Neurocomputing, Vols. 52-54, pp. 925-931, 2003.
#' @references T. Kreuz, J.S. Haas, A. Morelli, H.D.I. Abarbanel, and A. Politi. Measuring spike train synchrony, Journal of Neuroscience Methods, Vol. 165(1), pp. 151-161, 2007.
#' @references A.R.C. Paiva, I. Park, and J.C. Principe. A reproducing kernel Hilbert space framework for spike train signal processing, Neural Computation, Vol. 21(2), pp. 424-449, 2009.
#' @examples
#' ## An example to show that the prediction error of the magnitude based on
#' ## 1-nearest neighbor predictor can be reduced by taking marks into account.
#' ## It will take about 5 minutes if you run this example
#' \dontrun{
#' library(mmpp)
#' data(Miyagi20030626)
#' ## split the original MPP by using 3-hour time window
#' sMiyagi <- splitMPP(Miyagi20030626,h=60*60*3,scaleMarks=TRUE)$S
#'
#' ## target of the prediction is the maximum magnitude in the window
#' y <- NULL
#' for(i in 1:length(sMiyagi)){
#'   y <- c(y, max(sMiyagi[[i]]$magnitude))
#' }
#'
#' y <- y[-1]
#' sMiyagi[[length(sMiyagi)]] <- NULL
#'
#' ## number of whole partial MPPs splitted by a 3-hour time window
#' N <- length(sMiyagi)
#' ## training samples are past one week data
#' Ntr <- 24*7/3
#' ## number of different prediction methods
#' Nd <- 10
#'
#' err <- matrix(0, N-Ntr, Nd)
#' colnames(err) <- c("f SPP","iip SPP","cooc smooth SPP","cooc count SPP","iei SPP",
#'                    "f MPP","iip MPP","cooc smooth MPP","cooc count MPP","iei MPP")
#'
#' ## predict the max magnitude in the next 3-hour based on the similarity
#' ## between the current partial point process and the 7-days past partial point process
#' cat("running prediction experiment")
#' for(t in 1:(N-Ntr)){
#'   cat(".")
#'   qid <- Ntr+t
#'   q <- sMiyagi[[qid]]
#'  
#'   ## simple PP
#'   ## fmetric with tau=1
#'   sim2query <- NULL
#'   for(i in 1:Ntr){
#'     sim2query <- c(sim2query,fmetric(q$time,sMiyagi[[qid-i]]$time))
#'   }
#'   err[t,1] <- abs(y[qid]-y[t:(Ntr+t-1)][which.max(sim2query)])
#'  
#'   ## iipmetric with tau=1
#'   sim2query <- NULL
#'   for(i in 1:Ntr){
#'     sim2query <- c(sim2query,iipmetric(q$time,sMiyagi[[qid-i]]$time))
#'   }
#'   err[t,2] <- abs(y[qid]-y[t:(Ntr+t-1)][which.max(sim2query)])
#'
#'   ## coocmetric (smooth) with tau=1
#'   sim2query <- NULL
#'   for(i in 1:Ntr){
#'     sim2query <- c(sim2query,coocmetric(q$time,sMiyagi[[qid-i]]$time,type="smooth"))
#'   }
#'   err[t,3] <- abs(y[qid]-y[t:(Ntr+t-1)][which.max(sim2query)])
#'  
#'   ## coocmetric (count)
#'   sim2query <- NULL
#'   for(i in 1:Ntr){
#'     sim2query <- c(sim2query,coocmetric(q$time,sMiyagi[[qid-i]]$time,type="count"))
#'   }
#'   err[t,4] <- abs(y[qid]-y[t:(Ntr+t-1)][which.max(sim2query)])
#'
#'   ## iei metric
#'   sim2query <- NULL
#'   for(i in 1:Ntr){
#'     sim2query <- c(sim2query,ieimetric(q$time,sMiyagi[[qid-i]]$time))
#'   }
#'   err[t,5] <- abs(y[qid]-y[t:(Ntr+t-1)][which.max(sim2query)])
#'
#'   ## marked PP with latitude, longitude, depth, and magnitude
#'   ## fmetric with tau=1
#'   sim2query <- NULL
#'   for(i in 1:Ntr){
#'     sim2query <- c(sim2query,fmetric(q,sMiyagi[[qid-i]]))
#'   }
#'   err[t,6] <- abs(y[qid]-y[t:(Ntr+t-1)][which.max(sim2query)])
#'  
#'   ## iipmetric with tau=1
#'   sim2query <- NULL
#'   for(i in 1:Ntr){
#'     sim2query <- c(sim2query,iipmetric(q,sMiyagi[[qid-i]]))
#'   }
#'   err[t,7] <- abs(y[qid]-y[t:(Ntr+t-1)][which.max(sim2query)])
#'  
#'   ## coocmetric (smooth) with tau=1
#'   sim2query <- NULL
#'   for(i in 1:Ntr){
#'     sim2query <- c(sim2query,coocmetric(q,sMiyagi[[qid-i]],type="smooth"))
#'   }
#'   err[t,8] <- abs(y[qid]-y[t:(Ntr+t-1)][which.max(sim2query)])
#'  
#'   ## coocmetric (count)
#'   sim2query <- NULL
#'   for(i in 1:Ntr){
#'     sim2query <- c(sim2query,coocmetric(q,sMiyagi[[qid-i]],type="count"))
#'   }
#'   err[t,9] <- abs(y[qid]-y[t:(Ntr+t-1)][which.max(sim2query)])
#'
#'   ## iei metric
#'   sim2query <- NULL
#'   for(i in 1:Ntr){
#'     sim2query <- c(sim2query,ieimetric(q,sMiyagi[[qid-i]]))
#'   }
#'   err[t,10] <- abs(y[qid]-y[t:(Ntr+t-1)][which.max(sim2query)])
#'  
#' }
#' cat("done\n")
#' print(colMeans(err))
#' ##f SPP      iip SPP     cooc smooth SPP  cooc count SPP  iei SPP 
#' ##0.7002634  0.6839529   0.7263602        0.6632930       0.7905148 
#' ##f MPP      iip MPP     cooc smooth MPP  cooc count MPP  iei MPP 
#' ##0.6839529  0.6317594   0.6643804        0.6622056       0.7698548 
#' }
NULL

#' The aftershock data of 26th July 2003 earthquake of M6.2 at the northern Miyagi-Ken Japan
#' 
#' A reparameterization of the \code{main2006JUL26} data frame from the \code{SAPP} package.
#'
#' The variables are as follows:
#'
#' \itemize{
#'   \item time. time in seconds from the main shock.
#'   \item longitude. longitude of the seismic center.
#'   \item latitude. latitude of the seismic center.
#'   \item depth. depth of the seismic center.
#'   \item magnitude. magnitude of the earthquake.
#' }
#'
#' @docType data
#' @keywords datasets
#' @name Miyagi20030626
#' @usage data(Miyagi20030626)
#' @format A data object with 2305 seismic events.
#' @source \code{SAPP} R package available at \url{http://cran.r-project.org}
NULL


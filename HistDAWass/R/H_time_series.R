## functions for analysis and plotting
## smoothing techniques ----
## moving averages -----
#' Smoothing with moving averages of a histogram time series
#' @description (Beta verson of) Extends the moving average smoothing of a time 
#' series to a histogram time series, using L2 Wasserstein distance.
#' @param HTS A \code{HTS} object (a histogram time series).
#' @param k an integer value, the number of elements for moving averages
#' @param weights a vector of positive weights for a weighted moving average 
#' @return a list with the results of the smoothing procedure.
#' @slot k the number of elements for the average
#' @slot weights the vector of weights for smoothing
#' @slot AveragedHTS The smoothed HTS
#' @examples
#' mov.av.smoothed=HTS.moving.averages(HTS=RetHTS, k=5)
#' # a show method for HTS must be implemented you can see it using
#' # str(mov.av.smoothed$AveragedHTS)
#' @export
HTS.moving.averages=function(HTS, k=3, weights=rep(1,k)){
  # check and standardize weights
  if (min(weights)<0){stop("Weights must be positive")}
  weights=weights/sum(weights)
  if (k%%2==0){#the number of elements is even
    stop("for k even it is not yet implemented")
  }
  else{
    epocs=length(HTS@data)-k+1
    AveragedHTS=new("HTS",epocs=epocs)
       counts=0
       shift=(k+1)/2
       for (i in 1:epocs){
         AveragedHTS@data[[i]]@tstamp=HTS@data[[i+shift]]@tstamp
         AveragedHTS@data[[i]]@period=HTS@data[[i+shift]]@period
         tmp=weights[1]*HTS@data[[i]]
         for (j in 2:k){
           tmp=tmp+weights[j]*HTS@data[[i+j-1]]
         }
         #write in the new HTS
         AveragedHTS@data[[i]]@x=tmp@x
         AveragedHTS@data[[i]]@p=tmp@p
         AveragedHTS@data[[i]]@m=tmp@m
         AveragedHTS@data[[i]]@s=tmp@s
         
       }
       
  }
  return(list(k=k,weights=weights,AveragedHTS=AveragedHTS))
  
}
## exponential smoothing ----
#' Smoothing with exponential smoothing of a histogram time series
#' @description (Beta verson of) Extends theexponential smoothing of a time 
#' series to a histogram time series,using L2 Wasserstein distance.
#' @param HTS A \code{HTS} object (a histogram time series).
#' @param alpha a number between 0 and 1 for exponential smoothing 
#' @return a list with the results of the smoothing procedure.
#' @slot smoothing.alpha  the alpha parameter
#' @slot AveragedHTS The smoothed HTS
#' @examples
#' mov.expo.smooth=HTS.exponential.smoothing(HTS=RetHTS, alpha=0.8)
#' # a show method for HTS must be implemented you can see it using
#' # str(mov.expo.smooth$AveragedHTS)
#' @export
HTS.exponential.smoothing=function(HTS, alpha=0.9){
  if ((alpha<0)||(alpha>1)) {stop("The smoothing parameter alpha must be between 0 and 1")}
  epocs=length(HTS@data)-1
  SmoothHTS=new("HTS",epocs=epocs)
  s=HTS@data[[1]]
  for (i in 2:epocs){
    SmoothHTS@data[[i-1]]@tstamp=HTS@data[[i]]@tstamp
    SmoothHTS@data[[i-1]]@period=HTS@data[[i]]@period
    s=(1-alpha)*s+alpha*HTS@data[[i]]
    SmoothHTS@data[[i-1]]@x=s@x
    SmoothHTS@data[[i-1]]@p=s@p
    SmoothHTS@data[[i-1]]@m=s@m
    SmoothHTS@data[[i-1]]@s=s@s
  }
  return(list(smoothing.alpha=alpha, SmoothedHTS=SmoothHTS))
}
## prediction techniques ----
## predicting via k-nn ----
#' K-NN predictions of a histogram time series
#' @description (Beta verson of) Extends the K-NN algorithm for predicting a time 
#' series to a histogram time series, using L2 Wasserstein distance.
#' @param HTS A \code{HTS} object (a histogram time series).
#' @param position an integer, the data histogram to predict
#' @param k the number of neighbours (default=3) 
#' @return a \code{distributionH} object predicted from data.
#' @references Javier Arroyo, Carlos Mate, Forecasting histogram time series with k-nearest neighbours methods, 
#' International Journal of Forecasting, Volume 25, Issue 1, January-March 2009, Pages 192-207,
#'  ISSN 0169-2070, http://dx.doi.org/10.1016/j.ijforecast.2008.07.003.\cr
#' @details Histogram time series (HTS) describe situations where a distribution of values is available 
#' for each instant of time. These situations usually arise when contemporaneous or temporal aggregation 
#' is required. In these cases, histograms provide a summary of the data that is more informative than those 
#' provided by other aggregates such as the mean. 
#' Some fields where HTS are useful include economy, official statistics and environmental science. 
#' The function adapts the k-Nearest Neighbours (k-NN) algorithm to forecast HTS and, more generally,
#'  to deal with histogram data. The proposed k-NN relies on the L2 Wasserstein distance that is used 
#'  to measure dissimilarities between sequences of histograms and to compute the forecasts. 
#' 
#' @examples
#' prediction=HTS.predict.knn(HTS=RetHTS, position=108, k=3)
#' @export
HTS.predict.knn=function(HTS,position=length(HTS@data), k=3){
  if (length(HTS@data)>position){HTS@data=HTS@data[1:position]}
  #compute distances
  dista=vector(mode="numeric",position-1)
  for (i in 1:(position-1)){
    dista[i]=WassSqDistH(HTS@data[[position]],HTS@data[[i]])
  }
  firstksorted=sort(dista)[1:k]
  positions=vector(mode="numeric",k)
  vectorofdistr=new(Class = "MatH",nrows=k,ncols=1)
  for (v in 1:k){
  positions[v]=which(dista==(firstksorted[v]))+1
  vectorofdistr@M[v,1][[1]]=new("distributionH",
                                HTS@data[[positions[v]]]@x,
                                HTS@data[[positions[v]]]@p,
                                HTS@data[[positions[v]]]@m,
                                HTS@data[[positions[v]]]@s)
  }
  #compute the mean of histograms at position positions
  pred=WH.vec.mean(vectorofdistr)
  return(pred)
}
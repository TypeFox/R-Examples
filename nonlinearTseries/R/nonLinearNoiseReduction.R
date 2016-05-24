
# References: http://www.mpipks-dresden.mpg.de/~tisean/TISEAN_2.1/docs/chaospaper/node23.html#SECTION00061000000000000000
#' Nonlinear noise reduction
#' @description
#' Function for denoising a given time series using nonlinear analysis techniques. 
#' @details
#' This function takes a given time series and denoises it. The denoising
#' is achieved by averaging each Takens' vector in an m-dimensional space
#' with his neighbours (time lag=1). Each neighbourhood is specified with balls of a given radius
#' (max norm is used).
#' @param time.series The original time series to denoise.
#' @param embedding.dim Integer denoting the dimension in which we shall embed the \emph{time.series}.
#' @param radius The radius used to looking for neighbours in the phase space (see details).
#' @return A vector containing the denoised time series.
#' @references H. Kantz  and T. Schreiber: Nonlinear Time series Analysis (Cambridge university press)
#' @author Constantino A. Garcia
#' @rdname nonLinearNoiseReduction
#' @export nonLinearNoiseReduction
#' @useDynLib nonlinearTseries
nonLinearNoiseReduction=function(time.series, embedding.dim, radius){
  #build takens' vectors using time.lag = 1
  takens=buildTakens(time.series=time.series,embedding.dim=embedding.dim,time.lag=1) 
  
  denoised.time.series = .C("nonlinearNoiseReduction",timeSeries = as.double(time.series), 
                                  takens = as.double(takens),
                                  numberTakens = as.integer(nrow(takens)),
                                  embeddingD = as.integer(embedding.dim),
                                  eps = as.double(radius),
                                  numberBoxes = as.integer(400),
                                  PACKAGE="nonlinearTseries")
  
  return(denoised.time.series$timeSeries)
  
}

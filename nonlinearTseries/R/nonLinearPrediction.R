# References: http://www.mpipks-dresden.mpg.de/~tisean/TISEAN_2.1/docs/chaospaper/node23.html#SECTION00061000000000000000
#' Nonlinear time series prediction
#' @description
#' Function for predicting futures values of a given time series using previous values
#' and  nonlinear analysis techniques. 
#' @details
#' Using \emph{time.series} measurements, an embedding in \emph{embedding.dim}-dimensional phase space with
#' time lag \emph{time.lag} is used to predict the value following the given time series after \emph{prediction.steps}
#' sample steps. This is done by finding all the neighbours of the last Takens' vector in a  radius of size \emph{radius} (the
#' max norm is used). If no neighbours are found within a distance radius, the neighbourhood size
#' is increased until succesful using \emph{radius.increment}(\emph{radius} = \emph{radius} + \emph{radius.increment}).
#' @param time.series Previous values of the time series that the algorithm will use to make the prediction.
#' @param embedding.dim Integer denoting the dimension in which we shall embed the \emph{time.series}.
#' @param time.lag Integer denoting the number of time steps that will be use to construct the 
#' Takens' vectors.
#' @param radius The radius used to looking for neighbours in the phase space (see details).
#' @param radius.increment The increment used when no neighbours are found (see details).
#' @param prediction.step Integer denoting the number of time steps ahead for the forecasting.
#' @return The predicted value \emph{prediction.step} time steps ahead.
#' @references H. Kantz  and T. Schreiber: Nonlinear Time series Analysis (Cambridge university press)
#' @examples
#' \dontrun{
#' h=henon(n.sample=5000,start=c(0.324,-0.8233))
#' predic=nonLinearPrediction(time.series=h$x[10:2000],embedding.dim=2,
#'                            time.lag=1,
#'                            prediction.step=3,radius=0.03,
#'                            radius.increment=0.03/2)
#' cat("real value: ",h$x[2003],"Vs Forecast:",predic)
#' }
#' @author Constantino A. Garcia
#' @rdname nonLinearPrediction
#' @export nonLinearPrediction
#           http://www.mpipks-dresden.mpg.de/~tisean/TISEAN_2.1/docs/chaospaper/node18.html#SECTION00052000000000000000
nonLinearPrediction=function(time.series,embedding.dim,time.lag,prediction.step,radius,radius.increment){
  nfound=0
  av=0
  l=length(time.series)
  #vector of lag delays used to build the takens' vectors
  jumpsvect = seq((embedding.dim - 1) * time.lag,0, -time.lag)
  #reference takensVector
  takensVector=time.series[l-jumpsvect]
  #first position that we can use for construct a 'reverse' takens' vector: t(n)=[t(n-(m-1)*time.lag),t(n-(m-2)*time.lag),...,t(n)]
  beg=(embedding.dim-1)*time.lag+1
  #last position having into account that we want to use it to predict prediction.steps
  #steps forward
  en=l-prediction.step
  #whilen  no neighbours are found, we increase the size of the neighbourhood
  while (nfound==0){
    #build takens vectors and check if there exist some neighbour of takensVector
    for (i in beg:en){
        if (isNeighbour(takensVector,time.series[i-jumpsvect],embedding.dim,radius)){
          #average of predictions
          av=av+time.series[[i+prediction.step]]
          nfound=nfound+1
        }
    }
    #increment the size of the neighbourhood
    radius=radius+radius.increment   
  }
  #return prediction
  return(av/nfound)
  
}
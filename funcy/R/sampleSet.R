#
#  Copyright (C) 2011-2015 Christina Yassouridis
#  
#

##Center functions
funcyLib <- function(){
    Fun1 <- function(x){x^2}
    Fun2 <- function(x){sqrt(x)}
    Fun3 <- function(x){sin(2*pi*x)}
    Fun4 <- function(x){x^3}
    Fun5 <- function(x){-x^2}
    Fun6 <- function(x){(x-1)}
    return(fcts=list(Fun1, Fun2, Fun3, Fun4, Fun5, Fun6))
}

##Generating data samples
sampleFuncy <- function(obsNr=100,
                        k=4,
                        timeNr=20,
                        timeNrMax=NULL,
                        timeNrMin=NULL,
                        timeInterval=c(0,1),
                        nrGridPts=30,
                        sd = 0.3,
                        reg=TRUE
                        ){
    
    if(is.null(timeNr) & reg==TRUE)
        stop("Please determine the number of evaluation points in timeNr!")
    else if(is.null(timeNrMax) & is.null(timeNrMin) & reg==FALSE)
        stop("timeNrMin and timeNrMax have to be specified.")
   
    meanFcts <- funcyLib()[sample(1:6,k)]
    if(reg)
        timeNrMin <- timeNrMax <- timeNr
    T <- timesamples(obsNr, timeNrMax, timeNrMin, timeInterval, nrGridPts, reg)
    N <- T$timeNrIndividual
    time <- T$time
  
    x <- matrix(0, nrow=obsNr, ncol=timeNrMax)
    eval <- NULL
    class <- rep(0,k)
    cl <- 1
    for(l in 1: k){
        Until <- floor(obsNr*l/k)
        if(l==k) Until <- obsNr
        
        for(i in (floor(obsNr*(l-1)/k)+1):Until)  {
            x[i,1:N[i]]<-meanFcts[[l]](time[i,1:N[i]])+rnorm(meanFcts[[l]](time[i,1:N[i]]), sd=sd)
            eval <- c(eval,x[i,1:N[i]])
            class[i] <- cl
        }
        cl <- cl+1
    }

  object <- new("sampleFuncy")
  object@clusters=class
  object@reg=reg
  if(reg==FALSE){
      curveIndex <- rep(1:obsNr,N)
      timevector <- unlist(lapply(seq.int(length(N)), function(x)
          time[x,1:N[x]]))
      object@data <- cbind(curveIndx=curveIndex, eval=eval,time=timevector)
  }else
      object@data=t(x)
  
  return(object)
}

##Points on irregular time interval and different time points
timesamples <- function(obsNr,  timeNrMax, timeNrMin, timeInterval, nrGridPts, reg = FALSE){
    T <- matrix(0, nrow = obsNr, ncol = timeNrMax)
    t <- matrix(0, nrow = obsNr, ncol = timeNrMax)
    N <- rep(0,obsNr)
    seq <- seq(from = timeInterval[1], to=timeInterval[2], length.out = nrGridPts)
    if(!reg & (timeNrMin!= timeNrMax)) number <- sample(timeNrMin:timeNrMax, obsNr, replace=TRUE)
    else  number <- rep(timeNrMax, obsNr)

    for(j in 1:obsNr){
        if(!reg)  timeIndex <- sort(sample(1:nrGridPts,number[j]))
        else timeIndex <-floor(seq(from= 1, to = nrGridPts, by =(nrGridPts-1)/(timeNrMax-1)))
        T[j,1:length(timeIndex)] <- timeIndex
        t[j,1:sum(T[j,]!=0)] <- seq[T[j,]]
    }

    N <- apply(T, 1, function(x) max(which(x>0)))

    return(list(timeIndex=T, time=t, timeNrIndividual=N))
}

setGeneric("Cluster",
           function(object)
               standardGeneric("Cluster")
           )
           
setMethod("Cluster", "sampleFuncy",
          function(object)
              object@clusters
          )

setGeneric("Data",
           function(...) utils::data(...))

setMethod("Data","sampleFuncy",
         function(object,...)
             object@data
         )





#'EstTime
#'
#'\itemize{
#'  \item Simulates k state-transition tables (see: \code{\link{StateTrans}}) based one state-transtion table containing expected population frequencies using the \code{\link{EstFreq}} function. 
#'  \item Repeats simulations for a number of defined time intervals or until a termination criterion is reached.
#'  \item The estimated number of cells with low or zero frequencies are computed.
#'  \item The proportion of low/zero frequencies is plotted against the number of time points/intervalls. 
#'}
#'
#'@param x a matrix containing the assumed probabilities (see \code{\link{EstFreq}} for further detail!)
#'@param t optional: a vector of time intervalls for which frequencies should be simulated. Can be used to decrease simulation time!
#'@param crit optional: but must be sepcified if t is not; 
#'Simulations will end if the relative frequency of zero cells is less than crit
#'@param zero optional: if FALSE simulations will end if number of low frequencies 
#'(instead of zero frequencies) is less than crit
#'@param k Number of simulations (at least 20.000 is recommended)
#'@param min.cell a single integer defines what counts as a low frequency (lower than 5 by convention)
#'@param pos position of the output's legend. 
#'Options are: "bottomleft", "bottomright", "upperleft", and "upperright".
#'@param smoothed logical value. If true, output lines will be smoothed!
#'@param show.it single integer that defines which steps of iteration protocol should be shown. Only active if t is not defined else iteration protocoll is replaced with a progression bar
#'@param max.it single integer that defines the maximum number of iterations if t is not specified.
#'
#'@return EstTime object; a list of three vectors
#'if printed, provides a plot of expected number of low and zero cell frequencies
#'
#'@details First vector represents time points, second vector provides rel. frequency of cases with zero's, third vector rel. frequency of cases with low cell frequencies.
#'
#'
#'@examples
#'\dontrun{
#'my.trans.table<-matrix(c(0.57, 0.13,0.05,0.05,0.05, 0.05,0.05,0.05),4,2)
#'my.cellproblems<-EstTime(my.trans.table, k=500) # low k only for examplification purposes!
#'my.cellproblems
#'}
#'@export



EstTime<-function(x,t=NA,crit=0.05,zero=T,min.cell=5, k=20000, pos="bottomleft", smoothed=T, show.it=10, max.it=10000){

  zero<-numeric(length(t))
  low<-numeric(length(t))
  
  if(!(any(is.na(t)))){
    mes<-paste("\n\nThis function will run a total of", format(length(t)*k, scientific=FALSE), "simulations! \nComputation may take some minutes!")
    cat(mes)
    Sys.sleep(1)

    cat("\n\nProgress:")

    pb <- utils::txtProgressBar(min = 0, max = length(t), style = 3)


    for(i in 1:length(t)){
      a<-EstFreq(x, t[i], min.cell, k)
      zero[i]<-1-mean(a[[1]]==0)
      low[i]<-1-mean(a[[2]]<min.cell)
      Sys.sleep(0.001)
      utils::setTxtProgressBar(pb, i)
    }
    t<-t
  } # If ZERO is used!
  if(any(is.na(t))&&zero){
    mes<-paste("\n\nThis function will run an unknown number of simulations! \nComputation may take some minutes!\n")
    cat(mes)
    Sys.sleep(1)

    t1<-seq(1, max.it, 5)
    timer<-integer(1)

    for(i in 1:length(t1)){
      a<-EstFreq(x, t1[i], min.cell, k)
      zero.break<-1-mean(a[[1]]==0)
      timer<-i
      if(i%%show.it==0) cat(paste("Preprocessing Cycle: ",timer, "zeros", round(zero.break,2), "crit: zero<.", crit, "\n"))
      if(zero.break<crit) break
      Sys.sleep(0.001)
    }

    t<-round(1*t1[timer]/3,0):max.it

    for(i in 1:length(t)){
      a<-EstFreq(x, t[i], min.cell, k)
      zero[i]<-1-mean(a[[1]]==0)
      low[i]<-1-mean(a[[2]]<min.cell)
      if(i%%show.it==0) cat(paste("Main-Process, time number: ",t[i],"zeros:",round(zero[i],2), "low freq:", round(low[i],2), "crit: zero<.",crit, "\n" ))
      if(zero[i]<crit) break
      Sys.sleep(0.001)
    }
   t<-t[1:length(zero)]
   cat(paste("End of simulation. Crit for zero found at", max(t)))
  }
  #### IF LOW is used! // if ZERO is FALSE!
  if(any(is.na(t))&&!zero){
    mes<-paste("\n\nThis function will run an unknown number of simulations! \nComputation may take some minutes!\n")
    cat(mes)
    Sys.sleep(1)

    t1<-seq(1, max.it, 5)
    timer<-integer(1)

    for(i in 1:length(t1)){
      a<-EstFreq(x, t1[i], min.cell, k)
      low.break<-1-mean(a[[2]]<min.cell)
      timer<-i
      if(i%%show.it==0) cat(paste("Preprocessing Cycle: ",timer, "low freq.", round(low.break,2), "crit: zero<.", crit, "\n"))
      if(low.break<crit) break
      Sys.sleep(0.001)
    }

    t<-round(1*t1[timer]/3,0):max.it

    for(i in 1:length(t)){
      a<-EstFreq(x, t[i], min.cell, k)
      zero[i]<-1-mean(a[[1]]==0)
      low[i]<-1-mean(a[[2]]<min.cell)
      if(i%%show.it==0) cat(paste("Main-Process, time number: ",t[i],"zeros:",round(zero[i],2), "low freq:", round(low[i],2), "crit: low<.",crit, "\n" ))
      if(low[i]<crit) break
      Sys.sleep(0.001)
    }
    t<-t[1:length(zero)]
    cat(paste("End of simulation. Crit for low freq. found at", max(t)))
  }

  output<-list(t,zero, low)
  attributes(output)$pos<-pos
  class(output)<-"EstTime"
  return(output)
}



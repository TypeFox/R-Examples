ProbBins <- function(ano, fcyr, thr, quantile=TRUE, posdates = 3,
                     posdim = 2, compPeriod= "Full period") { 
  # define dimensions
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  nbdim <- length(dim(ano))
  nfcyr<-length(fcyr)

  if (nbdim < 7){
    ano <- Enlarge(ano, 7)
  }

  #remove posdates of posdim and supress duplicates
  posdim <- setdiff(posdim, posdates)
  nbpos <- length(posdim)
  #permute dimensions in ano
  ano <- aperm(ano, c(posdates, posdim, setdiff(seq(1,7,1), c(posdates, posdim))))
  dimano <- dim(ano)
  nsdates=dimano[1]
  #calculate the number of elements on posdim dimension in ano
  nmemb=1
  if (nbpos > 0){
    for (idim in 2:(nbpos+1)){
      nmemb=nmemb*dimano[idim]
    }
  }
  
  # separate forecast and hindcast
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  fore <- array(ano[fcyr, , , , , , ], dim = c(nfcyr,
                                               dimano[2:7]))
  # the members and startdates are combined in one dimension
  sample_fore <- array(fore, dim=c(nfcyr*nmemb, dimano[(nbpos+2):7]))
  
  if(compPeriod=="Full period"){
    hind <- ano
    sample <- array(hind, dim=c(nsdates*nmemb, dimano[(nbpos+2):7]))
  }
  
  if (compPeriod=="Without fcyr"){
    hind <- array(ano[-fcyr, , , , , , ], dim = c(nsdates-nfcyr,
                                                 dimano[2:7]))
    sample <- array(hind, dim=c((nsdates-nfcyr)*nmemb, dimano[(nbpos+2):7]))
  }
  
  #quantiles for each grid point and experiment
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  if (quantile==TRUE){
    qum <- apply(sample, seq(2,7-nbpos,1), FUN=quantile,probs=thr,na.rm=TRUE,names=FALSE,type=8)                                        
  }else{
    qum<-array(thr,dim=c(length(thr), dimano[(nbpos+2):7]))
  }

  # This function assign the values to a category which is limited by the thresholds
  # It provides binary information
  
  counts <- function (dat, nbthr){
    thr <- dat[1:nbthr]
    data <-  dat[nbthr+1:(length(dat)-nbthr)]
    prob <- array(NA, dim=c(nbthr+1, length(dat)-nbthr))
    prob[1,]=1*(data <= thr[1])
    if(nbthr!=1){
      for (ithr in 2:(nbthr)){
        prob[ithr,]=1*((data > thr[ithr - 1]) & (data <= thr[ithr]))
      }
    }
    prob[nbthr+1,]=1*(data > thr[nbthr])
    return(prob) 
  }
      
  # The thresholds and anomalies are combined to use apply
  data <- abind(qum, sample_fore, along = 1)
  
  # PBF:Probabilistic bins of a forecast.
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # This array contains zeros and ones that indicate the category where your forecast is. 

  PBF <- array(apply(data, seq(2,7-nbpos,1), FUN=counts, nbthr=length(thr)),
               dim=c(length(thr)+1, nfcyr, nmemb, dimano[(nbpos+2):nbdim]))
  
  return(PBF)
  
  
  if (compPeriod == "cross-validation"){
    for (iyr in 1:fcyr){
      ProbBins(ano,fcyr=iyr,thr=thr,posdates=posdates,posdim=posdim)
    } 
  }
}

#' @export logLik.pers
#' @title S3 logLik for Object of class "pers"
#' @description S3 logLik method to extract the log-likelihood for object of class\code{"pers"} 
#' @param object object of class\code{"pers"}
#' @param sat a "logical" with default set to \code{sat=FALSE} to return the Log-Likelihood of the data for the unrestricted modell based on parameters estimated with function \code{\link{pers}}. If set to \code{sat=TRUE} the Log-Likelihood of the saturated model is returned instead.
#' @param ... not used jet.

########################### hier die logLik method fuer pers #############################
logLik.pers<-function(object, sat=FALSE, ...){
  
  # compute m_v - number of categories per item as vector - used in both cases
  #m_v <-sapply(1:nrow(object$pair$threshold), function(i) {length(na.omit(object$pair$threshold[i,]))+1})
  m_v <- object$pair$m # new 03-12-2015 m is in any case part of "pair"
  
  if(sat==FALSE){ 
    # estimated model # NEW and correct since 21.11.2015 checked against WinMira (dichotom)
    #str.pattern(object$pair$resp)
   P <- pvx.super(theta_v=object) # probabilities for hole dataset

    #     pat_ind <- rownames(object$pair$resp)%in%(rownames(unique(object$pair$resp))) # index for unique pattern 
#     pat_fre <- (table(str.pattern(object$pair$resp))) # frequencies for pattern
#     pat_lik <- (rowSums(log(P[pat_ind,]),na.rm=T)) # likelihood for (unique) pattern
#     Log_Likelihood <- sum( pat_lik * pat_fre) # likelihood data (conditional?)
#     # old things prior 21.1.2015
    #     #P <- matrix(rep(.1,5000),ncol = 5,nrow = 1000) #test
    #     sum(log(P),na.rm=T)
      Log_Likelihood <- sum(log(P),na.rm=T) # evtl. : P[complete.cases(P),]    
#     pat_fre[1] <- 0
#     pat_fre[30] <- 0
#     
#     sum( pat_lik * pat_fre)
       
    # neu auslesen der WLE geschätzten likelihood
       # Log_Likelihood <- sum(unique(object$pers$WLL)) #sum(as.numeric(names(table(object$pers$WLL))))

    
    # object$pair$threshold
    df <- sum(m_v-1)+sum(m_v-1)-1 # changed 27-3-2015 (acc. WINMIRA) (sum(m_v-1)-1)*2 # number of free model parameters (df) (number of scorgroups - 1 ): sum(m_v-1)  +  (numper of free itemparameter): sum(m_v-1)-1  
    #     df_k <- sum(m_v-1)-1
    #     df_k_1 <- sum(m_v-1)-1
    nall <- dim(object$pair$resp)[1]
    nobs <- nall ### ev. check this !!! OK
    #structure(-213.902975399879, nall = 108L, nobs = 108, df = 8, class = "logLik")
    result <- structure(Log_Likelihood, nall = nall, nobs = nobs, df = df, class = "logLik") 
  }
  
  if(sat==TRUE){
    # saturated model
    # tested against WinMira 20-04-2013 --> OK!!!  
    x<-object$pair$resp
    #df.sat<-(m^k)-1 # degrees of fredom for saturated model
    df <- prod(m_v)-1# changed 27-3-2015 (acc. WINMIRA) prod(m_v) # changed for polytomies
    b<-dim(x)[2]
    l<-dim(x)[1]
    ppaste <- function(z){paste(z, collapse = "")} # creates pattern strings
    zaehl<-(as.matrix(table(apply(x, 1, ppaste))))[ ,1] # compute values for numerator and exponent
    nen<-rep(l,times=length(zaehl))     # compute values for denominator
    lik.sat<-sum(log((zaehl/nen)^zaehl)) # compute log likelihood of saturated model
    result <- structure(lik.sat, nall = l, nobs = l, df = df, class = "logLik")   
  }
  return(result)
}
# Filename: logLik_hapassoc.R
# Version : 

# HapAssoc- Inference of trait-haplotype associations in the presence of uncertain phase
# Copyright (C) 2003  K.Burkett, B.McNeney, J.Graham

# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

########################################################################


logLik.hapassoc <- function(object,...){
  dfs <- dim(object$var)[1]
  loglik <- object$loglik
  cat(paste("Log Lik: ",round(loglik,4)," (df=",dfs,")\n",sep=""))
  invisible(list(logLik=loglik, df=dfs))
}


anova.hapassoc <- function(object,redfit,display=TRUE,...){

  fullfit <- object
  full.like <- fullfit$loglik
  red.like <- redfit$loglik

  dfs <- dim(fullfit$var)[1]-dim(redfit$var)[1]
  LRTstat <- 2*(full.like-red.like)
  pval <- pchisq(LRTstat,df=dfs,lower.tail=FALSE)

  if (display==TRUE){
    cat("\n\thapassoc: likelihood ratio test\n\n")
    cat(paste("Full model:",fullfit$model[2],fullfit$model[1],
              fullfit$model[3],"\n"))
    cat(paste("Reduced model:",redfit$model[2],redfit$model[1],
              redfit$model[3],"\n\n"))
    cat(paste("LR statistic =",round(LRTstat,4),", df =",dfs,", p-value = ",
              round(pval,4),"\n"))
  }
  
  invisible(list(LRTstat=LRTstat,df=dfs,pvalue=pval))
}



#-------------------------------------------------------------------#

## log.likelihood is now called in hapassoc


loglikelihood <- function(haplos.list, EMresults){

  
  # Get previous results
  haplos<-haplos.list$haploMat
  freq <- EMresults$freq
  IDs <- EMresults$ID

  # Get the P(X) part. These two commands came from function hapassoc
  haplo.probs <- rep(1, nrow(haplos)) + isMultiHetero(haplos.list)
  haplo.probs <- haplo.probs * freq[haplos[, 1],] * freq[haplos[,2],]
  names(haplo.probs) <- NULL

  # Determine vector P(Y|X) for all individuals using the function
  # probY and multiply by haplo.probs
  like.i <- probY(EMresults)*haplo.probs

  # Calculate log likelihood over all individuals by first summing
  # within ID, then taking log.
  value <- sum(log(tapply(like.i,IDs,sum)))

  return(value)
}

#--------------------------------------------------------------------#


## Use the R dnorm, dbinom, dpois and dgamma functions to
## get these values. Use the "fits", as they are the fitted means


probY <- function(myfit){

  y <- myfit$response
  fits <- myfit$fits

  if(myfit$family$family=="binomial"){
    probY<- dbinom(y,1,prob=fits) }

  if(myfit$family$family=="poisson"){ 
    probY <- dpois(y,lambda=fits)}

  if(myfit$family$family=="gaussian") {
    probY<- dnorm(y,mean=fits, sd=sqrt(myfit$dispersion))}

  if(myfit$family$family=="Gamma"){
    a = 1/myfit$dispersion
    s = fits/a
    probY <- dgamma(y,shape=a,scale=s)}

  return(probY)
}



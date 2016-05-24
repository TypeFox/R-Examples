# Filename: summaryHap.R
# Version : $Id: summaryHap.R,v 1.13 2008/07/18 22:43:35 mcneney Exp $

# HapAssoc- Estimation of trait-haplotype associations in the presence of uncertain phase
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

summary.hapassoc<-function(object, ...) {
names(object)[names(object) == "freq"] <- "gamma"
if(object$converged==FALSE) { warning("hapassoc failed to converge\n"); return() }
family<-object$family$family
if(family=="Gamma"){
 #For a Gamma model, use a moment rather than ML estimate of phi.
 #The ML estimate is overly sensitive to round-off errors in small 
 #responses and to departures from the Gamma model (Mc&N pp295-296).
  dispersion=momentPhiGamma(object) 
  object$var<-object$var*dispersion/(object$dispersion)
} else {
  dispersion=object$dispersion
}

numbeta<-length(object$beta)
if(family=="Gamma"|family=="gaussian") #remove row and col for dispersion
   object$var<-object$var[-(numbeta+1),-(numbeta+1)]
numfreqs<-length(object$gamma)

if(numbeta>1)
  coef.table<-cbind(object$beta,sqrt(diag(object$var[1:numbeta,1:numbeta])))
else
  coef.table<-cbind(object$beta,sqrt(object$var[1:numbeta,1:numbeta]))

#Now calculate the se for the last frequency estimator
varfreqs<-object$var[(numbeta+1):(numbeta+numfreqs-1),(numbeta+1):(numbeta+numfreqs-1)]
cvec<-rep(-1,numfreqs-1)
varlast<-t(cvec)%*%varfreqs%*%cvec
freq.table<-cbind(object$gamma, c(sqrt(diag(varfreqs)),sqrt(varlast)))

#Compute z-scores for regression coefficients
coef.table<-cbind(coef.table,coef.table[,1]/coef.table[,2])

#Now two-sided p-values. 
coef.table<-cbind(coef.table,1-pchisq(coef.table[,3]^2,df=1))
dimnames(coef.table)<-list(names(object$beta),
                           c("Estimate","Std. Error","zscore","Pr(>|z|)"))
dimnames(freq.table)<-list(
      paste("f.",dimnames(object$gamma)[[1]],sep=""),
      c("Estimate","Std. Error"))

sumList <- list(call=object$call,subjects=length(unique(object$ID)),
            coefficients=coef.table,frequencies=freq.table,
            dispersion=dispersion,loglik=object$loglik,
            family=object$family$family)
class(sumList) <- "summary.hapassoc"
return(sumList)
}
## Other functions called in summary.hapassoc

########################################################################

print.summary.hapassoc<-function(x, ...){
  cat("Call:\n")
  print(x[[1]])

  cat("\nNumber of subjects used in analysis:",x[[2]],"\n")

  cat("\nCoefficients:\n")
  print(x[[3]])

  cat("\nEstimated haplotype frequencies:\n")
  print(x[[4]])

  cat(paste("\n(Dispersion parameter for ",x$family," family taken to be ",x$dispersion,")\n",sep=""))

  cat(paste("\nLog-likelihood:",x$loglik,"\n\n"))
}

momentPhiGamma<-function(object) {
  ans<-0
  uniqueID<-unique(object$ID)
  for(i in 1:length(uniqueID)) {
    ind<-(object$ID==uniqueID[i]) #identify pseudo-individuals for uniqueID[i]
    myy<-(object$response[ind])[1] #response is the same for all such ps-indiv
    mywt<-object$wts[ind]; myfit<-object$fits[ind]
    mubar<-sum(mywt*myfit)
    a<-sum(mywt*(myfit-mubar)^2)
    b<-sum(mywt*myfit^2)
    ans<-ans+((myy-mubar)^2-a)/b
  }
  n<-sum(object$wts); p<-length(object$beta)+length(object$gamma)
  return(ans/(n-p))
}

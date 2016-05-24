#  Part of the R package, http://www.R-project.org
#
#  Copyright (C) 1995-2012 The R Core Team
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  http://www.r-project.org/Licenses/

computeLoglikeFromPartition <- function(labels, expdata){
  
  if(is.matrix(expdata) == FALSE & is.data.frame(expdata) == FALSE) 
    stop(paste(sQuote("expdata"), "must be a matrix"))
  
  if(is.list(labels) == FALSE & is.vector(labels) == FALSE) 
    stop(paste(sQuote("labels"), "must be a list"))
  
  
  ## function computeLoglikeFromPartition
  ## input labels (list): block labels for each variable
  ## input expdata (matrix): data 
  ## output loglike (real): loglikehood of the model with the block diagonal covariance
  ## output df (integer): degree of freedom of the model
  ## output labels (list): labels provided as input
  
  ## require
  ## library(mvtnorm)
  
    covSSS <-matrix(rep(0,dim(expdata)[2]*dim(expdata)[2]),ncol=dim(expdata)[2], nrow=dim(expdata)[2])
    for(lab in 1:length(unique(labels))){
        covSSS[which(labels==lab),which(labels==lab)] <-1/(dim(expdata)[1])*t(expdata[,which(labels==lab)])%*%expdata[,which(labels==lab)]
    }
    dK <- matrix(table(labels))
    return(list(loglike=sum(dmvnorm(expdata, mean=rep(0,dim(expdata)[2]), sigma=covSSS,log=TRUE)), df=sum(dK*(dK-1)/2), labels=labels))
}



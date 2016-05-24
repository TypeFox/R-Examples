# Author : F.Rohart,  Australian Institute for Bioengineering and Nanotechnology, The University of Queensland, Brisbane, QLD
# created: pre 01-01-2013
# last modification: 14-03-2015
# Copyright (C) 2014
#
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 2
# of the License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.

predict.mht <- predict.mht.order <-
function(object,newx,level,...)
{
    if(missing(level)) level=0.95

	 if(missing(newx))
     {
         out=object$fitted.values
         pred.alpha=out
	  }else{
          
          p=ncol(newx)
          if(is.null(colnames(newx)))
          {
              colnames(newx)=paste("X", 1:p,sep="")
          }
          
          n=nrow(newx)
          means.X=object$data$means.X
          sigma.X=object$data$sigma.X
		  
          if(sum(sum(newx[,1]==1)==n))
          {
              newx=scale(newx,center=means.X,scale=sigma.X)
			  newx[which(is.na(newx))]=0
			  intercept=TRUE
              
		  }else{
			  message("intercept has been added")
			  
              newx=scale(newx,center=means.X,scale=sigma.X)
			  newx[which(is.na(newx))]=0
			  newx=cbind(1,newx)
              colnames(newx)=c("Intercept",colnames(newx)[-1])
			  intercept=FALSE
		  }
          
          if(!identical(colnames(object$data$X),colnames(newx))) stop("colnames of object$data$X and newx must be the same")
		  
        alpha=as.numeric(colnames(object$coefficients))

      
        pred.alpha=array(0,c(nrow(newx),3,length(alpha)))
        for(i in 1:length(alpha)) #calculation of the CI/PI for each alpha
        {
          fitData=data.frame(object$data$X[,object$relevant_var[i,],drop=FALSE])
          names(fitData)=colnames(object$data$X[,object$relevant_var[i,],drop=FALSE])
          reg=lm(object$data$Y~.-1,data=fitData)
          predData=data.frame(newx[,object$relevant_var[i,],drop=FALSE])
          names(predData)=colnames(newx[,object$relevant_var[i,],drop=FALSE])
          pred=predict(reg,predData,interval="prediction")
          pred.alpha[,,i]=pred
        }
        if(!is.null(rownames(newx)))
        {dimnames(pred.alpha)[[1]]=rownames(newx)}
        
        dimnames(pred.alpha)[[2]]=colnames(pred)
        dimnames(pred.alpha)[[3]]=paste("alpha=",alpha)
        
        pred.alpha=pred.alpha[,1:3,1:length(alpha),drop=FALSE]
        


      }
      
pred.alpha
}

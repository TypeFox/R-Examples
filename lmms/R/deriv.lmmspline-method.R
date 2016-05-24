# Jasmin Straube, Queensland Facility of Advanced Bioinformatics
# Part of this script was borrowed from parallel and snow package.
#
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU Moleculesral Public License
# as published by the Free Software Foundation; either version 2
# of the License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Moleculesral Public License for more details.
#
# You should have received a copy of the GNU Moleculesral Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.

#' Derivative information for \code{lmmspline} objects
#' 
#' Calculates the derivative information for \code{lmmspline} objects with a \code{"p-spline"} or \code{"cubic p-spline"} basis.
#' 
#' @importFrom stats deriv
#' @param expr An object of class \code{lmmspline}.
#' @param ... Additional arguments which are passed to \code{deriv}.
#' @return deriv returns an object of class \code{lmmspline} containing the following components:
#' \item{predSpline}{ \code{data.frame} containing the predicted derivative values based on the linear model object or the linear mixed effect model object.}
#' \item{modelsUsed}{\code{numeric} vector indicating the model used to fit the data. 0 = linear model, 1 = linear mixed effect model spline (LMMS) with defined basis ("cubic" by default), 2 = LMMS taking subject-specific random intercept, 3 = LMMS with subject specific intercept and slope.}
#' \item{model}{\code{list} of models used to model time profiles.}
#' \item{derivative}{\code{logical} value indicating if the predicted values are the derivative information.}
#' @examples
#' \dontrun{
#' data(kidneySimTimeGroup)
#' # run lmmSpline on the samples from group 1 only
#' G1 <- which(kidneySimTimeGroup$group=="G1")
#' testLMMSplineTG<- lmmSpline(data=kidneySimTimeGroup$data[G1,],
#'                   time=kidneySimTimeGroup$time[G1],
#'                   sampleID=kidneySimTimeGroup$sampleID[G1],
#'                   basis="p-spline",keepModels=T)
#' testLMMSplineTGDeri <- deriv(testLMMSplineTG)
#' summary(testLMMSplineTGDeri)}

#' @export
deriv.lmmspline <- function(expr, ...){

models <- expr@models
if(length(models)==0)
  stop('You will need to keep the models to get the derivative information.')

basis <- expr@basis


if(sum(expr@basis%in%c('p-spline','cubic p-spline'))==0)
  stop('Objects must be modelled with p-spline or cubic p-spline basis.')


derivLme <- function(fit){ 
  #random slopes
  
  if(class(fit)=='lm'){
    beta.hat <- rep(fit$coefficients[2],length(unique(fit$model$time)))
    return(beta.hat)
    
  }else if(class(fit)=='lme'){
    u <- unlist(fit$coefficients$random$all)
    beta.hat <- fit$coefficients$fixed[2]
    Zt <-  fit$data$Zt[!duplicated(fit$data$Zt),]>0
    deriv.all <-    beta.hat + rowSums(Zt%*%t(u)) 
    return(deriv.all)
  }
}

#penalized cubic

derivLmeCubic <- function(fit){ 
  #random slopes
  if(class(fit)=='lm'){
    beta.hat <- rep(fit$coefficients[2],length(unique(fit$model$time)))
    return(beta.hat)
    
  }else if(class(fit)=='lme'){
    u <- unlist(fit$coefficients$random$all)
    beta.hat <- fit$coefficients$fixed[2]
    PZ <-  fit$data$Zt[!duplicated(fit$data$Zt),]
    PZ <-PZ^(1/3)
    deriv.all <-    beta.hat + rowSums((PZ*PZ)%*%(t(u)*3)) 
    return(deriv.all)
  }
  
}


deri <- NULL
new.data <- lapply(1:length(models),function(i){
  if(basis=='p-spline')
    deri <- derivLme(models[[i]])
  if(basis=='cubic p-spline')
    deri <- derivLmeCubic(models[[i]])
  return(deri)
})

if(class(models[[1]])=="lme"){
  time <- models[[1]]$data$time
}else{
  time <- models[[1]]$model$time
}


deri <- as.data.frame(matrix(unlist(new.data),nrow=length(models),ncol=(length(unlist(new.data))/length(models)),byrow=T))
rownames(deri) <- rownames(expr@predSpline)
colnames(deri) <- unique(sort(time))
expr@predSpline <- deri
expr@derivative <- TRUE
return(expr)
}

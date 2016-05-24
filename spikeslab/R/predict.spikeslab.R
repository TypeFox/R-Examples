####**********************************************************************
####**********************************************************************
####
####  SPIKE AND SLAB 1.1.5
####
####  Copyright 2013, Cleveland Clinic Foundation
####
####  This program is free software; you can redistribute it and/or
####  modify it under the terms of the GNU General Public License
####  as published by the Free Software Foundation; either version 2
####  of the License, or (at your option) any later version.
####
####  This program is distributed in the hope that it will be useful,
####  but WITHOUT ANY WARRANTY; without even the implied warranty of
####  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
####  GNU General Public License for more details.
####
####  You should have received a copy of the GNU General Public
####  License along with this program; if not, write to the Free
####  Software Foundation, Inc., 51 Franklin Street, Fifth Floor,
####  Boston, MA  02110-1301, USA.
####
####  ----------------------------------------------------------------
####  Written and Developed by:
####  ----------------------------------------------------------------
####    Hemant Ishwaran, Ph.D.
####    Director of Statistical Methodology
####    Professor, Division of Biostatistics
####    Clinical Research Building, Room 1058
####    1120 NW 14th Street
####    University of Miami, Miami FL 33136
####
####    email:  hemant.ishwaran@gmail.com
####    URL:    http://web.ccs.miami.edu/~hishwaran
####    --------------------------------------------------------------
####    J. Sunil Rao, Ph. D.
####    Professor and Director of the Division of Biostatistics, 
####    Department of Epidemiology & Public Health
####    Clinical Research Bldg, R-669
####    1120 NW 14th Street, Room 1056
####    Miami, FL 33136
####    email:  rao.jsunil@gmail.com
####    URL:    http://biostat.med.miami.edu/people/primary-faculty/sunil-rao
####  ----------------------------------------------------------------
####  Maintained by:
####    Udaya B. Kogalur, Ph.D.
####    Adjunct Staff
####    Dept of Quantitative Health Sciences
####    Cleveland Clinic Foundation
####    
####    Kogalur & Company, Inc.
####    5425 Nestleway Drive, Suite L1
####    Clemmons, NC 27012
####
####    email:  ubk@kogalur.com
####    URL:    http://www.kogalur.com
####    --------------------------------------------------------------
####
####**********************************************************************
####**********************************************************************

predict.spikeslab <- function(object, newdata = NULL, ...)
{

 ### --------------------------------------------------------------
 ### preliminary checks

  #check that object is compatible
  if (!inherits(object,"spikeslab"))
     stop("This function only works for objects of class `spikeslab'")
  
  #check whether object inherits cv type
  #make suitable alterations to merge mixing results
  if (sum(inherits(object, c("spikeslab", "cv"), TRUE) == c(1, 2)) == 2) {
    object <- object$spikeslab.obj
  }
  
  #check whether object inherits mixing type
  #make suitable alterations to merge mixing results
  do.mix <- FALSE
  yhat.mix <- NULL
  if (sum(inherits(object, c("spikeslab", "mixing"), TRUE) == c(1, 2)) == 2) {
     grr.mix.scale <- object$grr.mix.scale
     grr.mix.scale[is.na(grr.mix.scale)] <- 0
     object <- object$spikeslab.obj
     do.mix <- TRUE
  }
  
  #check for newdata
  if (is.null(newdata)) {
    newdata <- object$x
    colnames(newdata) <- object$names
  }
  #special allowance for newdata with one observation
  newdata <- rbind(newdata)
  #check for colnames
  if (is.null(colnames(newdata))) {
    if (ncol(newdata) == length(object$names)) {
       colnames(newdata) <- object$names
    } else {
       stop("Number of columns of testing data does not match training data")
    }
  }
  
 ### --------------------------------------------------------------
 ### ensure coherence between test and training data

 #restrict predictors in test data to those in the training data
 newdata <- as.data.frame(rbind(newdata))
 if (is.null(object$terms)) {
   old.xnames <- object$names
 }
 else {
   old.xnames <- attr(object$terms, "term.labels")
 }
 new.xnames <- names(newdata)
 if (sum(!is.element(old.xnames, new.xnames)) > 0)
   stop("predictors in test data don't match training data")
 newdata <- newdata[, match(old.xnames, new.xnames), drop = FALSE]
 
 #build the test design matrix
 if (is.null(object$terms)) {
   x.new <- newdata
 }
 else {
   Terms <- delete.response(object$terms)
   old.na.action <- options()$na.action
   na.keep <- function(x){x}
   options(na.action = na.keep)
   x.new <- model.matrix(Terms, newdata)
   options(na.action=old.na.action)
 }

 #ensure that test data matches training data one more time
 #needed because of factors
 cov.names <- unlist(dimnames(x.new)[2])
 x.new <- x.new[, match(object$names, cov.names), drop = FALSE]

 #remove NA's 
 x.new <- as.matrix(rbind(x.new))
 has.NA <- which(apply(x.new,1,function(x){any(is.na(x))}))
 if (length(has.NA) > 0){
   x.new <- as.matrix(x.new[-has.NA,, drop = FALSE])
 }
 rownames(x.new) <- colnames(x.new) <- NULL
  
 ### --------------------------------------------------------------
 ### extract x/y centering information
 ### get coefficient values
 ### correct for division by zero
 ### NA coefficient estimates ---> 0
  
 y.center <- object$y.center
 x.center <- object$x.center
 bma.scale <- object$bma.scale
 bma.scale[is.na(bma.scale)] <- 0
 gnet.scale <- object$gnet.scale
 gnet.scale[is.na(gnet.scale)] <- 0
 gnet.path <- object$gnet.path$path
 gnet.path[is.na(gnet.path)] <- 0

 ### --------------------------------------------------------------
 ### predicted values

 yhat.bma <- c(y.center - sum(bma.scale * x.center) + x.new %*% bma.scale)
 yhat.gnet <- c(y.center - sum(gnet.scale * x.center) + x.new %*% gnet.scale)
 yhat.gnet.path <- apply(gnet.path, 1, function(sbeta) {
                         c(y.center - sum(sbeta * x.center) + x.new %*% sbeta)})
 if (do.mix) {
   yhat.mix <- c(y.center - sum(grr.mix.scale * x.center) + x.new %*% grr.mix.scale)
 }

 ### --------------------------------------------------------------
 #### return the goodies

 return(list(yhat.bma = yhat.bma,
             yhat.gnet = yhat.gnet,
             yhat.mix = yhat.mix,
             yhat.gnet.path = yhat.gnet.path))

}

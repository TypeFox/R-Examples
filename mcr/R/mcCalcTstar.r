###############################################################################
##
## mcCalcTstar.R
##
## Function for computing of resampling T-Star statistic which we need
## to calculate tBoot confidence intervals for response and bias
## at some clinical decision point.
##
## Copyright (C) 2011 Roche Diagnostics GmbH
##
## This program is free software: you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## any later version.
##
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with this program.  If not, see <http://www.gnu.org/licenses/>.
##
###############################################################################

#' Compute Resampling T-statistic.
#' 
#' Compute Resampling T-statistic.
#' for Calculation of t-Bootstrap Confidence Intervals.
#'
#' @param .Object object of class "MCResultResampling".
#' @param x.levels a numeic vector of clinical desision points of interest.
#' @param iter.max maximal number of iterations  for calculation of weighted deming regression.
#' @param threshold threshold for calculation of weighted deming regression.
#' @return Tstar numeric vector containing resampling pivot statistic.
#' @references  Carpenter J., Bithell J. 
#'              Bootstrap confidence intervals: when, which, what? A practical guide for medical statisticians.
#'              Stat Med, 19 (9), 1141-1164 (2000).
mc.calcTstar <- function(.Object, x.levels, iter.max=30, threshold = 0.000001)
{
    stopifnot(is.numeric(x.levels))

    npoints <- length(.Object@data[,"x"])
    nsamples <- .Object@nsamples
    nnested <- .Object@nnested
    regmeth <- .Object@regmeth
    cimeth <- .Object@cimeth
    bootcimeth <- .Object@bootcimeth

    b0 <- .Object@glob.coef[1]
    b1 <- .Object@glob.coef[2]
    B0 <- .Object@B0
    B1<- .Object@B1
    sigmaB0 <- .Object@sigmaB0
    sigmaB1<- .Object@sigmaB1
    xw <- .Object@xmean
    MX <- .Object@MX

    # 1) response.glob  for global estimation of response at xc
    # 2) response       for resampling estimation of response at xc
    # 3) sigma.response for standard error of response estimation at xc
    #
    #  !! if we don't have any analytical estimate of standard error we
    #  need to use "nestedbootstrap" and obtain an resampling estimate at xc
    #
    # xc - medical desision point

    response.glob <- response <- sigma.response <- matrix(nrow=nsamples,ncol=length(x.levels))
    for (k in 1:nsamples) 
        response.glob[k,]<-b0 + b1*x.levels

    if(cimeth != "nestedbootstrap")
    {
        for (i in 1:nsamples) 
        {
            response[i,] <-B0[i] + B1[i]*x.levels
            Sigma2 <- sigmaB0[i]^2+sigmaB1[i]^2*x.levels*(x.levels-2*MX[i])
            if(any(Sigma2<0))
                stop("The calculation of CI for Response and Bias with tBoot method is collapsed.")
            else sigma.response[i,]<-sqrt(Sigma2)
        }
    } 
    else 
    {
        if(is.element(regmeth,c("Deming","WDeming"))) 
            error.ratio <- .Object@error.ratio

        data<-getData(.Object)
        X<-data[  ,"X"]
        Y<-data[  ,"Y"]

        if(regmeth=="WDeming" & (min(X)<0 | min(Y)<0)) stop("Weighted Deming regression for negative values is not available.")
      	if(regmeth=="PaBa" & (min(X)<0 | min(Y)<0)) stop("Passing-Bablok regression for negative values is not available.")

      	# Choice of arguments for regfun
      	if(regmeth == "LinReg")
      		callfun.reg <- function(idx){ mc.linreg(X[idx],Y[idx]) } 	
        else if(regmeth == "WLinReg")
            callfun.reg <- function(idx){ mc.wlinreg(X[idx],Y[idx]) } 
        else if (regmeth == "Deming")
      		callfun.reg <- function(idx){ mc.deming(X[idx],Y[idx],error.ratio=error.ratio) } 
        else if (regmeth == "WDeming")
      		callfun.reg <- function(idx) { mc.wdemingConstCV(X[idx],Y[idx],error.ratio=error.ratio,iter.max=iter.max,threshold=threshold) } 
        else if (regmeth == "PaBa") 
        {
      		## Determine if slope 1 or -1 is expected on complete data set
      		paba.posCor <- cor.test(X,Y,method="kendall")$estimate>=0
      		## Compute slope matrix once for all further computations
      		paba.angM <- mc.calcAngleMat(X,Y,posCor=paba.posCor)
      		## Regression function
      		callfun.reg<- function(idx)
            {
      			mc.res <- mc.paba(paba.angM[idx,idx],X[idx],Y[idx],posCor=paba.posCor,calcCI=FALSE)
      			return(list(b0=mc.res["Intercept","EST"],b1=mc.res["Slope","EST"],xw=as.numeric(NA)))
      		}
      	}     
        else 
            stop("Unknown regression function!")

        ## Bootstrap
        for (i in 1:nsamples) 
        {
            ## Draw bootstrap sample
            index <- sample(1:npoints,size=npoints,replace=TRUE)

            ## Regression on bootstrap sample
            d <- callfun.reg(index)
            ## Parameter point estimates
            response[i,]<-d$b0+d$b1*x.levels

      	    ## Parameter SE estimation

            ## Bootstrap
            response.nested <- matrix(nrow=nnested,ncol=length(x.levels))

          	for (ii in 1:nnested) 
            {
                ## Draw nested bootstrap sample
                nindex <- sample(1:npoints,size=npoints,replace=TRUE)
                ## Regression on nested bootstrap sample
                d <- callfun.reg(index[nindex])
                response.nested[ii,]<-d$b0+d$b1*x.levels
            } #end nest

            ## Estimate SD from nested bootstrap results
            sigma.response[i,]<-apply(response.nested,2,sd)
        }

    } # end else

    Tstar<-(response-response.glob)/sigma.response
    return(Tstar)
}
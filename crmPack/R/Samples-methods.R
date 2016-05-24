#####################################################################################
## Author: Daniel Sabanes Bove [sabanesd *a*t* roche *.* com],
##         Wai Yin Yeung [w *.* yeung1 *a*t* lancaster *.* ac *.* uk]
## Project: Object-oriented implementation of CRM designs
##
## Time-stamp: <[Samples-methods.R] by DSB Mon 11/05/2015 17:46>
##
## Description:
## Methods for processing the MCMC samples.
##
## History:
## 25/03/2014   file creation
## 10/07/2015   Adding more methods for pseudo models
#####################################################################################

##' @include McmcOptions-methods.R
##' @include Model-methods.R
##' @include fromQuantiles.R
{}

## --------------------------------------------------
## Extract certain parameter from "Samples" object to produce
## plots with "ggmcmc" package
## --------------------------------------------------

##' Get specific parameter samples and produce a data.frame
##'
##' Here you have to specify with \code{pos} which
##' parameter you would like to extract from the \code{\linkS4class{Samples}}
##' object
##'
##' @param x the \code{\linkS4class{Samples}} object
##' @param pos the name of the parameter
##' @param envir for vectorial parameters, you can give the indices of the
##' elements you would like to extract. If \code{NULL}, the whole vector samples
##' will be returned
##' @param mode not used
##' @param inherits not used
##'
##' @return the data frame suitable for use with \code{\link[ggmcmc]{ggmcmc}}
##'
##' @example examples/Sample-methods-get.R
##' @export
##' @keywords methods
setMethod("get",
          signature=
              signature(x="Samples",
                        pos="character",
                        envir="ANY",
                        mode="ANY",
                        inherits="ANY"),
          def=
          function(x,
                   pos,
                   envir=NULL,
                   mode=NULL,
                   inherits=NULL){

              ## check the parameter name
              stopifnot(is.scalar(pos),
                        pos %in% names(x@data))

              ## get the samples for this parameter
              d <- x@data[[pos]]
              ## this can be either a vector or a matrix

              ## how many parameters do we have?
              nPars <- NCOL(d)

              ## what are the names of all parameter
              ## elements?
              elements <-
                  if(nPars == 1L)
                      pos
                  else
                      paste(pos,
                            "[", seq_len(nPars), "]",
                            sep="")

              ## in case we have a vector parameter
              if(nPars > 1L)
              {
                  ## what are the indices to be returned?
                  indices <-
                      if(is.null(envir))
                      {
                          seq_along(elements)
                      } else {
                          stopifnot(is.numeric(envir),
                                    all(envir %in% seq_along(elements)))
                          as.integer(envir)
                      }

                  ## subset the data matrix and par names appropriately
                  d <- d[, indices, drop=FALSE]
                  elements <- elements[indices]

                  ## and also reduce the number of parameters
                  nPars <- length(indices)
              }

              ## now we can build
              ret <- data.frame(Iteration=seq_len(NROW(d)),
                                Chain=1L,
                                Parameter=
                                factor(rep(elements, each=NROW(d)),
                                       levels=elements),
                                value=as.numeric(d))

              ## add the attributes
              ret <- structure(ret,
                               nChains=1L,
                               nParameters=nPars,
                               nIterations=x@options@iterations,
                               nBurnin=x@options@burnin,
                               nThin=x@options@step,
                               description=elements,
                               parallel=FALSE)
              return(ret)
          })


## --------------------------------------------------
## Get fitted curves from Samples
## --------------------------------------------------

##' Fit method for the Samples class
##'
##' Note this new generic function is necessary because the \code{\link{fitted}}
##' function only allows the first argument \code{object} to appear in the
##' signature. But we need also other arguments in the signature.
##'
##' @param object the \code{\linkS4class{Samples}} object
##' @param model the \code{\linkS4class{Model}} object
##' @param data the \code{\linkS4class{Data}} object
##' @param \dots unused
##' @return the data frame with required information (see method details)
##'
##' @export
##' @keywords methods
setGeneric("fit",
           def=
           function(object,
                    model,
                    data,
                    ...){
               ## there should be no default method,
               ## therefore just forward to next method!
               standardGeneric("fit")},
           valueClass="data.frame")


## --------------------------------------------------
## Get fitted dose-tox curve from Samples
## --------------------------------------------------

##' @param points at which dose levels is the fit requested? default is the dose
##' grid
##' @param quantiles the quantiles to be calculated (default: 0.025 and
##' 0.975)
##' @param middle the function for computing the middle point. Default:
##' \code{\link{mean}}
##'
##' @describeIn fit This method returns a data frame with dose, middle, lower
##' and upper quantiles for the dose-toxicity curve
##' @example examples/Sample-methods-fit.R
##' 
setMethod("fit",
          signature=
          signature(object="Samples",
                    model="Model",
                    data="Data"),
          def=
          function(object,
                   model,
                   data,
                   points=data@doseGrid,
                   quantiles=c(0.025, 0.975),
                   middle=mean,
                   ...){
              ## some checks
              stopifnot(is.probRange(quantiles),
                        is.numeric(points))

              ## first we have to get samples from the dose-tox
              ## curve at the dose grid points.
              probSamples <- matrix(nrow=sampleSize(object@options),
                                    ncol=length(points))

              ## evaluate the probs, for all samples.
              for(i in seq_along(points))
              {
                  ## Now we want to evaluate for the
                  ## following dose:
                  probSamples[, i] <- prob(dose=points[i],
                                           model,
                                           object)
              }

              ## extract middle curve
              middleCurve <- apply(probSamples, 2L, FUN=middle)

              ## extract quantiles
              quantCurve <- apply(probSamples, 2L, quantile,
                                  prob=quantiles)

              ## now create the data frame
              ret <- data.frame(dose=points,
                                middle=middleCurve,
                                lower=quantCurve[1, ],
                                upper=quantCurve[2, ])

              ## return it
              return(ret)
          })

## --------------------------------------------------
## Get fitted dose-tox and dose-biomarker curves from Samples
## --------------------------------------------------

##' @describeIn fit This method returns a data frame with dose, and middle,
##' lower and upper quantiles, for both the dose-tox and dose-biomarker (suffix
##' "Biomarker") curves, for all grid points (Note that currently only the grid
##' points can be used, because the DualEndpointRW models only allow that)
##' 
##' @example examples/Sample-methods-fit-DualEndpoint.R
setMethod("fit",
          signature=
          signature(object="Samples",
                    model="DualEndpoint",
                    data="DataDual"),
          def=
          function(object,
                   model,
                   data,
                   quantiles=c(0.025, 0.975),
                   middle=mean,
                   ...){
              ## some checks
              stopifnot(is.probRange(quantiles))

              ## first obtain the dose-tox curve results from the parent method
              start <- callNextMethod(object=object,
                                      model=model,
                                      data=data,
                                      points=data@doseGrid,
                                      quantiles=quantiles,
                                      middle=middle,
                                      ...)

              ## now obtain the dose-biomarker results

              ## get the biomarker level samples
              ## at the dose grid points.
              biomLevelSamples <- matrix(nrow=sampleSize(object@options),
                                         ncol=data@nGrid)

              ## evaluate the biomLevels, for all samples.
              for(i in seq_len(data@nGrid))
              {
                  ## Now we want to evaluate for the
                  ## following dose:
                  biomLevelSamples[, i] <- biomLevel(dose=data@doseGrid[i],
                                                     xLevel=i,
                                                     model,
                                                     object)
              }

              ## extract middle curve
              middleCurve <- apply(biomLevelSamples, 2L, FUN=middle)

              ## extract quantiles
              quantCurve <- apply(biomLevelSamples, 2L, quantile,
                                  prob=quantiles)

              ## now create the data frame
              biomResults <- data.frame(middleBiomarker=middleCurve,
                                        lowerBiomarker=quantCurve[1, ],
                                        upperBiomarker=quantCurve[2, ])

              ## return both, pasted together
              return(cbind(start, biomResults))
          })


## --------------------------------------------------
## Approximate posterior with (log) normal distribution
## --------------------------------------------------

##' Approximate posterior with (log) normal distribution
##'
##' It is recommended to use \code{\link{set.seed}} before, in order
##' to be able to reproduce the resulting approximating model exactly.
##'
##' @param object the \code{\linkS4class{Samples}} object
##' @param model the \code{\linkS4class{Model}} object
##' @param data the \code{\linkS4class{Data}} object
##' @param \dots additional arguments (see methods)
##' @return the approximation model
##'
##' @export
##' @keywords methods
setGeneric("approximate",
           def=
           function(object, model, data, ...){
               ## there should be no default method,
               ## therefore just forward to next method!
               standardGeneric("approximate")},
           valueClass="Model")



##' @param points optional parameter, which gives the dose values at which
##' the approximation should rely on (default: 5 values equally spaced from
##' minimum to maximum of the dose grid)
##' @param refDose the reference dose to be used (default: median of
##' \code{points})
##' @param logNormal use the log-normal prior? (not default) otherwise, the
##' normal prior for the logistic regression coefficients is used
##' @param verbose be verbose (progress statements and plot)? (default)
##'
##' @describeIn approximate Here the \dots argument can transport additional arguments for
##' \code{\link{Quantiles2LogisticNormal}}, e.g. in order to control the
##' approximation quality, etc.
##' 
##' @example examples/Sample-methods-approximate.R
setMethod("approximate",
          signature=
          signature(object="Samples"),
          def=
          function(object,
                   model,
                   data,
                   points=
                   seq(from=min(data@doseGrid),
                       to=max(data@doseGrid),
                       length=5L),
                   refDose=median(points),
                   logNormal=FALSE,
                   verbose=TRUE,
                   ...){

              ## get the required quantiles at these dose levels:
              quants <- fit(object,
                            model,
                            data,
                            points=points,
                            quantiles=c(0.025, 0.975),
                            middle=median)

              ## get better starting values if it is already a logistic normal
              ## model
              if(is(model, "LogisticNormal") && (! logNormal))
              {
                  means <- sapply(object@data,
                                  mean)
                  cov <- cov(as.data.frame(object@data))

                  parstart <- c(means[1], means[2],
                                sqrt(cov[1, 1]), sqrt(cov[2, 2]),
                                cov2cor(cov)[1, 2])
              } else if(is(model, "LogisticLogNormal") && logNormal) {
                  datTrafo <- with(object@data,
                                   cbind(alpha0,
                                         log(alpha1)))

                  means <- colMeans(datTrafo)
                  cov <- cov(datTrafo)

                  parstart <- c(means[1], means[2],
                                sqrt(cov[1, 1]), sqrt(cov[2, 2]),
                                cov2cor(cov)[1, 2])
              } else {
                  parstart <- NULL
              }

              ## run the approx function
              quantRes <- Quantiles2LogisticNormal(dosegrid=quants$dose,
                                                   refDose=refDose,
                                                   lower=quants$lower,
                                                   upper=quants$upper,
                                                   median=quants$middle,
                                                   verbose=verbose,
                                                   parstart=parstart,
                                                   logNormal=logNormal,
                                                   ...)

              if(verbose)
              {
                  matplot(x=points,
                          quantRes$required,
                          type="l", col="blue", lty=1)
                  matlines(x=points,
                           quantRes$quantiles,
                           col="red", lty=1)
                  legend("bottomright",
                         legend=c("original", "approximation"),
                         col=c("blue", "red"),
                         lty=1,
                         bty="n")
              }

              ## return the model
              return(quantRes$model)
          })

## --------------------------------------------------
## Plot dose-tox fit from a model
## --------------------------------------------------


##' Plotting dose-toxicity model fits
##'
##' @param x the \code{\linkS4class{Samples}} object
##' @param y the \code{\linkS4class{Model}} object
##' @param data the \code{\linkS4class{Data}} object
##' @param xlab the x axis label
##' @param ylab the y axis label
##' @param showLegend should the legend be shown? (default)
##' @param \dots not used
##' @return This returns the \code{\link[ggplot2]{ggplot}}
##' object for the dose-toxicity model fit
##'
##' @example examples/Sample-methods-plot.R
##' @export
##' @importFrom ggplot2 qplot scale_linetype_manual
setMethod("plot",
          signature=
          signature(x="Samples",
                    y="Model"),
          def=
          function(x, y, data, ...,
                   xlab="Dose level",
                   ylab="Probability of DLT [%]",
                   showLegend=TRUE){

              ## check args
              stopifnot(is.bool(showLegend))

              ## get the fit
              plotData <- fit(x,
                              model=y,
                              data=data,
                              quantiles=c(0.025, 0.975),
                              middle=mean)

              ## make the plot
              gdata <-
                  with(plotData,
                       data.frame(x=rep(dose, 3),
                                  y=c(middle, lower, upper) * 100,
                                  group=
                                  rep(c("mean", "lower", "upper"),
                                      each=nrow(plotData)),
                                  Type=
                                  factor(c(rep("Estimate",
                                               nrow(plotData)),
                                           rep("95% Credible Interval",
                                               nrow(plotData) * 2)),
                                         levels=
                                         c("Estimate",
                                           "95% Credible Interval"))))

              ret <- ggplot2::qplot(x=x,
                                    y=y,
                                    data=gdata,
                                    group=group,
                                    linetype=Type,
                                    colour=I("red"),
                                    geom="line",
                                    xlab=xlab,
                                    ylab=ylab,
                                    ylim=c(0, 100))

              ret <- ret +
                  ggplot2::scale_linetype_manual(breaks=
                                                 c("Estimate",
                                                   "95% Credible Interval"),
                                                 values=c(1,2), guide=ifelse(showLegend,
                                                 "legend", FALSE))

              return(ret)
          })


## --------------------------------------------------
## Special method for dual endpoint model
## --------------------------------------------------


##' Plotting dose-toxicity and dose-biomarker model fits
##'
##' When we have the dual endpoint model,
##' also the dose-biomarker fit is shown in the plot
##'
##' @param x the \code{\linkS4class{Samples}} object
##' @param y the \code{\linkS4class{DualEndpoint}} object
##' @param data the \code{\linkS4class{DataDual}} object
##' @param extrapolate should the biomarker fit be extrapolated to the whole
##' dose grid? (default)
##' @param showLegend should the legend be shown? (not default)
##' @param \dots additional arguments for the parent method
##' \code{\link{plot,Samples,Model-method}}
##' @return This returns the \code{\link[ggplot2]{ggplot}}
##' object with the dose-toxicity and dose-biomarker model fits
##'
##' @example examples/Sample-methods-plot-DualEndpoint.R
##' @export
setMethod("plot",
          signature=
          signature(x="Samples",
                    y="DualEndpoint"),
          def=
          function(x, y, data, extrapolate=TRUE, showLegend=FALSE, ...){

              stopifnot(is.bool(extrapolate))

              ## call the superclass method, to get the toxicity plot
              plot1 <- callNextMethod(x, y, data, showLegend=showLegend, ...)

              ## only look at these dose levels for the plot:
              xLevels <-
                  if(extrapolate)
                      seq_along(data@doseGrid)
                  else
                      1:max(data@xLevel)

              ## get the plot data for the biomarker plot
              functionSamples <- matrix(nrow=sampleSize(x@options),
                                        ncol=length(xLevels))

              ## evaluate the biomLevels, for all samples.
              for(i in seq_along(xLevels))
              {
                  ## Now we want to evaluate for the
                  ## following dose:
                  functionSamples[, i] <-
                      biomLevel(dose=data@doseGrid[xLevels[i]],
                                xLevel=xLevels[i],
                                model=y,
                                samples=x)
              }

              ## extract mean curve
              meanCurve <- colMeans(functionSamples)

              ## extract quantiles
              quantiles <- c(0.025, 0.975)
              quantCurve <- apply(functionSamples, 2L, quantile,
                                  prob=quantiles)

              ## now create the data frame
              plotData <- data.frame(dose=data@doseGrid[xLevels],
                                     mean=meanCurve,
                                     lower=quantCurve[1, ],
                                     upper=quantCurve[2, ])

              ## make the second plot
              gdata <-
                  with(plotData,
                       data.frame(x=rep(dose, 3),
                                  y=c(mean, lower, upper),
                                  group=
                                  rep(c("mean", "lower", "upper"),
                                      each=nrow(plotData)),
                                  Type=
                                  factor(c(rep("Estimate",
                                               nrow(plotData)),
                                           rep("95% Credible Interval",
                                               nrow(plotData) * 2)),
                                         levels=
                                         c("Estimate",
                                           "95% Credible Interval"))))

              plot2 <- ggplot2::qplot(x=x,
                                      y=y,
                                      data=gdata,
                                      group=group,
                                      linetype=Type,
                                      colour=I("blue"),
                                      geom="line",
                                      xlab="Dose level",
                                      ylab="Biomarker level")

              plot2 <- plot2 +
                  ggplot2::scale_linetype_manual(breaks=
                                                 c("Estimate",
                                                   "95% Credible Interval"),
                                                 values=c(1,2),
                                                 guide=ifelse(showLegend,
                                                 "legend", FALSE))

              ## arrange both plots side by side
              ret <- gridExtra::arrangeGrob(plot1, plot2, ncol=2)
              return(ret)
          })


## -------------------------------------------------------------------------------------
## Get fitted dose-tox curve from Samples for 'LogisticIndepBeta' model class
## ------------------------------------------------------------------------------------
##' @describeIn fit This method return a data frame with dose, middle lower and upper quantiles 
##' for the dose-DLE curve using DLE samples for \dQuote{LogisticIndepBeta} model class
##' @example examples/Samples-method-fitDLE.R
setMethod("fit",
          signature=
            signature(object="Samples",
                      model="LogisticIndepBeta",
                      data="Data"),
          def=
            function(object,
                     model,
                     data,
                     points=data@doseGrid,
                     quantiles=c(0.025, 0.975),
                     middle=mean,
                     ...){
              ## some checks
              stopifnot(is.probRange(quantiles),
                        is.numeric(points))
            
              ## first we have to get samples from the dose-tox
              ## curve at the dose grid points.
              probSamples <- matrix(nrow=sampleSize(object@options),
                                    ncol=length(points))
              
              ## evaluate the probs, for all samples.
              for(i in seq_along(points))
              {
                ## Now we want to evaluate for the
                ## following dose:
                probSamples[, i] <- prob(dose=points[i],
                                         model,
                                         object)
              }
              
              ## extract middle curve
              middleCurve <- apply(probSamples, 2L, FUN=middle)
              
              ## extract quantiles
              quantCurve <- apply(probSamples, 2L, quantile,
                                  prob=quantiles)
              
              ## now create the data frame
              ret <- data.frame(dose=points,
                                middle=middleCurve,
                                lower=quantCurve[1, ],
                                upper=quantCurve[2, ])
              
              ## return it
              return(ret)
            })

## -------------------------------------------------------------------------------------
## Get fitted dose-efficacy curve from Samples for 'Effloglog' model class
## ------------------------------------------------------------------------------------

##' @describeIn fit This method returns a data frame with dose, middle, lower, upper quantiles for 
##' the dose-efficacy curve using efficacy samples for \dQuote{Effloglog} model class
##' @example examples/Samples-method-fitEff.R
setMethod("fit",
          signature=
            signature(object="Samples",
                      model="Effloglog",
                      data="DataDual"),
          def=
            function(object,
                     model,
                     data,
                     points=data@doseGrid,
                     quantiles=c(0.025, 0.975),
                     middle=mean,
                     ...){
              ## some checks
              stopifnot(is.probRange(quantiles),
                        is.numeric(points))
              
              ## first we have to get samples from the dose-tox
              ## curve at the dose grid points.
              ExpEffSamples <- matrix(nrow=sampleSize(object@options),
                                      ncol=length(points))
              
              ## evaluate the probs, for all samples.
              for(i in seq_along(points))
              {
                ## Now we want to evaluate for the
                ## following dose:
                ExpEffSamples[, i] <- ExpEff(dose=points[i],
                                             model,
                                             object)
              }
              
              ## extract middle curve
              middleCurve <- apply(ExpEffSamples, 2L, FUN=middle)
              
              ## extract quantiles
              quantCurve <- apply(ExpEffSamples, 2L, quantile,
                                  prob=quantiles)
              
              ## now create the data frame
              ret <- data.frame(dose=points,
                                middle=middleCurve,
                                lower=quantCurve[1, ],
                                upper=quantCurve[2, ])
              
              ## return it
              return(ret)
            })
## ==========================================================================================
## --------------------------------------------------------------------
## Get fitted dose-efficacy based on the Efficacy Flexible model
## -------------------------------------------------------------
##' @describeIn fit This method returns a data frame with dose, middle, lower and upper 
##' quantiles for the dose-efficacy curve using efficacy samples for \dQuote{EffFlexi} 
##' model class
##' @example examples/Samples-method-fitEffFlexi.R
setMethod("fit",
          signature=
            signature(object="Samples",
                      model="EffFlexi",
                      data="DataDual"),
          def=
            function(object,
                     model,
                     data,
                     points=data@doseGrid,
                     quantiles=c(0.025, 0.975),
                     middle=mean,
                     ...){
              ## some checks
              stopifnot(is.probRange(quantiles),
                        is.numeric(points))
              
              ## first we have to get samples from the dose-tox
              ## curve at the dose grid points.
              ExpEffSamples <- matrix(nrow=sampleSize(object@options),
                                      ncol=length(points))
              
              ## evaluate the probs, for all samples.
              for(i in seq_along(points))
              {
                ## Now we want to evaluate for the
                ## following dose:
                ExpEffSamples[, i] <- ExpEff(dose=points[i],
                                             model,
                                             object)
              }
              
              ## extract middle curve
              middleCurve <- apply(ExpEffSamples, 2L, FUN=middle)
              
              ## extract quantiles
              quantCurve <- apply(ExpEffSamples, 2L, quantile,
                                  prob=quantiles)
              
              ## now create the data frame
              ret <- data.frame(dose=points,
                                middle=middleCurve,
                                lower=quantCurve[1, ],
                                upper=quantCurve[2, ])
              
              ## return it
              return(ret)
            })
## ==============================================================
## ----------------------------------------------------------------
## Get fitted values at all dose levels from gain samples 
## -----------------------------------------------------------------
##' Get the fiited values for the gain values at all dose levels based on 
##' a given pseudo DLE model, DLE sample, a pseudo efficacy model, a Efficacy sample 
##' and data. This method returns a data frame with dose, middle, lower and upper quantiles 
##' of the gain value samples
##' 
##' @param DLEmodel the DLE pseudo model of \code{\linkS4class{ModelTox}} class object
##' @param DLEsamples the DLE samples of \code{\linkS4class{Samples}} class object
##' @param Effmodel the efficacy pseudo model of \code{\linkS4class{ModelEff}} class object
##' @param Effsamples the efficacy samples of \code{\linkS4class{Samples}} class object
##' @param data the data input of \code{\linkS4class{DataDual}} class object
##' @param \dots additional arguments for methods
##' 
##' @export
##' @keywords methods
setGeneric("fitGain",
           def=
             function(DLEmodel,
                      DLEsamples,
                      Effmodel,
                      Effsamples,
                      data,
                      ...){
               ## there should be no default method,
               ## therefore just forward to next method!
               standardGeneric("fitGain")},
           valueClass="data.frame")

##' @describeIn fitGain This method returns a data frame with dose, middle, lower, upper quantiles for 
##' the gain values obtained given the DLE and the efficacy samples
##' @param points at which dose levels is the fit requested? default is the dose
##' grid
##' @param quantiles the quantiles to be calculated (default: 0.025 and
##' 0.975)
##' @param middle the function for computing the middle point. Default:
##' \code{\link{mean}}
##' @example examples/Samples-method-fitGain.R
setMethod("fitGain",
          signature=
            signature(DLEmodel="ModelTox",
                      DLEsamples="Samples",
                      Effmodel="ModelEff",
                      Effsamples="Samples",
                      data="DataDual"),
          def=
            function(DLEmodel,
                     DLEsamples,
                     Effmodel,
                     Effsamples,
                     data,
                     points=data@doseGrid,
                     quantiles=c(0.025, 0.975),
                     middle=mean,
                     ...){
              ## some checks
              stopifnot(is.probRange(quantiles),
                        is.numeric(points))
              
              ## first we have to get samples from the gain
              ## at the dose grid points.
              GainSamples <- matrix(nrow=sampleSize(DLEsamples@options),
                                    ncol=length(points))
              
              ## evaluate the probs, for all gain samples.
              for(i in seq_along(points))
              {
                ## Now we want to evaluate for the
                ## following dose:
                GainSamples[, i] <- gain(dose=points[i],
                                         DLEmodel=DLEmodel,
                                         DLEsamples=DLEsamples, 
                                         Effmodel=Effmodel,
                                         Effsamples=Effsamples)
              }
              
              ## extract middle curve
              middleCurve <- apply(GainSamples, 2L, FUN=middle)
              
              ## extract quantiles
              quantCurve <- apply(GainSamples, 2L, quantile,
                                  prob=quantiles)
              
              ## now create the data frame
              ret <- data.frame(dose=points,
                                middle=middleCurve,
                                lower=quantCurve[1, ],
                                upper=quantCurve[2, ])
              
              ## return it
              return(ret)
            })
## ---------------------------------------------------------------------------------
## Plot the fitted dose-DLE curve with pseudo DLE model with samples
## -------------------------------------------------------------------------------
##' Plot the fitted dose-DLE curve using a \code{\linkS4class{ModelTox}} class model with samples
##' 
##' @param x the \code{\linkS4class{Samples}} object
##' @param y the \code{\linkS4class{ModelTox}} model class object
##' @param data the \code{\linkS4class{Data}} object
##' @param xlab the x axis label
##' @param ylab the y axis label
##' @param showLegend should the legend be shown? (default)
##' @param \dots not used
##' @return This returns the \code{\link[ggplot2]{ggplot}}
##' object for the dose-DLE model fit
##' 
##' @example examples/Samples-method-plotModelTox.R
##' @export
##' @keywords methods
##' @importFrom ggplot2 qplot scale_linetype_manual
setMethod("plot",
          signature=
            signature(x="Samples",
                      y="ModelTox"),
          def=
            function(x, y, data, ...,
                     xlab="Dose level",
                     ylab="Probability of DLT [%]",
                     showLegend=TRUE){
              
              ## check args
              stopifnot(is.bool(showLegend))
              
         
              ## get the fit
              plotData <- fit(x,
                              model=y,
                              data=data,
                              quantiles=c(0.025, 0.975),
                              middle=mean)
              
              ## make the plot
              gdata <-
                with(plotData,
                     data.frame(x=rep(dose, 3),
                                y=c(middle, lower, upper) * 100,
                                group=
                                  rep(c("mean", "lower", "upper"),
                                      each=nrow(plotData)),
                                Type=
                                  factor(c(rep("Estimate",
                                               nrow(plotData)),
                                           rep("95% Credible Interval",
                                               nrow(plotData) * 2)),
                                         levels=
                                           c("Estimate",
                                             "95% Credible Interval"))))
              
              ret <- ggplot2::qplot(x=x,
                                    y=y,
                                    data=gdata,
                                    group=group,
                                    linetype=Type,
                                    colour=I("red"),
                                    geom="line",
                                    xlab=xlab,
                                    ylab=ylab,
                                    ylim=c(0, 100))
              
              ret <- ret +
                ggplot2::scale_linetype_manual(breaks=
                                                 c("Estimate",
                                                   "95% Credible Interval"),
                                               values=c(1,2), guide=ifelse(showLegend,
                                                                           "legend", FALSE))
              
              return(ret)
            })


## --------------------------------------------------------------------------------------------
## Plot the fitted dose-efficacy curve using a pseudo efficacy model with samples
## -------------------------------------------------------------------------------------------
##' Plot the fitted dose-effcacy curve using a model from \code{\linkS4class{ModelEff}} class
##' with samples
##' 
##' @param x the \code{\linkS4class{Samples}} object
##' @param y the \code{\linkS4class{ModelEff}} model class object
##' @param data the \code{\linkS4class{Data}} object
##' @param xlab the x axis label
##' @param ylab the y axis label
##' @param showLegend should the legend be shown? (default)
##' @param \dots not used
##' @return This returns the \code{\link[ggplot2]{ggplot}}
##' object for the dose-efficacy model fit
##' 
##' @example examples/Samples-method-plotModelEff.R 
##' @export
##' @keywords methods
##' @importFrom ggplot2 qplot scale_linetype_manual
setMethod("plot",
          signature=
            signature(x="Samples",
                      y="ModelEff"),
          def=
            function(x, y, data, ...,
                     xlab="Dose level",
                     ylab="Expected Efficacy",
                     showLegend=TRUE){
              
              ## check args
              stopifnot(is.bool(showLegend))
              
              ## get the fit
              plotData <- fit(x,
                              model=y,
                              data=data,
                              quantiles=c(0.025, 0.975),
                              middle=mean)
              
              ## make the plot
              gdata <-
                with(plotData,
                     data.frame(x=rep(dose, 3),
                                y=c(middle, lower, upper) ,
                                group=
                                  rep(c("mean", "lower", "upper"),
                                      each=nrow(plotData)),
                                Type=
                                  factor(c(rep("Estimate",
                                               nrow(plotData)),
                                           rep("95% Credible Interval",
                                               nrow(plotData) * 2)),
                                         levels=
                                           c("Estimate",
                                             "95% Credible Interval"))))
              
              ret <- ggplot2::qplot(x=x,
                                    y=y,
                                    data=gdata,
                                    group=group,
                                    linetype=Type,
                                    colour=I("blue"),
                                    geom="line",
                                    xlab=xlab,
                                    ylab=ylab,
                                    xlim=c(0,max(data@doseGrid)))
              
              ret <- ret +
                ggplot2::scale_linetype_manual(breaks=
                                                 c("Estimate",
                                                   "95% Credible Interval"),
                                               values=c(1,2), guide=ifelse(showLegend,
                                                                           "legend", FALSE))
              
              return(ret)
            })

## ----------------------------------------------------------------------------------------
## Plot of fitted dose-DLE curve based on a pseudo DLE model without sample 
##-------------------------------------------------------------------------------------
##' Plot of the fitted dose-tox based with a given pseudo DLE model and data without samples
##' 
##' @param x the data of \code{\linkS4class{Data}} class object
##' @param y the model of the \code{\linkS4class{ModelTox}} class object
##' @param xlab the x axis label
##' @param ylab the y axis label
##' @param showLegend should the legend be shown? (default)
##' @param \dots not used
##' @return This returns the \code{\link[ggplot2]{ggplot}}
##' object for the dose-DLE model plot
##' 
##' @example examples/Samples-method-plotModelToxNoSamples.R 
##' @export
##' @keywords methods
##' @importFrom ggplot2 qplot scale_linetype_manual
setMethod("plot",
          signature=
            signature(x="Data",
                      y="ModelTox"),
          def=
            function(x,y,
                     xlab="Dose level",
                     ylab="Probability of DLE",
                     showLegend=TRUE,...){
              ##check args
              
              stopifnot(is.bool(showLegend))
              
              ##Make sure the right model estimates are use with the given data
              y <- update(object=y,data=x)
              
              
              ##create data frame
              
              plotData <- data.frame(dose=x@doseGrid,
                                     probDLE=prob(dose=x@doseGrid,
                                                  model=y))
              ##Look for TD30 and TD35
              TD30 <-dose(prob=0.30,
                          model=y)
              TD35 <-dose(prob=0.35,
                          model=y)
              
              ##make the plot
              gdata <- with(plotData,
                            data.frame(x=dose,
                                       y=probDLE,
                                       group=rep("Estimated DLE",each=nrow(plotData)),
                                       Type=factor(rep("Estimated DLE",nrow(plotData)),levels="Estimated DLE")))
              
              plot1 <- ggplot2::qplot(x=x,
                                      y=y,
                                      data=gdata,
                                      group=group,
                                      linetype=Type,
                                      colour=I("red"),
                                      geom="line",
                                      xlab=xlab,
                                      ylab=ylab,
                                      ylim=c(0,1))
              
              plot1 <- plot1 + ggplot2::scale_linetype_manual(breaks="Estimated DLE",
                                                              values=c(1,2),
                                                              guide=ifelse(showLegend,"legend",FALSE))
              
              
              plot1 <- plot1 +
                geom_line(size=1.5,colour="red")
              
              return(plot1)
            })


## ---------------------------------------------------------------------------------------------
## Plot the fitted dose-efficacy curve given a pseudo efficacy model without samples
## ----------------------------------------------------------------------------------
##' Plot of the fitted dose-efficacy based with a given pseudo efficacy model and data without samples
##' 
##' @param x the data of \code{\linkS4class{DataDual}} class object
##' @param y the model of the \code{\linkS4class{ModelEff}} class object
##' @param xlab the x axis label
##' @param ylab the y axis label
##' @param showLegend should the legend be shown? (default)
##' @param \dots not used
##' @return This returns the \code{\link[ggplot2]{ggplot}}
##' object for the dose-efficacy model plot
##' 
##' @example examples/Samples-method-plotModelEffNoSamples.R 
##' @export
##' @keywords methods
##' @importFrom ggplot2 qplot scale_linetype_manual 

setMethod("plot",
          signature=
            signature(x="DataDual",
                      y="ModelEff"),
          def=
            function(x,y,...,
                     xlab="Dose level",
                     ylab="Expected Efficacy",
                     showLegend=TRUE){
              ##check args
              
              stopifnot(is.bool(showLegend))
              y <- update(object=y,data=x)
              
              ##create data frame
              
              plotEffData<- data.frame(dose=x@doseGrid,
                                       ExpEff=ExpEff(dose=x@doseGrid,
                                                     model=y))
              
              ##make the second plot
              ggdata<-with(plotEffData,
                           data.frame(x=dose,
                                      y=ExpEff,
                                      group=rep("Estimated Expected Efficacy",each=nrow(plotEffData)),
                                      Type=factor(rep("Estimated Expected Efficacy",nrow(plotEffData)),levels="Estimated Expected Efficacy")))
              
              ##Get efficacy plot
              plot2 <- ggplot(data=ggdata, aes(x=x,y=y), group=group) +
                xlab("Dose Levels")+
                ylab(paste("Estimated Expected Efficacy")) + xlim(c(0,max(x@doseGrid))) +
                geom_line(colour=I("blue"), size=1.5)
              
              plot2 <- plot2 +
                geom_line(size=1.5,colour="blue")
              
              
              return(plot2)
            })

## ----------------------------------------------------------------------------------------------------------
## Plot the gain curve using a pseudo DLE and a pseudo Efficacy model with samples
## ----------------------------------------------------------------------------------------------------
##' Plot the gain curve in addition with the dose-DLE and dose-efficacy curve using a given DLE pseudo model,
##' a DLE sample, a given efficacy pseudo model and an efficacy sample
##' 
##' @param DLEmodel the dose-DLE model of \code{\linkS4class{ModelTox}} class object
##' @param DLEsamples the DLE sample of \code{\linkS4class{Samples}} class object
##' @param Effmodel the dose-efficacy model of \code{\linkS4class{ModelEff}} class object
##' @param Effsamples the efficacy sample of of \code{\linkS4class{Samples}} class object
##' @param data the data input of \code{\linkS4class{DataDual}} class object
##' @param \dots not used
##' @return This returns the \code{\link[ggplot2]{ggplot}}
##' object for the plot
##' 
##' @example examples/Samples-method-plotGain.R
##' @export
##' @keywords methods
setGeneric("plotGain",
           def=
             function(DLEmodel,
                      DLEsamples,
                      Effmodel,
                      Effsamples,
                      data,...){
               standardGeneric("plotGain")})
##' @describeIn plotGain Standard method
setMethod("plotGain",
          signature=
            signature(DLEmodel="ModelTox",
                      DLEsamples="Samples",
                      Effmodel="ModelEff",
                      Effsamples="Samples"),
          def=
            function(DLEmodel,DLEsamples,Effmodel,Effsamples,data,...){
              
              ##Get fitted values for probabilities of DLE at all dose levels 
              plotDLEData <- fit(DLEsamples,
                                 model=DLEmodel,
                                 data=data,
                                 quantiles=c(0.025, 0.975),
                                 middle=mean)
              
              
              
              ##Get fitted values for mean efficacy values at all dose levels 
              plotEffData <- fit(Effsamples, 
                                 model=Effmodel,
                                 data=data,
                                 quantiles=c(0.025, 0.975),
                                 middle=mean)
              
              ##Get fitted values for gain values at all dose levels 
              plotGainData <- fitGain(DLEmodel=DLEmodel,
                                      DLEsamples=DLEsamples,
                                      Effmodel=Effmodel,
                                      Effsamples=Effsamples,
                                      data=data)
              
              ##For each of the dose levels, take the mean for the probabilties of DLE, mean efiicacy values 
              ## and gain values. Hence combine them into a data frame
              
              plotData<-data.frame(dose=rep(data@doseGrid,3),
                                   values=c(plotDLEData$middle,
                                            plotEffData$middle,
                                            plotGainData$middle))
              ## only the line plots for the mean value of the DLE, efficacy and gain samples 
              ##at all dose levels
              gdata<-with(plotData,
                          data.frame(x=dose,
                                     y=values,
                                     group=c(rep("p(DLE)",length(data@doseGrid)),
                                             rep("Mean Expected Efficacy",length(data@doseGrid)),
                                             rep("Gain",length(data@doseGrid))),
                                     Type=factor("Estimate",levels="Estimate")
                                     
                          ))
              
              plot1 <- ggplot(data=gdata, aes(x=x,y=y))+geom_line(aes(group=group,color=group),size=1.5)+
                ggplot2::scale_colour_manual(name="curves",values=c("green3","blue","red"))+
                xlab("Dose Level")+ xlim(c(0,max(data@doseGrid)))+
                ylab(paste("Values")) + ylim(c(min(gdata$y),max(gdata$y)))
              
              
              
              return(plot1)
            })

##----------------------------------------------------------------------------------------------------
## Plot the gain curve using a pseudo DLE and a pseudo Efficacy model without samples
## ----------------------------------------------------------------------------------------------------
##' Plot the gain curve in addition with the dose-DLE and dose-efficacy curve using a given DLE pseudo model,
##' and a given efficacy pseudo model
##' 
##' @describeIn plotGain Standard method
##' 
##' @example examples/Samples-method-plotGainNoSamples.R
##' @export
##' @keywords methods
setMethod("plotGain",
          signature=
            signature(DLEmodel="ModelTox",
                      DLEsamples="missing",
                      Effmodel="ModelEff",
                      Effsamples="missing"),
          def=
            function(DLEmodel,Effmodel,data,...){
            
              ##Make sure the model estimates are corresponds to the input data
              DLEmodel <- update(object=DLEmodel,data=data)
              Effmodel <- update(object=Effmodel,data=data)
              
              plotData<-data.frame(dose=rep(data@doseGrid,3),
                                   values=c(prob(dose=data@doseGrid,
                                                 model=DLEmodel),
                                            ExpEff(dose=data@doseGrid,
                                                   model=Effmodel),
                                            gain(dose=data@doseGrid,
                                                 DLEmodel=DLEmodel,
                                                 Effmodel=Effmodel)))
              gdata<-with(plotData,
                          data.frame(x=dose,
                                     y=values,
                                     group=c(rep("p(DLE)",length(data@doseGrid)),
                                             rep("Expected Efficacy",length(data@doseGrid)),
                                             rep("Gain",length(data@doseGrid))),
                                     Type=factor("Estimate",levels="Estimate")
                                     
                          ))
              
              ##plot1 <- ggplot(data=gdata, aes(x=x,y=y))+geom_line(aes(group=group,color=group),size=1.5)
              
              plot1 <- ggplot(data=gdata, aes(x=x,y=y))+geom_line(aes(group=group,color=group),size=1.5)+
                ggplot2::scale_colour_manual(name="curves",values=c("blue","green3","red"))+
                xlab("Dose Level")+ xlim(c(0,max(data@doseGrid)))+
                ylab(paste("Values")) + ylim(c(min(gdata$y),max(gdata$y)))
              
              
              
              TD30 <- dose(prob=0.3,model=DLEmodel)
              
              Gainfun<-function(DOSE){
                -gain(DOSE,DLEmodel=DLEmodel,Effmodel=Effmodel)
              }
              Gstar<-(optim(min(data@doseGrid),Gainfun,method = "L-BFGS-B",lower=min(data@doseGrid),upper=max(data@doseGrid))$par)
              MaxGain<--(optim(min(data@doseGrid),Gainfun,method = "L-BFGS-B",lower=min(data@doseGrid),upper=max(data@doseGrid))$value)
              
              
              if ((TD30 < min(data@doseGrid))|(TD30 > max(data@doseGrid))) {
                plot1<-plot1
                print(paste("TD30",paste(TD30," not within dose Grid")))} else {plot1 <-plot1 + geom_point(data=data.frame(x=TD30,y=0.3),aes(x=x,y=y),colour="violet", shape=16, size=8) +
                  annotate("text",label="p(DLE=0.3)",x=TD30+1,y=0.2,size=5,colour="violet")}
              
              
              
              if ((Gstar < min(data@doseGrid))|(Gstar > max(data@doseGrid))) {
                plot1<-plot1
                print(paste("Gstar=",paste(Gstar," not within dose Grid")))} else {plot1 <- plot1 + geom_point(data=data.frame(x=Gstar,y=MaxGain),aes(x=x,y=y),colour="green3", shape=17, size=8) +
                  annotate("text",label="Max Gain",x=Gstar,y=MaxGain-0.1,size=5,colour="green3")}
              
              return(plot1)
            })
##==========================================================================================

## -------------------------------------------------------------------------------
## Plot of the DLE and efficacy curve sides by side with samples
## -----------------------------------------------------------------------------
##' Plot of the DLE and efficacy curve side by side given a DLE pseudo model,
##' a DLE sample, an efficacy pseudo model and a given efficacy sample
##' 
##' @param DLEmodel the pseudo DLE model of \code{\linkS4class{ModelTox}} class object
##' @param DLEsamples the DLE samples of \code{\linkS4class{Samples}} class object
##' @param Effmodel the pseudo efficacy model of \code{\linkS4class{ModelEff}} class object
##' @param Effsamples the Efficacy samples of \code{\linkS4class{Samples}} class object
##' @param data the data input of \code{\linkS4class{DataDual}} class object
##' @param extrapolate should the biomarker fit be extrapolated to the whole
##' dose grid? (default)
##' @param showLegend should the legend be shown? (not default)
##' @param \dots additional arguments for the parent method
##' \code{\link{plot,Samples,Model-method}}
##' @return This returns the \code{\link[ggplot2]{ggplot}}
##' object with the dose-toxicity and dose-efficacy model fits
##' 
##' @example examples/Samples-method-plotDualResponses.R
##' 
##' @export
##' @keywords methods
setGeneric("plotDualResponses",
           def=
             function(DLEmodel,
                      DLEsamples,
                      Effmodel,
                      Effsamples,
                      data,...){
               standardGeneric("plotDualResponses")})

##' @describeIn plotDualResponses function todo
setMethod("plotDualResponses",
          signature=
            signature(DLEmodel="ModelTox",
                      DLEsamples="Samples",
                      Effmodel="ModelEff",
                      Effsamples="Samples"),
          def=
            function(DLEmodel, DLEsamples, Effmodel,Effsamples,data, extrapolate=TRUE, showLegend=FALSE,...){
              
              stopifnot(is.bool(extrapolate))
              ## Get Toxicity plot
              ## get the fit
              
              plotDLEData <- fit(DLEsamples,
                                 model=DLEmodel,
                                 data=data,
                                 quantiles=c(0.025, 0.975),
                                 middle=mean)
              
              ## make the plot
              gdata <-
                with(plotDLEData,
                     data.frame(x=rep(dose, 3),
                                y=c(middle, lower, upper) * 100,
                                group=
                                  rep(c("mean", "lower", "upper"),
                                      each=nrow(plotDLEData)),
                                Type=
                                  factor(c(rep("Estimate",
                                               nrow(plotDLEData)),
                                           rep("95% Credible Interval",
                                               nrow(plotDLEData) * 2)),
                                         levels=
                                           c("Estimate",
                                             "95% Credible Interval"))))
              
              ret1 <- ggplot2::qplot(x=x,
                                     y=y,
                                     data=gdata,
                                     group=group,
                                     linetype=Type,
                                     colour=I("red"),
                                     geom="line",
                                     xlab="Dose Levels",
                                     ylab="Probability of DLE [%]",
                                     ylim=c(0, 100))
              
              ret1 <- ret1 + ggplot2::scale_linetype_manual(breaks=c("Estimate",
                                                                     "95% Credible Interval"),
                                                            values=c(1,2), guide=ifelse(showLegend,
                                                                                        "legend", FALSE))
              ##only look at these dose levels for the plot:
              
              xLevels<- if (extrapolate) {
                seq_along(data@doseGrid)} else {1:max(data@xLevel)}
              
              ##get the plot data for the efficacy
              functionSamples <- matrix(nrow=sampleSize(Effsamples@options),
                                        ncol=length(xLevels))
              ##evaluate the efficacy for all samples
              for (i in seq_along(xLevels))
              { 
                ##Now we want to evaluate for the following dose
                functionSamples[,i] <- ExpEff(dose=data@doseGrid[xLevels[i]],
                                              model=Effmodel,
                                              samples=Effsamples)
              }
              ##extract mean curve
              meanCurve <- colMeans(functionSamples)
              
              ##extract quantiles
              quantiles=c(0.025, 0.975)
              quantCurve <- apply(functionSamples, 2L, quantile, prob=quantiles)
              
              ##now create the data frame
              plotEffData <-data.frame(dose=data@doseGrid[xLevels],
                                       mean=meanCurve,
                                       lower=quantCurve[1,],
                                       upper=quantCurve[2,])
              ##make the second plot
              ggdata<-with(plotEffData, data.frame(x=rep(dose,3),
                                                   y=c(mean,lower,upper),
                                                   group=
                                                     rep(c("mean","lower","upper"),
                                                         each=nrow(plotEffData)),
                                                   Type=
                                                     factor(c(rep("Estimate",
                                                                  nrow(plotEffData)),
                                                              rep("95% Credible Interval",
                                                                  nrow(plotEffData) * 2)),
                                                            levels=
                                                              c("Estimate",
                                                                "95% Credible Interval"))))
              
              plot2 <- ggplot2::qplot(x=x,
                                       y=y,
                                       data=ggdata,
                                       group=group,
                                       linetype=Type,
                                       colour=I("blue"),
                                       geom="line",
                                       xlab="Dose level",
                                       ylab="Expected Efficacy")
              
              plot2 <- plot2 +
                ggplot2::scale_linetype_manual(breaks=
                                                  c("Estimate",
                                                    "95% Credible Interval"),
                                                values=c(1,2),
                                                guide=ifelse(showLegend,
                                                             "legend", FALSE))
              
              ## arrange both plots side by side
              ret <- gridExtra::arrangeGrob(ret1, plot2, ncol=2)
              return(ret)
            })

##------------------------------------------------------------------------------
## Plot of the DLE and efficacy curve sides by side without  samples
## -----------------------------------------------------------------------------
##' Plot of the dose-DLE and dose-efficacy curve side by side given a DLE pseudo model 
##' and a given pseudo efficacy model without DLE and efficacy samples
##' 
##' @describeIn plotDualResponses Plot the DLE and efficacy curve side by side given a DLE model
##' and an efficacy model without any samples
##' 
##' @example examples/Samples-method-plotDualResponsesNoSamples.R
##'  
##' @export
##' @keywords methods 
setMethod("plotDualResponses",
          signature=
            signature(DLEmodel="ModelTox",
                      DLEsamples="missing",
                      Effmodel="ModelEff",
                      Effsamples="missing"),
          def=
            function(DLEmodel,Effmodel,data,...){
              
              ## Get Toxicity plot
              ## get the fit
              
              
              ##Make sure the model estimates are corresponds to the input data 
              DLEmodel <- update(object=DLEmodel,data=data)
              Effmodel <- update(object=Effmodel,data=data)
              
              
              plotDLEData <- data.frame(dose=data@doseGrid,
                                        probDLE=prob(dose=data@doseGrid,
                                                     model=DLEmodel))
              ## make the plot
              gdata <- with(plotDLEData,
                            data.frame(x=dose,
                                       y=probDLE,
                                       group=rep("Estimated DLE",each=nrow(plotDLEData)),
                                       Type=factor(rep("Estimated DLE",nrow(plotDLEData)),levels="Estimated DLE")))
              
              plot1 <- ggplot(data=gdata, aes(x=x,y=y), group=group) +
                xlab("Dose Levels")+
                ylab(paste("Probability of DLE")) + ylim(c(0,1)) + xlim(c(0,max(data@doseGrid))) +
                geom_line(colour=I("red"), size=1.5)
              
              
              plot1 <- plot1 +
                geom_line(size=1.5,colour="red")
              
              ##only look at these dose levels for the plot:
              
              ##get the plot data for the efficacy
              plotEffData<- data.frame(dose=data@doseGrid,
                                       ExpEff=ExpEff(dose=data@doseGrid,
                                                     model=Effmodel))
              
              ##make the second plot
              ggdata<-with(plotEffData,
                           data.frame(x=dose,
                                      y=ExpEff,
                                      group=rep("Estimated Expected Efficacy",each=nrow(plotEffData)),
                                      Type=factor(rep("Estimated Expected Efficacy",nrow(plotEffData)),levels="Estimated Expected Efficacy")))
              
              ##Get efficacy plot
              plot2 <- ggplot(data=ggdata, aes(x=x,y=y), group=group) +
                xlab("Dose Levels")+
                ylab(paste("Estimatimated Expected Efficacy")) + xlim(c(0,max(data@doseGrid))) +
                geom_line(colour=I("blue"), size=1.5)
              
              plot2 <- plot2 +
                geom_line(size=1.5,colour="blue")
              
              ## arrange both plots side by side
              ret <- gridExtra::arrangeGrob(plot1, plot2, ncol=2)
              return(ret)
            })
## =======================================================================================================









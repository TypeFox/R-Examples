###########################################
## ALL FUNCTIONS THAT HAVE BEEN REPLACED ##
###########################################

##' Deprecated functions, use replacements!
##'
##' Functions have been rename/replaced as:
##' \describe{
##'   \item{block.mult}{\code{\link{blockMult}}}
##' 
##'   \item{calc.F.times.X}{\code{\link{calc.FX}}}
##'   \item{calc.smooth.trends}{\code{\link{calcSmoothTrends}} and
##'                             \code{\link{updateTrend}}}
##'   \item{calc.tF.mat.F}{\code{\link{calc.tFXF}}}
##'   \item{calc.tF.times.mat}{\code{\link{calc.tFX}}, see also
##'     \code{\link{expandF}}.}
##'   \item{combineMesaData}{\code{\link{c.STmodel}}}
##'   \item{compute.ltaCV}{\code{\link{computeLTA}}}
##'   \item{cond.expectation}{\code{\link{predict.STmodel}}}
##'   \item{construct.LUR.basis}{\code{\link{createLUR}}}
##'   \item{construct.ST.basis}{\code{\link{createST}}}
##'   \item{create.data.matrix}{\code{\link{createDataMatrix}}}
##'   \item{create.data.model}{\code{\link{createSTmodel}}, see also
##'     \code{\link{updateCovf}} and \code{\link{processLocation}}.}
##'   \item{CVbasics}{Included in \code{\link{estimateCV.STmodel}}.}
##'   \item{CVresiduals.qqnorm}{\code{\link{qqnorm.STdata}},
##'     \code{\link{qqnorm.STmodel}}, or \code{\link{qqnorm.predCVSTmodel}}}
##'   \item{CVresiduals.scatter}{\code{\link{scatterPlot.STdata}},
##'     \code{\link{scatterPlot.STmodel}}, or
##'     \code{\link{scatterPlot.predCVSTmodel}}}
##' 
##'   \item{default.LUR.list}{\code{\link{processLUR}}}
##'   \item{default.ST.list}{\code{\link{processST}}}
##'   \item{detrend.data}{\code{\link{detrendSTdata}}}
##'   \item{dot.prod}{\code{\link{dotProd}}}
##'   \item{drop.observations}{\code{\link{dropObservations}}}
##' 
##'   \item{fit.mesa.model}{\code{\link{estimate.STmodel}}}
##' 
##'   \item{gen.gradient}{\code{\link{genGradient}}}
##'   \item{gen.hessian}{\code{\link{genHessian}}}
##'   \item{get.params}{\code{\link{loglikeSTgetPars}}}
##' 
##'   \item{loglike}{\code{\link{loglikeST}}}
##'   \item{loglike.dim}{\code{\link{loglikeSTdim}}}
##'   \item{loglike.grad}{\code{\link{loglikeSTGrad}}}
##'   \item{loglike.hessian}{\code{\link{loglikeSTHessian}}}
##'   \item{loglike.naive}{\code{\link{loglikeSTnaive}}}
##'   \item{loglike.naive.grad}{\code{\link{loglikeSTnaiveGrad}}}
##'   \item{loglike.naive.hessian}{\code{\link{loglikeSTnaiveHessian}}}
##'   \item{loglike.var.names}{\code{\link{loglikeSTnames}}}
##' 
##'   \item{make.sigma.B}{\code{\link{makeSigmaB}}, see
##'      \code{\link{parsCovFuns}}, \code{\link{namesCovFuns}}, and
##'      \code{\link{evalCovFuns}} for new covariance specifications.}
##'   \item{make.sigma.B.full}{\code{\link{makeSigmaB}} and
##'                            \code{\link{calc.FXtF2}}}
##'   \item{make.sigma.nu}{\code{\link{makeSigmaNu}}}
##'   \item{make.sigma.nu.cross.cov}{\code{\link{makeSigmaNu}}}
##' 
##'   \item{printMesaDataNbrObs}{\code{\link{print.STdata}},
##'                              \code{\link{summary.STdata}},
##'                              \code{\link{print.STmodel}}, and
##'                              \code{\link{summary.STmodel}}}
##'   \item{plotCV}{\code{\link{plot.predCVSTmodel}}}
##'   \item{plotMesaData}{\code{\link{plot.STdata}} and
##'                       \code{\link{plot.STmodel}}}
##'   \item{plotMonitoringLoc}{\code{\link{plot.STdata}} and
##'                            \code{\link{plot.STmodel}}}
##'   \item{plotPrediction}{\code{\link{plot.predictSTmodel}}}
##' 
##'   \item{remove.ST.mean}{\code{\link{removeSTcovarMean}}}
##'   \item{run.MCMC}{\code{\link{MCMC.STmodel}}}
##' 
##'   \item{setupSTdataset}{\code{\link{createSTdata}}}
##'   \item{simulateMesaData}{\code{\link{simulate.STmodel}}}
##'   \item{summaryStatsCV}{\code{\link{summary.predCVSTmodel}}}
##'   \item{SVD.miss}{\code{\link{SVDmiss}}}
##'   \item{SVD.smooth}{\code{\link{SVDsmooth}}}
##'   \item{SVD.smooth.cv}{\code{\link{SVDsmoothCV}}, see also
##'     \code{\link{plot.SVDcv}} and \code{\link{print.SVDcv}}.}
##' 
##'   \item{tstat}{Included in \code{\link{predict.STmodel}}}
##' }
##'
##' @title Deprecated functions, use replacements!
##' 
##' @param ... Unused, for compability.
##' 
##' @return Does not return.
##' 
##' @author Johan Lindström
##'
##' @export
make.sigma.B <- function(...){
  stop("USE makeSigmaB")
}

##' @rdname make.sigma.B
##' @export
make.sigma.B.full <- function(...){
  stop("USE makeSigmaB + calc.FXtF2")
}

##' @rdname make.sigma.B
##' @export
make.sigma.nu <- function(...){
  stop("USE makeSigmaNu")
}

##' @rdname make.sigma.B
##' @export
make.sigma.nu.cross.cov <- function(...){
  stop("USE makeSigmaNu")
}

##' @rdname make.sigma.B
##' @export
calc.tF.times.mat <- function(...){
 stop("USE calc.tFX")
}

##' @rdname make.sigma.B
##' @export
calc.F.times.X <- function(...){
    stop("USE calc.FX")
}

##' @rdname make.sigma.B
##' @export
calc.tF.mat.F <- function(...){
  stop("USE calc.tFXF")
}

##' @rdname make.sigma.B
##' @export
block.mult <- function(...){
  stop("USE blockMult")
}

##' @rdname make.sigma.B
##' @export
dot.prod <- function(...){
  stop("USE dotProd")
}

##' @rdname make.sigma.B
##' @export
SVD.miss <- function(...){
  stop("USE SVDmiss")
}

##' @rdname make.sigma.B
##' @export
SVD.smooth <- function(...){
  stop("USE SVDsmooth")
}

##' @rdname make.sigma.B
##' @export
SVD.smooth.cv <- function(...){
  stop("USE SVDsmoothCV")
}

##' @rdname make.sigma.B
##' @export
calc.smooth.trends <- function(...){
  stop("USE calcSmoothTrends and/or updateTrend")
}

##' @rdname make.sigma.B
##' @export
setupSTdataset <- function(...){
  stop("USE createSTdata")
}

##' @rdname make.sigma.B
##' @export
printMesaDataNbrObs <- function(...){
  stop("USE print.STdata+summary.STdata+print.STmodel+summary.STmodel")
}

##' @rdname make.sigma.B
##' @export
plotMonitoringLoc <- function(...){
  stop("USE plot.STdata+plot.STmodel")
}
  
##' @rdname make.sigma.B
##' @export
plotMesaData <- function(...){
  stop("USE plot.STdata+plot.STmodel")
}

##' @rdname make.sigma.B
##' @export
create.data.matrix <- function(...){
  stop("USE createDataMatrix")
}

##' @rdname make.sigma.B
##' @export
remove.ST.mean <- function(...){
  stop("USE removeSTcovarMean")
}

##' @rdname make.sigma.B
##' @export
detrend.data <- function(...){
  stop("USE detrendSTdata")
}

##' @rdname make.sigma.B
##' @export
create.data.model <- function(...){
  stop("USE createSTmodel")
}

##' @rdname make.sigma.B
##' @export
default.LUR.list <- function(...){
  stop("USE processLUR")
}

##' @rdname make.sigma.B
##' @export
default.ST.list <- function(...){
  stop("USE processST")
}

##' @rdname make.sigma.B
##' @export
construct.LUR.basis <- function(...){
  stop("USE createLUR")
}

##' @rdname make.sigma.B
##' @export
construct.ST.basis <- function(...){
  stop("USE createST")
}

##' @rdname make.sigma.B
##' @export
loglike.dim <- function(...){
  stop("USE loglikeSTdim")
}

##' @rdname make.sigma.B
##' @export
loglike.var.names <- function(...){
  stop("USE loglikeSTnames")
}

##' @rdname make.sigma.B
##' @export
get.params <- function(...){
  stop("USE loglikeSTgetPars")
}

##' @rdname make.sigma.B
##' @export
loglike <- function(...){
  stop("USE loglikeST")
}

##' @rdname make.sigma.B
##' @export
loglike.naive <- function(...){
  stop("USE loglikeSTnaive")
}

##' @rdname make.sigma.B
##' @export
gen.gradient <- function(...){
  stop("USE genGradient")
}

##' @rdname make.sigma.B
##' @export
gen.hessian <- function(...){
  stop("USE genHessian")
}

##' @rdname make.sigma.B
##' @export
loglike.grad <- function(...){
  stop("USE loglikeSTGrad")
}

##' @rdname make.sigma.B
##' @export
loglike.hessian <- function(...){
  stop("USE loglikeSTHessian")
}

##' @rdname make.sigma.B
##' @export
loglike.naive.grad <- function(...){
  stop("USE loglikeSTnaiveGrad")
}

##' @rdname make.sigma.B
##' @export
loglike.naive.hessian <- function(...){
  stop("USE loglikeSTnaiveHessian")
}

##' @rdname make.sigma.B
##' @export
fit.mesa.model <- function(...){
  stop("USE estimate.STmodel")
}

##' @rdname make.sigma.B
##' @export
cond.expectation <- function(...){
  stop("USE predict.STmodel")
}

##' @rdname make.sigma.B
##' @export
simulateMesaData <- function(...){
  stop("USE simulate.STmodel")
}

##' @rdname make.sigma.B
##' @export
combineMesaData <- function(...){
  stop("USE c.STmodel")
}

##' @rdname make.sigma.B
##' @export
drop.observations <- function(...){
  stop("USE dropObservations")
}

##' @rdname make.sigma.B
##' @export
plotPrediction <- function(...){
  stop("USE plot.predictSTmodel")
}

##' @rdname make.sigma.B
##' @export
tstat <- function(...){
  stop("INCLUDED in STmodel.predict")
}

##' @rdname make.sigma.B
##' @export
compute.ltaCV <- function(...){
  stop("USE computeLTA")
}

##' @rdname make.sigma.B
##' @export
CVbasics <- function(...){
  stop("INCLUDED in estimateCV.STmodel")
}

##' @rdname make.sigma.B
##' @export
summaryStatsCV <- function(...){
  stop("USE summary.predCVSTmodel")
}

##' @rdname make.sigma.B
##' @export
run.MCMC <- function(...){
  stop("USE MCMC.STmodel")
}

##' @rdname make.sigma.B
##' @export
plotCV <- function(...){
  stop("USE plot.predCVSTmodel")
}

##' @rdname make.sigma.B
##' @export
CVresiduals.qqnorm <- function(...){
  stop("USE qqnorm.predCVSTmodel")
}

##' @rdname make.sigma.B
##' @export
CVresiduals.scatter <- function(...){
  stop("USE scatterPlot.predCVSTmodel")
}

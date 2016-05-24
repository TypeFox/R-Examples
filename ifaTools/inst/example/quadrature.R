library(ifaTools)
library(ggplot2)
library(plyr)

if (0) {
  control <- expand.grid(
    # conditions
    qpoints = seq(3,27,1), qwidth = seq(3,5,1), seed=1:5, modelFactors = 1:2,
    # outcomes
    fit= NA, code = NA, grad=NA, posdef = NA, condnum=NA, pnorm=NA
  )
}

# generate random data for a 2 factor latent structure
numItems <- 12

mkdata <- function(seed, dataFactors) {
  set.seed(seed)
  
  # rpf.grm creates a graded response model
  spec <- list()  # a map from columns to item models
  spec[1:numItems] <- rpf.grm(dataFactors, outcomes=5)

  trueIMat <- matrix(1.5, nrow=4+dataFactors, ncol=numItems)
  trueIMat[1+dataFactors,] <- seq(-2,2, length.out = numItems/2)
  for (rx in (1+dataFactors):(4+dataFactors)) {
    trueIMat[rx,] <- trueIMat[rx-1,] - .3
  }
  rownames(trueIMat) <- c(paste0('f',1:dataFactors), paste0("b",1:4))
  colnames(trueIMat) <- paste0('i', 1:numItems)
  names(spec) <- colnames(trueIMat)
  
  perFactor <- 500
  
  data <- rpf.sample(perFactor * dataFactors, spec, trueIMat)
  
  data <- compressDataFrame(data)
}

mkmodel <- function(seed, data, modelFactors) {
  set.seed(seed)

  spec <- list()  # a map from columns to item models
  spec[1:numItems] <- rpf.grm(1, outcomes=5)
  
  computePlan <- mxComputeSequence(list(
    mxComputeEM('model.expectation', 'scores', mxComputeNewtonRaphson(),
                verbose=0, information='oakes1999',
                infoArgs=list(fitfunction='fitfunction')),
    mxComputeHessianQuality(),
    mxComputeOnce('fitfunction', 'gradient'),
    mxComputeReportDeriv()))
  
  imat <- mxMatrix(name="item", free=TRUE,
                   values=mxSimplify2Array(lapply(spec, rpf.rparam)))
  imat$labels[,] <- paste0('p',1:prod(dim(imat)))
  colnames(imat) <- colnames(data)[-which(colnames(data) == "freq")]
  
  template <- mxModel(model="model", imat,
                      mxMatrix(values=seed, nrow=1, ncol=1, name="seed"),
                      mxData(observed=data, type="raw", numObs = sum(data[['freq']]),
                             sort = FALSE),
                      mxExpectationBA81(ItemSpec=spec, weightColumn = "freq"),
                      mxFitFunctionML(),
                      computePlan)

  addExploratoryFactors(template, modelFactors-1)
}

setQuadrature <- function(model, pts, width) {
  model$model$expectation$qpoints <- as.integer(pts)
  model$model$expectation$qwidth <- width
  model
}

refModel <- NULL

for (cx in head(which(is.na(control$fit)),1):nrow(control)) {
  cntl <- as.list(control[cx,])
  if (!is.na(cntl$fit)) next
  
  if (is.null(refModel) ||
      refModel$model$expectation$ItemSpec[[1]]$factors != cntl$modelFactors) {
    data <- mkdata(1, cntl$modelFactors)
    refModel <- mkmodel(1, data, cntl$modelFactors)
    refModel <- setQuadrature(refModel, 41, 5)
    refModel <- mxRun(refModel)
  }
  
  model <- mkmodel(cntl$seed, data, cntl$modelFactors)
  model <- setQuadrature(model, cntl$qpoints, cntl$qwidth)
  
  fit <- mxRun(model, suppressWarnings = TRUE, silent = TRUE)
  control[cx,'fit'] <- fit$output$fit
  control[cx,'code'] <- fit$output$status$code
  control[cx,'posdef'] <- fit$output$infoDefinite
  control[cx,'condnum'] <- fit$output$conditionNumber
  control[cx,'grad'] <- norm(fit$output$gradient, "2")
  
  if (cntl$modelFactors > 1) {
    control[cx,'pnorm'] <-
      suppressWarnings(
        norm(varimax(toFactorLoading(fit$model$item$values[1:cntl$modelFactors,]))$loadings[] -
               varimax(toFactorLoading(refModel$model$item$values[1:cntl$modelFactors,]))$loadings[], "2"))
  } else {
    control[cx,'pnorm'] <-
      suppressWarnings(
        norm(toFactorLoading(fit$model$item$values[1:cntl$modelFactors,]) -
               toFactorLoading(refModel$model$item$values[1:cntl$modelFactors,]), "2"))
  }
  
#  print(control[cx,])
  
  if (1) {
    gcontrol <- subset(control, !is.na(fit))
    gcontrol <- ddply(gcontrol, ~qpoints + qwidth + modelFactors, summarize, pnorm=mean(pnorm))
    gcontrol <- ddply(gcontrol, ~modelFactors, transform, pnorm=log(pnorm))
    #  gcontrol$fit <- scale(control$fit - median(control$fit))
    if (nrow(gcontrol) > 3) {
      pl <- ggplot(gcontrol,
                   aes(x=qpoints, y=qwidth, fill=pnorm)) + geom_tile() +
        facet_wrap(~modelFactors)
      print(pl)
    }
  }  
}

if (0) {
  save(control, file="~/ifa/mbManual/quad.rda")
}

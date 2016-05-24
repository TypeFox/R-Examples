library(OpenMx)
library(ifaTools)
#library(GPArotation)
library(mvtnorm)
library(grid)

if (0) {
  result <- NULL
}

# --------------------------------------------------------------- generate data

mkTemplate <- function(seed) {
  set.seed(seed)
  
  # suppose our items are related to language development
  factorNames <- c('spelling', 'phonics', 'reading')
  
  # generate random data for a 2 factor latent structure
  numItems <- 12
  spec <- list()  # a map from columns to item models
  
  # rpf.grm creates a graded response model
  dataFactors <- ifelse(seed %% 2 == 0, 2, 1)
  spec[1:numItems] <- rpf.grm(dataFactors, outcomes=8)
  
  trueIMat <- matrix(1.5, nrow=rpf.numParam(spec[[1]]), ncol=numItems)
  trueIMat[1+dataFactors,] <- seq(0, 0.3*7, length.out = numItems/2)
  for (rx in (2+dataFactors):nrow(trueIMat)) {
    trueIMat[rx,] <- trueIMat[rx-1,] - 0.3
  }
  
  rownames(trueIMat) <- c(factorNames[1:dataFactors],
                          paste0("b",(dataFactors+1):nrow(trueIMat)))
  colnames(trueIMat) <- paste0('i', 1:numItems)
  names(spec) <- colnames(trueIMat)
  
  # half the items mostly assess spelling and the other half of the items assess phonics
  if (dataFactors > 1) {
    trueDiscrimination0 <- apply(trueIMat[1:dataFactors,], 2, norm, "2")
    trueIMat['spelling',1:(numItems/2)] <- 0
    trueIMat['phonics',(1+numItems/2):numItems] <- 0
    
    # we need to rescale the intercept to keep it in the original place
    trueIMat[(dataFactors+1):nrow(trueIMat),] <-
      apply(trueIMat[1:dataFactors,], 2, norm, "2") * trueIMat[(dataFactors+1):nrow(trueIMat),] / trueDiscrimination0
  }
  
  draws <- 1500

    # simulate some data
  if (0) {
    basis1 <- runif(2)
    basis1 <- basis1 / norm(basis1, "2")
    basis2 <- c()
    for (retry in 1:10) {
      basis2 <- runif(2)
      basis2 <- basis2 / norm(basis2, "2")
      if (basis1 %*% basis2 < .8) break
    }
  } else {
    if (dataFactors == 1) {
      theta <- rnorm(draws)
    } else {
      theta <- rmvnorm(draws, sigma=diag(dataFactors))
    }
  }
  
  data <- rpf.sample(t(theta), spec, trueIMat)

  if (0) {
    require(ggplot2)
    # Note severe signal degredation even using the true item parameters!
    latent <- data.frame(theta, factor=thetaFactor)
    ggplot(latent, aes(x=X1, y=X2, color=factor)) +
      geom_point(size=3, alpha=.15) + coord_fixed()
    tgrp <- list(spec=spec, param=trueIMat, data=data)
    sc <- as.data.frame(EAPscores(tgrp)[,1:2])
    ggplot(sc, aes(x=spelling, y=phonics)) +
      geom_point(size=3, alpha=.15) + coord_fixed()
  }
  
  # --------------------------------------------------------------- 1 factor
  
  # Needed to compute the number of statistics (unique rows)
  data <- compressDataFrame(data)
  
  computePlan <- mxComputeSequence(list(
    mxComputeEM('model.expectation', 'scores', mxComputeNewtonRaphson(),
                verbose=0, information='oakes1999', #tolerance=1e-8,
                infoArgs=list(fitfunction='fitfunction')),
    mxComputeHessianQuality(),
    mxComputeOnce('fitfunction','gradient'),
    mxComputeReportDeriv()))
  
  spec[1:numItems] <- lapply(spec, rpf.modify, factors=1)
  
  imat <- mxMatrix(name="item", free=TRUE,
                   values=mxSimplify2Array(lapply(spec, rpf.rparam)))
  imat$labels[,] <- paste0('p',1:prod(dim(imat)))
  
  template <- mxModel(model="model", imat,
                      mxData(observed=data, type="raw", numObs = draws),
                      mxExpectationBA81(ItemSpec=spec, weightColumn = "freq"),
                      mxFitFunctionML(),
                      computePlan)
  template
}

# --------------------------------------------------------------- power simulation

if (0) {
  start <- 1
  if (!is.null(result)) {
    start <- 1+nrow(result)
  }
  
  for (repl in start:100) {
    template <- mkTemplate(repl)
    
    m1 <- addExploratoryFactors(template, 0)
    m1 <- mxRun(m1, suppressWarnings = TRUE)
    m2 <- addExploratoryFactors(template, 1)
    m2 <- mxRun(m2, suppressWarnings = TRUE)
    #   m3 <- addExploratoryFactors(template, 2)
    #   m3 <- rampRun(m3)
    allm <- list(m1,m2)
    
    condnum <- sapply(allm, function (m) {
      if (is.null(m$output$conditionNumber)) { NA }
      else { log(m$output$conditionNumber) } })
    grad <- sapply(allm, function (m) {
      if (is.null(m$output$gradient)) { NA }
      else { norm(m$output$gradient, "2") }
    })
    
    r1 <- data.frame(seed=repl,
                     status=max(sapply(allm, function (m) m$output$status$code)),
                     grad1=grad[1],
                     grad2=grad[2],
                     cn1=condnum[1],
                     cn2=condnum[2],
                     '2v1'=mxCompare(m2,c(m1))[2,'p'])
    print(r1)
    result <- rbind(result, r1)
  }
}

# --------------------------------------------------------------- factor rotation

plotTwoFactors <- function(slope) {
  lvm <- varimax(toFactorLoading(slope))$loadings
  if (any(abs(lvm[lvm<0]) > .001)) stop("Got negative loadings")
  lvm[lvm<0] <- 0
  df <- as.data.frame(lvm[,1:2])
  df$name <- rownames(df)
  pl <- ggplot(df, aes_string(x=rownames(slope)[1],
                              y=rownames(slope)[2], label="name")) + geom_text()
  pm <- promax(lvm[,1:2])$rotmat
  for (dx in 1:ncol(pm)) {
    d1 <- .5 * pm[,dx] / sqrt(sum(pm[,dx]^2))
    pl <- pl + geom_segment(x=.5, y=.5, xend=d1[1] + .5, yend=d1[2] + .5,
                            arrow = arrow(length = unit(.5,"cm")))
  }
  pl + xlim(0,1) + ylim(0,1)
}

template <- mkTemplate(1)  # single factor

m1 <- addExploratoryFactors(template, 0)
m1 <- mxRun(m1, suppressWarnings = TRUE)
m2 <- addExploratoryFactors(template, 1)
m2 <- mxRun(m2, suppressWarnings = TRUE)
mxCompare(m2,m1)
plotTwoFactors(m2$model$item$values[1:2,])

template <- mkTemplate(2)  # 2 factors

m1 <- addExploratoryFactors(template, 0)
m1 <- mxRun(m1, suppressWarnings = TRUE)
m2 <- addExploratoryFactors(template, 1)
m2 <- mxRun(m2, suppressWarnings = TRUE)
mxCompare(m2,m1)
plotTwoFactors(m2$modelItem$item$values[1:2,])

if (0) {
  m0 <- mxRefModels(mxRename(template, "model0"), run=TRUE)
  
  mxCompare(m1,c(m0$Independence))  # likelihood ratio test
  m0G <- as.IFAgroup(m0$Independence)
  m1G <- as.IFAgroup(m1$modelItem)
  m2G <- as.IFAgroup(m2$modelItem)
  m3G <- as.IFAgroup(m3$modelItem)
  
  SitemFit(m1G)
  SitemFit(m2G)
  
  # We don't trust these tests, but here's how to do it
  multinomialFit(m1G, m0G)
  multinomialFit(m2G, m0G)
  multinomialFit(m3G, m0G)
  
  ChenThissen1997(m1G)
  ChenThissen1997(m2G)
}

if (0) {
  require("mirt")
  mdata <- expandDataFrame(as.data.frame(lapply(data, unclass)), freqName = "freq")
  mod <- mirt(mdata, 2, TOL=c(.001, .01))
  summary(mod, rotate="none")
  
  clist <- coef(mod)
  clist$GroupPars <- NULL
  FL <- toFactorLoading(mxSimplify2Array(clist)[1:2,])
}

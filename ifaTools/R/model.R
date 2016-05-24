placeInContainer <- function(model) {
    if (is(model$fitfunction, "MxFitFunctionMultigroup")) {
	    return(model)
    }

    compute <- mxComputeDefault()
    if (!is.null(model$compute)) compute <- model$compute

    mxModel(model=paste0(model$name, "Container"), model,
	    mxFitFunctionMultigroup(model$name), compute)
}

extractFromContainer <- function(container, expectationClass) {
	if (!is(container$fitfunction, "MxFitFunctionMultigroup")) {
		stop(paste(container$name, "has", container$fitfunction,
			   "instead of MxFitFunctionMultigroup"))
	}
	match <- list()
	for (sm in names(container$submodels)) {
		if (is(container[[sm]]$expectation, expectationClass)) {
			match <- c(match, container[[sm]])
		}
	}
	if (length(match) > 1) {
		stop(paste("More than one model within", container$name,
			   "has an", expectationClass))
	}
	match[[1]]
}

#' Adds exploratory factors to a single factor model
#' 
#' @param model a single factor (possibly multigroup) model
#' @param toAdd the number of factors to add
#' @param ...  Not used.  Forces remaining arguments to be specified by name.
#' @param addUniquenessPrior whether to add a uniqueness prior to the model (default TRUE)
#'
#' @importFrom methods is
#' @export
addExploratoryFactors <- function(model, toAdd, ..., addUniquenessPrior=TRUE) {
    garbageArguments <- list(...)
    if (length(garbageArguments) > 0) {
        stop("Values for the '...' argument are invalid; use named arguments")
    }

    model <- placeInContainer(model)
    imodel <- extractFromContainer(model, "MxExpectationBA81")

    ex <- imodel$expectation

  if (any(sapply(ex$ItemSpec, function(s) s$factors) != 1)) {
    stop(paste("Model", imodel$name, "has items that load",
               "on more or less than 1 factor"))
  }

  item <- imodel[[ex$item]]
  if (is.null(item)) {
    stop(paste("Cannot find item matrix in model", imodel$name))
  }

    unlabeled <- sum(is.na(item$labels[item$free]))
    if (unlabeled) {
	    item$labels[item$free & is.na(item$labels)] <- paste0("p", 1:unlabeled)
	    imodel <- mxModel(imodel, item)
    }

  if (toAdd > 0) {
    ex$ItemSpec <- lapply(ex$ItemSpec, rpf.modify, factors=1+toAdd)
    nitem <- mxMatrix(nrow=nrow(item)+toAdd, ncol=ncol(item),
                      name=item$name, dimnames=list(NULL, colnames(item)))
    for (layer in c("values", "labels", "free", "lbound", "ubound")) { 
      nitem[[layer]][1,] <- item[[layer]][1,]
      nitem[[layer]][(2+toAdd):nrow(nitem),] <- item[[layer]][2:nrow(item),]
    }
    rownames(nitem) <- c(rownames(item)[1],
                         paste0("explore", 1:toAdd),
                         rownames(item)[2:nrow(item)])
    
    for (fx in 1:toAdd) {
      lab <- paste0("E",fx)
      mask <- item$free[1,]
      mask[which(mask)[1:fx]] <- FALSE
      nitem$free[1+fx,] <- mask
      nitem$values[1+fx, mask] <- nitem$values[1, mask]
      nitem$labels[1+fx, mask] <- paste0(lab, "p", 1:sum(mask))
    }
    
    # Probably don't bother with more than 3 factors. Its too
    # slow and inaccurate.
    ex$qpoints <- min(ex$qpoints,
                      switch(as.character(1+toAdd),
                             '2'=31, '3'=19, '4'=15, 9))
    ex$qwidth <- min(ex$qwidth,
                     switch(as.character(1+toAdd),
                            '2'=5, '3'=4, '4'=3, 3))
    
    imodel <- mxModel(imodel, ex, nitem)
    item <- nitem
    
    mMat <- imodel[[ex$mean]]
    if (!is.null(mMat)) {
      nmMat <- mxMatrix(nrow=1, ncol=1+toAdd, values=0, free=TRUE,
                        name=mMat$name)
      colnames(nmMat) <- c(ifelse(length(rownames(mMat)),
                                  rownames(mMat), colnames(mMat)),
                           paste0("explore", 1:toAdd))
      # preserve labels? TODO
      imodel <- mxModel(imodel, nmMat)
    }
    
    cMat <- imodel[[ex$cov]]
    if (!is.null(cMat)) {
      ncMat <- mxMatrix(type="Symm", nrow=1+toAdd, ncol=1+toAdd,
                        values=diag(1+toAdd), free=TRUE,
                        dimnames=list(colnames(nmMat), colnames(nmMat)),
                        name=cMat$name)
      for (rx in 1:toAdd) {
        for (cx in (rx+1):(1+toAdd)) {
          lab <- paste(imodel$name, cMat$name,
                       paste0("r", rx, "c", cx), sep="_")
          ncMat$labels[rx,cx] <-
            ncMat$labels[cx,rx] <- lab
        }
      }
      imodel <- mxModel(imodel, ncMat)
    }
  }

    model <- mxModel(model, imodel)

    if (addUniquenessPrior) {
	    name <- imodel$name
	    up <- uniquenessPrior(imodel, 1+toAdd, name=paste(name,"UniquenessPrior",sep=""))
	    model <- mxModel(model, up)
	    model$fitfunction$groups <- c(model$fitfunction$groups, up$name)
    }

    mxRename(model, paste0(model$name,1+toAdd))
}

#' Replicate a model for each group of data
#'
#' The reference group is fixed to a zero mean and identity covariance matrix.
#' 
#' @param tmpl an OpenMx model
#' @param fullData the complete data including the column indicating group membership
#' @param mMat an MxMatrix for latent means
#' @param covMat an MxMatrix for latent covariance
#' @param ...  Not used.  Forces remaining arguments to be specified by name.
#' @param splitCol the name of the column used to indicate group membership
#' @param refGroup the name of the reference group
#' @param split whether to split the data (defaults to TRUE)
#' @param compressData whether to apply compressDataFrame (defaults to TRUE)
#' 
#' @export
replicateModelBy <- function(tmpl, fullData, mMat, covMat, ...,
			     splitCol="population", refGroup="general", split=TRUE,
			     compressData=TRUE) {
    garbageArguments <- list(...)
    if (length(garbageArguments) > 0) {
        stop("Values for the '...' argument are invalid; use named arguments")
    }

  dataCol <- setdiff(colnames(fullData), splitCol)
  itemMatName <- tmpl$expectation$item
    if (compressData) {
	    tmpl$expectation$weightColumn <- 'freq'
    }
  ivalues <- tmpl[[itemMatName]]$values
  newLabels <- tmpl[[itemMatName]]$labels
  
  for (r in 1:nrow(newLabels)) {
    for (c in 1:ncol(newLabels)) {
      if (!is.na(newLabels[r,c]) || is.na(ivalues[r,c])) next
      newLabels[r,c] <- paste(itemMatName,"_r", r, "c", c, sep = "")
    }
  }
  
  tmpl[[itemMatName]]$labels <- newLabels
  
  mask <- fullData[[splitCol]] == refGroup
  refData <- subset(fullData, !split | mask, dataCol)
    if (compressData) {
	    refData <- compressDataFrame(refData)
    }
  refModel <- mxModel(model=tmpl, name=refGroup,
                      mxData(observed=refData, type="raw", sort=!compressData))
    if (compressData) {
	    refModel$data$numObs <- sum(refData[['freq']])
    }
    
  extraGroups <- setdiff(unique(fullData[[splitCol]]), refGroup)
  container <- mxModel(model="container", refModel)
  if (!split) {
    extraGroups <- NULL
  } else {
    container <- mxModel(container,
                         mxModel(model="latentFit",
                                 mxFitFunctionMultigroup(paste0(extraGroups, "Latent"))))
  }

  for (grp in extraGroups) {
    mask <- fullData[[splitCol]] == grp
	  gdata <- subset(fullData, mask, dataCol)
	  if (compressData) {
		  gdata <- compressDataFrame(gdata)
	  }
	  gmodel <- mxModel(model=tmpl, name=grp, mMat, covMat,
			    mxData(observed=gdata, type="raw", sort=!compressData))
	  if (compressData) {
		  gmodel$data$numObs <- sum(gdata[['freq']])
	  }
    lmodel <- mxModel(model=paste0(grp,"Latent"),
                      mxDataDynamic(type='cov', expectation=paste(grp,'expectation',sep='.')),
                      mxExpectationNormal(covariance=paste(grp,'cov',sep='.'),
                                          means=paste(grp,'mean',sep='.')),
                      mxFitFunctionML())
    container <- mxModel(model=container, gmodel, lmodel)
  }
  
  enames <- c(refGroup, extraGroups)
  mstep <- mxComputeNewtonRaphson()
  container <- mxModel(container,
                       mxFitFunctionMultigroup(c(enames, 'latentFit')),
		       mxComputeSequence(list(
			   mxComputeEM(paste(enames, 'expectation', sep="."), 'scores',
				       mstep, information="oakes1999", infoArgs=list(fitfunction='fitfunction')),
			   mxComputeStandardError(),
			   mxComputeHessianQuality())))
  container
}

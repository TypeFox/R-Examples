"getModel" <- function (data, modelspec, modeldiffs, datasetind, opt) 
{
  if(length(datasetind) != length(data))
    datasetind <- rep(1, length(data))
  modellist <- vector("list", length(data))
  resultlist <- vector("list", length(data))
  plugin <- c("psi.df", "x", "x2", "nt", "nl", "simdata", "inten", 
              "datafile", "satMat", "outMat")
  for(i in 1:length(data)) {
    resultlist[[i]] <- res() 
    modellist[[i]] <- modelspec[[datasetind[i] ]]
    for (sl in plugin) 
      slot(modellist[[i]], sl) <- slot(data[[i]], sl) # TODO: avoid data duplication!
  }
  #if ( (currModel@modelspec)[[1]]@dispmu && length((currModel@modelspec)[[1]]@parmu[[1]])>0 && length((currModel@modelspec)[[1]]@lambdac)==0){
  #  try(
  #    (currModel@modelspec)[[1]]@lambdac <- slot((currModel@data)[[1]],(currModel@modelspec)[[1]]@clpType)[[1]]
  #  )
  #}
  modellist <- addDscal(modellist, modeldiffs$dscal)
  if (length(modeldiffs$free) != 0) 
    modellist <- diffFree(modellist, modeldiffs$free)
  if (length(modeldiffs$remove) != 0) 
    modellist <- diffRemove(modellist, modeldiffs$remove)
  if (length(modeldiffs$add) != 0) 
    modellist <- diffAdd(modellist, modeldiffs$add)
  if (length(modeldiffs$change) != 0) {
    reCh <- diffChange(modellist, modeldiffs$change)
    modellist <- reCh$modellist
    modeldiffs$change <- reCh$diffsadd
  }
  if (length(modeldiffs$rel) != 0) 
    modellist <- diffRel(modellist, modeldiffs$rel)
  if (length(modeldiffs$dscal) == 0) 
    modeldiffs$dscal <- list()
  if(length(modeldiffs$weightList) > 0)
    modellist <- addOutlierWeights(modellist, modeldiffs$weightList) 
  modellist <- lapply(modellist, applyWeightingModel)
  modellist <- lapply(modellist, initModellist)
  modellist <- getPrelBetweenDatasets(modellist, modeldiffs$rel) 
  linkclp <- if(length(modeldiffs$linkclp) < 1) 
    list(1:length(modellist))
  else modeldiffs$linkclp
  grlist <- list()
  
  for(i in 1:length(linkclp)) {
    mlist <- list()
    mlabel <- vector()
    for(j in linkclp[[i]]) {
      mlist <- append(mlist, modellist[[j]] )
      mlabel <- append(mlabel, j)
    }
    gr <- getGroups(mlist, modeldiffs, mlabel)
    grlist <- append(grlist, gr)
  }
  if(length(modeldiffs$getXsuper) > 1) 
    getXsuper <- modeldiffs$getXsuper
  else getXsuper <- FALSE
  if(opt@writedata)
    writeData(data,opt)
  multimodel(modellist = modellist, data = data, datasetind = datasetind, 
             modelspec = modelspec, stderrclp = opt@stderrclp,
             modeldiffs = modeldiffs,
             nnlscrit = opt@nnlscrit, nnls = opt@nnls, optlist = list(opt),
             groups = grlist, getXsuper = getXsuper,
             trilinear = opt@trilinear,
             fit = fit(resultlist=resultlist), algorithm = opt@algorithm)
}

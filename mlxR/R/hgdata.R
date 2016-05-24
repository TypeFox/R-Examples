hgdata <- function(lv)
{
  group <- lv$group
  N <- length(group)
  ip <- vector("list",N)
  id <- vector("list",N)
  for (i in seq(1,N)){
    if (!isfield(group[[i]],'level'))
      group[[i]]$level <- 'longitudinal'
    
    if (!isfield(group[[i]],'size'))
      group[[i]]$size <- 1
    
    id[[i]] <- list()
    id[[i]]$dose <- lv$treatment[[i]]
    id[[i]]$response <- lv$output[[i]]
    id[[i]]$name <- i
  }
  dataIn=list(group=group,individual=id)
  if (!is.null(lv$parameter)){
    pp <- lv$parameter
    pp$id <- NULL
    cc <- processing_categorical(pp)
    if (!is.null(cc))
      dataIn$catvar <- cc
    dm <- data.matrix(pp)
    ip <- list(name=names(pp), value=matrix(dm,nrow=N,ncol=ncol(dm)))
    dataIn$individual_parameters=ip
  }
  if (isfield(lv,"depot"))
    dataIn$depot <- lv$depot
  if (isfield(lv,"regressor"))
    dataIn$regressor <- lv$regressor
  if (isfield(lv,"varlevel")){
    lvar <- lv$varlevel
#     lvv <- lvar$value
#     dvar <- diff(lvv[,which(lvar$colTypes=='OCC')])
#     for (k in seq(1,ncol(dvar)-1)){
#       i0 <- which(dvar[,k]==1)
#     }
    dataIn$variability <- lvar
  }
  return(dataIn)
}

#------------------------------------------------------
processing_categorical <- function(param)
{
  pf=param[sapply(param, is.factor)]
  npf <- names(pf)
  if (length(npf)>0){
    catvar <- list(name=NULL,categories=list())
    catvar$name <- c(catvar$name , npf)
    for (j in 1:length(npf)){
      catvar$categories = append(catvar$categories, list(levels(pf[[j]])))
    }
  }else{
    catvar <- NULL
  }
  return(catvar)
}
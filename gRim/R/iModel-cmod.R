##########################################################
##
## Continuous interaction model (graphical Gaussian model)
##
##########################################################

cmod <- function(formula, data, marginal=NULL, fit=TRUE, details=0){

  if (inherits(data,"data.frame")){
    tmp   <- cov.wt(data, method="ML")
    S     <- tmp$cov
    n.obs <- tmp$n.obs
  } else {
    S     <- data$cov
    n.obs <- data$n.obs
  }

  varNames <- colnames(S)
  ans      <- .pFormula2(formula, varNames, marginal)
  ##, v.sep = ":",  g.sep = "+", ignore.power.value=TRUE) 
  glist <- ans$glist
  ## Get varNames in the order matching to the data:
  varNames <- varNames[sort(match(ans$varNames, varNames))]
  
  datainfo <- list(S=S[varNames,varNames],
                   n.obs=n.obs, data=data)
    
  res <- list(glist          = glist,
              varNames       = varNames,
              datainfo       = datainfo,
              fitinfo        = NULL,
              isFitted       = FALSE
              )

  upd   <- .cModel_finalize(glist, varNames)  
  res[names(upd)] <- upd  
  class(res) <- c("cModel","iModel")
  
  if (fit){
    res <- fit(res)
  }
  res
}


.cModel_finalize <- function(glist, varNames){

  ugm   <- ugList(glist, result="matrix")
  glist <- maxCliqueMAT(ugm)[[1]]
  isd   <- length(mcsMAT(ugm))>0

  glistNUM <- lapply(glist,
                     function(ll) {
                       match(ll, varNames)
                     })

  ret      <- list(glist          = glist,
                   glistNUM       = glistNUM,
                   isDecomposable = isd,
                   isGraphical    = TRUE
                   )
  return(ret)
}


fit.cModel <- function(object, engine="ggmfit",start=NULL, ...){

  switch(engine,
         "ggmfit" ={
           ff<-ggmfit (object$datainfo$S, n.obs=object$datainfo$n.obs, glist=object$glist,
                       start=start, details=0,...)
         },
         "ggmfitr"={
           ff<-ggmfitr(object$datainfo$S, n.obs=object$datainfo$n.obs, glist=object$glist,
                       start=start, details=0,...)
         }
         )

  idev  <-  ff$n.obs * (log(ff$detK) + sum(log(diag(ff$S))))  ## ideviance to independence model  
  idim      <-  ff$nvar 
  sat.dim   <-  ((idim+1)*idim) / 2
  dim.unadj <-  sat.dim - ff$df

  idf       <-  (dim.unadj-idim)
  logL.sat  <-  ff$logL + ff$dev/2

  aic       <-  -2*ff$logL + 2*dim.unadj
  bic       <-  -2*ff$logL + log(ff$n.obs)*dim.unadj

  dimension <- c(mod.dim=dim.unadj, sat.dim=sat.dim, i.dim=idim,df=ff$df,idf=idf)
  
  ans   <- list(dev=ff$dev, ideviance=idev, logL.sat=logL.sat,
                aic=aic, bic=bic,
                dimension=dimension
                )

  ff$S <- ff$n.obs <- ff$dev <- ff$df <- NULL
  ans <- c(ff,ans)
  
  object$fitinfo  <- ans
  object$isFitted <- TRUE
  class(object)   <- c("cModel","iModel")
  object
}































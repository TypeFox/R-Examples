##########################################################
##
## Discrete interaction model (log-linear model)
##
##########################################################

## FIXME: Delete df from loglin() output.

dmod <- function(formula, data, marginal=NULL, interactions=NULL, fit=TRUE, details=0){

  if (!inherits(data, c("data.frame","table"))){
    stop("data must be either a dataframe or a table \n")
  }
  if (inherits(data,"data.frame")){
    is_data.frame <- TRUE
    varNames <- names(data)
  } else {
    is_data.frame <- FALSE
    varNames <- names(dimnames(data))
  }

  if (length(marginal)>0){
    zzz  <- unlist(lapply(marginal, pmatch, varNames))
    zzz  <- zzz[!is.na(zzz)]
    marginal <- varNames[zzz]
  }

  ans <- .pFormula2(formula, varNames, marginal, interactions)
  
  if (is_data.frame){
      data <- xtabs(~., data=data[ans$varNames])
  } else {    
    if (length(ans$varNames) != length(varNames)){
      data         <- as.table(tableMargin(data, ans$varNames))
    }
  }
  varNames     <- names(dimnames(data))
   
  res <- list(glist     = ans$glist,
              varNames  = varNames,
              datainfo  = list(data=data),
              fitinfo   = NULL,
              isFitted  = FALSE
              )

  upd <- .dModel_finalize(ans$glist, varNames)
  res[names(upd)] <- upd

  class(res) <- c("dModel","iModel")

  if (fit){
    res <- fit(res)
  }
  res
}


.dModel_finalize<- function(glist, varNames){

  zzz <- isGSD_glist(glist)

  ## Requires data, glist
  ##
  glistNUM <- lapply(glist,
                     function(ll) {
                       match(ll, varNames)
                     })
  ret      <- list(glistNUM       = glistNUM,
                   isDecomposable = zzz[2],
                   isGraphical    = zzz[1]
                   )
##   cat(".dModel_finalize\n")
##   print(ret)
  return(ret)
}








fitted.dModel <- function(object,...){
  if (inherits(object,"fitted")){
    object$fitinfo$fit
  }
}


fit.dModel <- function(object, engine="loglin", print=FALSE, ...){

  ## FIXME: At some point we should allow for data in the form of a dataframe
  ##
  switch(engine,
         "loglin"={llfit<-loglin(object$datainfo$data, object$glist, fit=TRUE,
                                 print=print, ...)
                   names(llfit)[1] <- 'dev'  ## loglin() calls this slot 'lrt'
                 }
         )

  ## Calculate df's and adjusted df's
  ## Requires data, glist
  ##

  if (!is.null(object$glistNUM))
    glistNUM <- object$glistNUM
  else
    glistNUM <- lapply(object$glist, function(ll) {
      match(ll, object$varNames)
    })

  sat.dim.unadj   <- prod(dim(object$datainfo$data)) - 1
  sat.dim.adj     <- sum(object$datainfo$data>0) - 1
  ind.dim.unadj   <- sum(dim(object$datainfo$data)-1)

  if (object$isDecomposable){
    rr <- ripMAT(glist2adjMAT(object$glist))
    dim.adj   <- .loglinDecDim(rr$cliques, rr$separators, tableinfo=object$datainfo$data, adjust=TRUE)
    dim.unadj <- .loglinDecDim(rr$cliques, rr$separators, tableinfo=object$datainfo$data, adjust=FALSE)
  } else {
    dim.adj   <- NA
    dim.unadj <- .loglinGenDim(glistNUM, dim(object$datainfo$data))
  }

  indep.model <- loglin(object$datainfo$data, as.list(object$varNames),iter=1, print=FALSE)

  df.adj       <- sat.dim.adj   - dim.adj
  df.unadj     <- sat.dim.unadj - dim.unadj


  ideviance    <-  -(llfit$dev-indep.model$lrt)
  idf          <-  -(llfit$df-indep.model$df)

  iii <- object$datainfo$data * llfit$fit > 0

  llfit$logL      <- sum(object$datainfo$data[iii] * log(llfit$fit[iii]/sum(llfit$fit)))
  llfit$ideviance <- ideviance



  extra1      <- list(dim.unadj = dim.unadj,
                      dim.adj   = dim.adj,
                      df.unadj  = df.unadj,
                      df.adj    = df.adj,
                      idf       = idf,
                      ideviance = ideviance)

  dimension <- c(mod.dim=dim.unadj, sat.dim=sat.dim.unadj, i.dim=ind.dim.unadj, df=df.unadj, idf=idf,
                 mod.dim.adj = dim.adj,
                 sat.dim.adj = sat.dim.adj,
                 df.adj      = df.adj  )


  llfit$df        <- NULL
  
  llfit$aic       <- -2*llfit$logL + 2*dimension['mod.dim']
  llfit$bic       <- -2*llfit$logL + log(sum(object$datainfo$data))*dimension['mod.dim']


  ## Calculate warning codes for sparsity
  sparseinfo <- .dModel_sparsity (llfit, object$isDecomposable,
                                  object$glist,
                                  object$datainfo$data)

  fitinfo <- llfit
  fitinfo$sparseinfo <- sparseinfo
  fitinfo$dimension  <- dimension

  object$fitinfo  <- fitinfo

  object$isFitted <- TRUE

  class(object)   <- c("dModel","iModel")
  object
}

.dModel_sparsity <- function(llfit, isDecomposable, glist, data){

    ## Calculate warning codes for sparsity
  ##
  fc       <- as.numeric(llfit$fit)
  all.gt.0 <- sum( fc < 0.00001 ) == 0
  all.gt.5 <- sum( fc < 5 ) == 0

  frac.gt.0 <- sum(fc>0)/length(fc)
  frac.gt.5 <- sum(fc>5)/length(fc)

  if (!isDecomposable)
    {
      ## cat ("Model is not decompsable...\n")
      if (all.gt.0){
        df.ok        <- TRUE
        sparse.df.ok <- FALSE
        all.gt.5 <- sum( fc < 5 ) == 0
        if (!all.gt.5){
          #cat("Warning: table is sparse and asymptotic chi2 distribution is questionable (2)\n")
          chi2.ok <- FALSE
        } else {
          chi2.ok <- TRUE
        }
      } else {
        #cat("degrees of freedom are unadjusted and can not be trusted (1)\n")
        #cat("Warning: asymptotic chi2 distribution is questionable (2)\n")
        df.ok        <- FALSE
        sparse.df.ok <- FALSE
        chi2.ok      <- FALSE
      }
    }
  else
    {
      ## cat("Model is decompsable...\n")
      if (all.gt.0){
        df.ok         <- TRUE
        sparse.df.ok  <- TRUE
      } else {
        df.ok         <- FALSE
        sparse.df.ok  <- TRUE
      }

      all.gt.5 <- 1 - sum( fc < 5 ) > 0
      if (all.gt.5){
        ## It is a dense table
        chi2.ok <- TRUE
      } else {
        ## It is not a dense table, but the table may be dense in the cliques
        ## Then find, in each clique, those marginal cells with positive counts. This leads to the df adjustment.
        ## For each marginal cell with positive counts, check if counts are > 5. If so, the chi2 is ok.
        sparsecode3 <- rep.int(0, length(glist))
        for (ii in 1:length(glist)){
          tmC <- tableMargin(data, glist[[ii]])
          tm0 <- tmC>0
          tm5 <- tmC[tm0]>5
          if (length(tm5) == length(tm0)){
            sparsecode3[ii] <- 1
          }
        }
        all.gt.5.when.gt.0 <- sum(sparsecode3)==length(glist)

        if (all.gt.5.when.gt.0){
          #cat("Warning: table is sparse and degrees of freedom have been adjusted to reflect sparsity of table (3)\n")
          chi2.ok <- TRUE
        } else {
          #cat("Warning: table is sparse and asymptotic chi2 distribution is questionable (4)\n")
          chi2.ok <- FALSE
        }
      }
    }

  sparseinfo <-
    c(chi2.ok=chi2.ok, df.ok=df.ok, sparse.df.ok=sparse.df.ok)

  sparseinfo
}




























  ## Possibly take out a slice of the table
##   if (!missing(context)){
##     connames  <- context[[1]]
##     conlevels <- context[[2]]
##     data <- as.table(tableSlice(data, margin=connames, level=conlevels))
##   } else {
##     connames <- conlevels <- NULL
##   }





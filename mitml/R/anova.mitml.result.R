anova.mitml.result <- function(object, ...){

  # *** select method for object class
  #
  cls <- class(object[[1]])

  # default
  lr.method <- "default"
  # lm 
  if(cls[1]=="lm") lr.method <- "lm"
  # merMod (lme4)
  if(any(grepl("^l?merMod$",cls[1]))){
    if(!requireNamespace("lme4", quietly=TRUE)) stop("The 'lme4' package must be installed to handle 'merMod' class objects.")
    lr.method <- "lmer"
  }
  # lme (nlme)
  if(any(grepl("^.?lme$",cls[1]))){
    if(!requireNamespace("nlme", quietly=TRUE)) stop("The 'nlme' package must be installed to handle 'lme' class objects.")
    lr.method <- "nlme"
  }

  # check use of D3
  if(lr.method%in%c("lm","lmer","nlme")){
    method <- "D3"
  }else{
    warning("The 'D3' method is currently not supported for models of class '",cls[1],"'. Switching to 'D2'.")
    method <- "D2"
  }
  
  # *** testModels
  #

  modlist <- list(object,...)

  # order models
  nm <- length(modlist)
  df <- integer(nm)
  for(mm in 1:nm){ 
    df[mm] <- switch(lr.method,
      lmer=attr(logLik(modlist[[mm]][[1]]),"df"),
      nlme=attr(logLik(modlist[[mm]][[1]]),"df"),
      lm=.tryResidualDf(modlist[[mm]][[1]]),
      default=.tryResidualDf(modlist[[mm]][[1]])
    )
  }
  if(lr.method%in%c("lm","default")) df <- abs(df-max(df))
  modlist <- modlist[ order(df,decreasing=T) ]

  # stepwise comparison
  outlist <- as.list(1:(nm-1))
  for(mm in 1:(nm-1)){
    if(method=="D3"){
      outlist[[mm]] <- testModels(model=modlist[[mm]], null.model=modlist[[mm+1]],
                                  method="D3")
    }else{
      outlist[[mm]] <- testModels(model=modlist[[mm]], null.model=modlist[[mm+1]],
                                  method="D2", use="likelihood")
    }
  }

  # get model formulas
  fml <- character(nm)
  for(mm in 1:nm){
    fml[mm] <- tryCatch(
      gsub("[[:space:]]", "", Reduce(paste, deparse( formula(modlist[[mm]][[1]]) ))),
      error=function(f) NULL
    )
  }

  # check for REML
  reml <- any(sapply(outlist, function(x) x$reml))

  out <- list(
    call=match.call(),
    test=outlist,
    formula=fml,
    method=method,
    use="likelihood",
    reml=reml
  )

  class(out) <- "mitml.anova"
  out

}


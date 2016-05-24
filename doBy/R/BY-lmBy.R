
lmBy <- function(formula, data, id=NULL, ...){
  cl   <- match.call()
  mff  <- parseGroupFormula(formula)
  groupData  <- splitBy(mff$groupFormula, data=data)

  mmm <- mff$model
  mm  <- lapply(groupData, function(wd) {
    zzz<-lm(mmm, data=wd, ...)
    zzz$call[[2]]<- mmm
    zzz
  })

  if (!is.null(id)){
    id.vars <- unique(c(all.vars(id), all.vars(mff$groupFormula)))
  } else {
    id.vars <- all.vars(mff$groupFormula)
  }
  id.data <- do.call(rbind, lapply(groupData, function(wd) {wd[1,id.vars,drop=FALSE]}))

  attr(mm,  "call")     <- cl
  attr(mm,  "dataList") <- groupData
  attr(mm,  "idData")   <- id.data	
  
  class(mm) <- "lmBy"
  mm
}


print.lmBy <- function(x, ...){
  ##lapply(c(x), print)
  print(c(x))
  return(invisible(x))
}


summary.lmBy <- function(object, ...){
  res <- lapply(object, summary)
  class(res) <- "summary_lmBy"
  res
}

print.summary_lmBy <- function(x, ...){
  lapply(x, print)
  return(invisible(x))
}

coef.summary_lmBy <- function(object, simplify=FALSE, ...){
  ans <- lapply(object, coef)
  if (simplify){
    cc <- do.call(rbind, ans)
    cn <- colnames(cc)
    rn <- rownames(cc)

    
    nn <- names(ans)
    rn <- rownames(ans[[1]])
    ff <- factor(rep(nn, each=length(rn)))
    
    rownames(cc) <- NULL
    ans <- data.frame(ff, rn, as.data.frame(cc))
    colnames(ans) <- c("stratum", "parameter", cn)
  }
  ans
}


getBy <- function(object, name=c()){
  if (missing(name)) 
    stop("'name' must not be missing")
  switch(class(object),
         "lmBy"={
           ii <- match(name, c("dataList","idData"))
           if (is.na(ii))
             stop(sprintf("%s not available", name))
           attr(object,name)	
         })
  
}

coef.lmBy <- function(object, augment=FALSE, ...){
  ans <- do.call(rbind, lapply(object, coef))
  if (augment){
    ans <- cbind(ans, getBy(object,"idData"))
  }
  ans
}


fitted.lmBy <- function(object, augment=FALSE,...){
  ans <- lapply(object, fitted)
  if (augment) {
    ans <- mapply(function(a,b){data.frame(.fit=a,b)}, ans, getBy(object, "dataList"),
                  SIMPLIFY=FALSE)
  }
  ans
}


residuals.lmBy <- function(object, augment=FALSE,...){
  ans <- lapply(object, residuals)
  if (augment) {
    ans <- mapply(function(a,b){data.frame(.fit=a,b)}, ans, getBy(object, "dataList"),
                  SIMPLIFY=FALSE)
  }
  ans
}




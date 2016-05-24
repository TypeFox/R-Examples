parseGroupFormula <- function(form)
{
  if (class(form) != "formula" || length(form) != 3)
    stop("formula must be a two-sided formula object")
  rhs <- form[[3]]
  if (class(rhs) != "call" || rhs[[1]] != as.symbol('|'))
    stop("rhs of formula must be a conditioning expression")
  form[[3]] <- rhs[[2]]
  groups <- rhs[[3]]
  grpFormula <- as.formula(paste("~", deparse(groups)))
  list(model = form, groups = groups, groupFormula=grpFormula)
}

## createFunBy <- function(formula, data, FUN, class=NULL, ...){
##   mm <- parseGroupFormula(formula)
##   groupData <- splitBy(mm$groupFormula, data=data)
##   res <- lapply(groupData, function(xx) FUN(mm$model, data=xx, ...))
##   ##attr(res, "splitByData") <- attributes(groupData)
##   if (!is.null(class))
##     class(res) <- class
##   res
## }

formulaFunBy <- function(formula, group, data, FUN, class=NULL, ...){
  grpData <- splitBy( group, data=data )
  res <- lapply( grpData, function(.xx){ FUN(formula, data=.xx, ...)})
  if (!is.null(class))
    class(res) <- class
  res
}

xyFunBy <- function(xy, group, data, FUN, class=NULL, ...){
  if ( !(class(xy) %in% c("character", "formula") ))
    stop(" 'xy' must be a character vector or a two-sided formula ")
  if (class(xy)=="formula"){
    vars <- all.vars(xy)
  } else {
    vars = xy
  }
  grpData <- splitBy( group, data=data )
  if (length(vars)==2){
    res <- lapply(grpData, function(dd){
      x <- dd[,vars[1]]
      y <- dd[,vars[2]]
      FUN(x, y, ...)})
  }
  else{
    res <- lapply(grpData, function(dd){
      x <- dd[,vars[1]]
      FUN(x, ...)})	
  }
  if (!is.null(class))
    class(res) <- class
  res
}





## funBy <- function(formula, data, FUN, ...){

##   cl <- match.call()
##   ep <- 1
##   lhs      <- formula[[2]]
##   lh.var   <- .lhsParse(lhs)
##   if (length(lh.var)==1)
##     lh.var <- list(lh.var)
##   lh.var   <- paste(unlist(lh.var))	
##   rh.var   <- .rhsVar(formula)
  
##   data.var <- names(data)
##   cls      <- lapply(data, class)
  
##   if (length(lh.var)==1 && lh.var=="."){		
##     num.var   <- data.var[  cls %in% c("numeric","integer")]
##     lh.var    <- setdiff(num.var, rh.var)
##     #cat(sprintf("deriving numeric lh.var: %s\n", toString(lh.var))) 
##     if (length(lh.var)==0)
##       stop("There are no numeric variables to scale")
##   } else {
##     #cat(sprintf("using specified numeric lh.var: %s\n", toString(lh.var))) 
##   }

##   print("here")  
##   if (length(rh.var)==1 && rh.var=="."){		
##     fac.var <- data.var[!(cls %in% c("numeric","integer"))]
##     if (length(fac.var)>0){
##       if(ep>=1) cat(sprintf("stratifying by non-numerics fac.var: %s\n", toString(fac.var)))
##       rh.var    <- fac.var 
##       rh.data   <- data[, rh.var, drop=FALSE]
##       rh.string <- apply(rh.data, 1, paste, collapse="@")
##       rh.unique <- unique(rh.string)
##     } else {
##       if(ep>=1) cat(sprintf("no non--numerics to stratify by\n"))
##       rh.string <- rep.int(1, NROW(data))
##       rh.unique <- 1 
##     }
##   } else {
##     if(ep>=1) cat(sprintf("stratifying by specified rh.var: %s\n", toString(rh.var)))    
##     rh.data   <- data[, rh.var, drop=FALSE]
##     rh.string <- apply(rh.data, 1, paste, collapse="@")
##     rh.unique <- unique(rh.string)
##   }


##   lh.idx   <- match(lh.var, names(data))
##   lh.data  <- as.matrix(data[,lh.var,drop=FALSE])

##   print(lh.data)
##   for (ii in seq_along(rh.unique)){	
##     idx  <- rh.unique[ii] == rh.string

##     sc   <- FUN(lh.data[idx,,drop=F], ...)
##     lh.data[idx,] <- sc
##   }
  
##   data[,lh.idx] <- lh.data
##   data
## }	


subsetBy <- function(formula, subset, data=parent.frame(), select, drop=FALSE, join=TRUE,...){
  ddd<-splitBy(formula, data=data)
  subsetMissing <- missing(subset)
  selectMissing <- missing(select)  
  e <- substitute(subset)
  ddd<-lapply(ddd, 
    function(x){
      if (subsetMissing) 
          r <- TRUE
      else {
          r <- eval(e, x, parent.frame())
          if (!is.logical(r)) 
              stop("'subset' must evaluate to logical")
          r <- r & !is.na(r)
      }
      if (selectMissing) 
          vars <- TRUE
      else {
          nl <- as.list(1:ncol(x))
          names(nl) <- names(x)
          vars <- eval(substitute(select), nl, parent.frame())
      }
      x[r, vars, drop = drop]
    }
  )
  if (join)
    return(do.call("rbind",ddd))
  else
    return(ddd)
}


## subsetBy2 <- function(formula, x=parent.frame(), subset, select, drop=FALSE, join=TRUE, ...){

##   cl <- match.call()
##   cl[[2]] <- NULL
##   cl[[1]] <- as.name("subset")
##   cl[[2]] <- as.name("_x")
##   cl[[match("join", names(cl))]] <- NULL
  

##   xlist <-splitBy(formula, data=x)
##   ans <- lapply(xlist, function(xx){eval(cl, envir=list(`_x`=xx), parent.frame()) })
##   if (join)
##     return(do.call("rbind",ans))
##   else
##     return(ans)

## }

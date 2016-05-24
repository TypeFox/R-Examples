splitBy <-function (formula, data = parent.frame(), drop=TRUE) { #, return.matrix=FALSE){

  data.var  <- names(data)

  if ( !(class(formula) %in% c("formula", "character")) ){
    stop("'formula' must be a right hand sided formula or a character vector")
  }

  if ( class(formula) == "formula" ){
    rhs       <- formula[[2]]
    rhs.var   <- all.vars( rhs )
  } else {
    rhs.var <- formula
  }

  cls<-lapply( data, class )
  factor.columns <- which(unlist(lapply(cls, function(x) any( x %in% "factor"))))
  #print(factor.columns)


  cls       <- lapply(data, class)
  num.idx   <- cls %in% c("numeric","integer")
  fac.var   <- data.var[ !num.idx ]

  rhs.fac   <- intersect( rhs.var, data.var )
  if ("." %in% rhs.var){ ## need all factors not mentioned elsewhere as grouping factors
    rhs.fac <- union( fac.var, rhs.fac)
  }
  #str(list(rhs.var=rhs.var, rhs.fac=rhs.fac, fac.var=fac.var))

  rh.trivial <- length( rhs.var ) == 0 #; cat(sprintf("rh.trivial=%d\n", rh.trivial))

  ## FIXME rhs.fac, rhs.var -- clean up!!!
  ## Use: data, rhs.fac, rh.trivial
  if ( rh.trivial ){
    grps <- rep.int(1, nrow(data))
    unique.grps <- 1
    rh.idx    <- 1
  } else {
    grps <- .get_rhs_string( data, rhs.fac, sep.string="|")
    unique.grps <- unique(grps)
    rh.idx    <- match(unique.grps, grps)
  }


  groupData <- vector("list", length(unique.grps))
  names(groupData) <- unique.grps
  for (ii in 1:length(unique.grps)){
    dd <- data[unique.grps[ ii ] == grps,,drop=FALSE]
    if (drop && length(factor.columns)>0){
      for (jj in 1:length(factor.columns)){
        dd[, factor.columns[ jj ]] <- factor( dd[, factor.columns[ jj ]] )
      }
    }
    groupData[[ ii ]] <- dd
  }

  idxvec <- vector("list", length(unique.grps))
  names(idxvec) <- unique.grps
  for (ii in 1:length(unique.grps)){
    idxvec[[ii]] <- which(grps==unique.grps[ii])
  }

  groupid <- data[rh.idx, rhs.fac, drop=FALSE]
  rownames(groupid) <- 1:nrow(groupid)

  attr(groupData,"groupid") <- groupid
  attr(groupData,"idxvec")  <- idxvec
  attr(groupData,"grps")    <- grps

  class(groupData) <- c("splitByData", "list")

  groupData
}

print.splitByData <- function(x,...){
#  print(attr(x,"groupid"))
  print(cbind(listentry=names(x), attr(x,"groupid")))
  return(invisible(x))
}

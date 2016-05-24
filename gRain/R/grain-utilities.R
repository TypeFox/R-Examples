randomCPT <- function(object, states=c("yes","no")){

    if (class(object) != "graphNEL")
        stop("'object' must be a graphNEL object\n")
    if (!is.DAG(object))
        stop("'object' is not a DAG\n")

    vpa <- vpar( object )
    n.states  <- length(states)
    cpt <- lapply(vpa, function(zz)
                  cptable(zz, values=runif( n.states^length(zz) ), levels=states))

    compileCPT( cpt )
}



.formula2char <- function(f) {
	unlist(rhsf2list(f))
}


.namesDimnames <- function(x)
    names(dimnames(x))


setSliceValue <- function(x, slice, complement=FALSE, value=0){
    margin <- names(slice)
    level  <- unlist(slice, use.names=FALSE)
    idx <- tableGetSliceIndex(x, margin = margin, level = level,
                              complement = complement)
    x[idx] <- value
    x
}




getgrain <- function(object, name=c("universe", "data", "dag", "ug", "cptlist",
                                 "origpot", "temppot", "equipot", "rip",
                                 ## "isInitialized",
                                 "isCompiled", "isPropagated",
                                 "evidence", "pEvidence",
                                 "control", "details")){

    switch(name,
           universe 		    = object$universe,
           data 			      = object$data,
           dag 				      = object$dag,
           ug 			    	  = object$ug,
           cptlist          = object$cptlist,

           origpot          = object$origpot,
           temppot          = object$temppot,
           equipot          = object$equipot,
           rip              = object$rip,

           ## isInitialised    = object$isInitialized,
           isCompiled       = object$isCompiled,
           isPropagated     = object$isPropagated,

           evidence         = object$evidence,
           pEvidence        = object$pEvidence,
           control          = object$control,
           details          = object$details
           )
}


.infoPrint <- function(details, limit=1, ...,  prefix='.'){
  if(details>=limit){
    cat(paste(rep(prefix, limit), collapse=''), "")
    cat(...)
  }
}

.infoPrint <- function(details, limit=1, ...,  prefix='.'){}


.infoPrint2 <- function(details, limit=1, fmt, ...,  prefix='.'){
  if (details>=limit)
    cat(paste(paste(rep(prefix, limit), collapse=''),  sprintf(fmt, ...), collapse=' '))
}

.colstr <- function(x, collapse=" ")
  paste(x, collapse=collapse)



.printList <- function(x){
  mapply(function(xx,ii) cat(" ",ii,paste(xx, collapse=' '),"\n"), x, 1:length(x))
  return()
}


printlist <- function(x,d=0) UseMethod("printlist")

printlist.default <- function(x,d=0){
  paste("(", paste(x,collapse=' '),")",sep='')
}

printlist.list <- function(x,d=0){
  tmp     <- unlist(lapply(x, printlist, d+2),recursive=FALSE)
  prefix  <- as.list(c("(",rep(" ",length(tmp)-1)))
  posfix  <- as.list(c(rep(" ",length(tmp)-1),")"))
  as.list(mapply(function(l,x,r) {paste(l,x,r,sep='')}, prefix, tmp, posfix))
}

splitVec <- function(val, lev) UseMethod("splitVec")

splitVec.default <- function(val, lev){
  m    <- matrix(val,ncol=lev)
  cval <- unlist(apply(m,2,list),recursive=FALSE)
  cval
}

splitVec.list <- function(val, lev){
  lapply(val, splitVec, lev)
}



## ###############################################################
##
## as.grain() methods
##
## added July 2014
##
## ###############################################################

## as.grain <- function(x, ...){
##     UseMethod("as.grain")
## }

## as.grain.CPTspec <- function(x, ...){
##     grain( x )
## }

## as.grain.graphNEL <- function(x, data, smooth=0, ...){
##     if (missing(data))
##         stop("'data' needed to create network from graph")
##     grain(x, data=data, smooth=smooth)
## }

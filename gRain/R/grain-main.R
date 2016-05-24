### #####################################################
###
### Creating grain objects
###
### #####################################################

grain <- function(x, data=NULL, control=list(), smooth=0, details=0,...){
  UseMethod("grain")
}

## A list of cpt's
##
grain.CPTspec <- function(x, data=NULL, control=list(), smooth=0, details=0,...){
  ##cat("grain.CPTspec\n")
  control  <- .setControl(control)
  ans  <- c(list(universe    = attr(x,"universe"),
                 data        = data,
                 dag         = attr(x,"dag"),    ## Needed to save network in Hugin format
                 cptlist     = c(x)              ## Needed to save network in Hugin format
                 ),
            .setExtraComponents(control, details))

  class(ans) <- c("CPTgrain","grain")
  return(ans)
}

grain.POTspec <- function(x, data=NULL, control=list(), smooth=0, details=0,...){
  ## cat("grain.POTspec\n")
  control  <- .setControl(control)
  ans  <- c(list(universe    = attr(x, "universe"),
                 data        = data,
                 equipot     = c(x),
                 ug          = attr(x, "ug"),
                 rip         = attr(x, "rip"),
                 dag         = attr(x, "dag"),    ## Needed to save network in Hugin format
                 cptlist     = attr(x, "cptlist") ## Needed to save network in Hugin format
                 ),
            .setExtraComponents(control, details))
  class(ans) <- c("POTgrain","grain")
  ans
}

## A graph + data (wrappers for calling grain.POTspec and grain.CPTspec)
grain.graphNEL <- function(x, data=NULL, control=list(), smooth=0, details=0,...){
    if (missing(data))
        stop("Data must be given to create grain from graph\n")
    if (!(is.array(data) || is.data.frame(data)))
        stop("Data must be an array or a dataframe\n")

  if (is.DAG(x)){
    ans <- grain(compileCPT(extractCPT(data, x, smooth=smooth)),
                 data=data, control=control, details=details)
  } else {
    if (is.TUG(x)){
      ans <- grain(compilePOT(extractPOT(data, x, smooth=smooth)),
                   data=data, control=control, details=details)
    } else {
      stop("graph 'x' is neither a directed acyclic graph or a triangulated undirected graph")
    }
  }
  return(ans)
}

grain.dModel <- function(x, data=NULL, control=list(), smooth=0, details=0,...){

    if (!x$isDecomposable)
        stop("Model must be decompsable graphical model\n")
    g <- ugList( terms( x ) )
    if (is.null(data))
        data <- x$datainfo$data
    grain(g, data=data, smooth=smooth, details=details, ...)
}


## Printing grain
##
print.grain <- function(x,...){
    cat("Independence network: Compiled:", x$isCompiled,
        "Propagated:", x$isPropagated, "\n")
    cat("  Nodes:"); str(unname(nodeNames(x)))
    if ( !is.null(x$evidence) ){
        cat("  Evidence:\n");
        print(as.data.frame( getEvidence(x)$summary) )
        if (!is.null((p <- pEvidence(x))))
            cat(sprintf("  pEvidence: %f\n", p))
    }
    return(invisible(x))
}

.setControl <- function(control){
  con <- list(timing=0)
  con[(namc <- names(control))] <- control
  con
}

.setExtraComponents <- function(control, details){
  list(
      ## FIXME Do we isInitialized? I doubt! Status: Now removed!
      ## isInitialized = FALSE,
      isCompiled    = FALSE,
      isPropagated  = FALSE,
      evidence      = NULL,
      control       = .setControl(control),
      details       = details
      )
}





## NOTICE:
##
## extractPOT() generates
## {p(C1), p(R2|S1), ..., p(Rn|Sn)}
## so these are clique potentials but they are not equilibrated.
## extractPOT() also generates  {p(v|pa(v)} which is not needed except to be
## able to save a network as a hugin net file.
##
## extractCPT() generates
## {p(v|pa(v))}
##
## compilePOT() and compileCPT() only makes 'internal' computations/setups
## and do not fundamentally change the above.
##


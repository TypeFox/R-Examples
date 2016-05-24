#' Verify input parameters for the graph
#'
#' Mainly for internal use.
#'
#' @param coord Coordinates of the locations
#' @param type Type of graph
#' @param par Parameter(s) for the graph
#' @param maxR Maximum range for edges, helps in large patterns.
#' @param doDists Precompute distances? Speeds up some graphs, takes up memory.
#' @param preGraph Precomputed graph, taken as a super-graph
#'
#' @export

sg_verify_parameters <- function(coord, type, par, maxR, doDists, preGraph) {
  # which graph
  i<-pmatch(type, names(SG_GRAPH_PARAMETERS))
  if(is.na(i))
    stop(paste0("'", type, "' is not a supported graph type. Pick one from:\n",
                             paste(names(SG_GRAPH_PARAMETERS), collapse=" ") )
           )

  # check par
  ####################################################################
  # graphs with some parameters
  #
  # special cases with marks:
  par_ok <-  TRUE
  if(i %in% c(3,4)) { # mass geometric or mark cross
    if(!is.numeric(par) | length(par) != nrow(coord)){
      par_should_be <- c(marks=paste0("vector of length ", nrow(coord),  " for the marks"))
      par_ok <- FALSE
    }
  }
  # RST needs coordinates
  else if(i == 8) {
    if(!is.numeric(par) | length(par) != ncol(coord)){
      par_should_be <- paste0("vector of length ", ncol(coord),  " for the center point coordinates.")
      par_ok <- FALSE
    }
  }
  # CCC
  else if(i==10){
    if( !is.factor(par) | length(par) != nrow(coord)){
      par_should_be <- paste0("vector of length ", nrow(coord),  " for the center point coordinates.")
      par_ok <- FALSE
    }
    else{
      par <- as.integer(par)
    }
  }
  else{
    ################
    # no marks
    par <- as.numeric(par)
    par_should_be <- unlist( SG_GRAPH_PARAMETERS[[i]] )
    par_ok <- length(par) == length(par_should_be)
    if(length(par) & length(par_should_be)==0) stop(paste0("'", type, "' does not take parameters."))
  }

  ###############
  if(!par_ok)
    stop(paste0("'", type, "' graph needs par=",  paste(par_should_be, collapse=","))
  )
  # check maxR
  if(maxR<0) stop("'maxR' < 0 given.")
  if(maxR>0){
    # not all support pre-R
    Rsup <- c(2,9)
    if(! i %in% Rsup) stop("Graphs that at the moment support maxR: ", paste0(names(SG_GRAPH_PARAMETERS)[Rsup], collapse=", "))
  }
  # check preGraph
  if(!is.null(preGraph)){
    if(!is(preGraph, "sg")) stop("preGraph not of class 'sg'.")
    if(preGraph$N != nrow(coord)) stop("preGraph and pattern are of different size.")
  }
  # clash
  if(maxR>0 & !is.null(preGraph)) stop("Can not handle both 'maxR' and 'preGraph'. Perhaps use cut.sg on preGraph?")
  ##############
  # ok
  par
}

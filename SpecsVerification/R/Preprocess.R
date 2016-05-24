.PreprocessEns <- function(x) {
  # check class of x (note that it can have multiple classes, hence any(...))
  if(!(any(class(x) %in% c("matrix", "data.frame", "numeric","logical")))) {
    stop("Can only handle ensembles of class matrix, data.frame, or numeric.")
  }
  # interpret vector-valued ensemble as one K-member ensemble, i.e. one row
  if (is.null(dim(x))) {
    x <- matrix(x, nrow=1)
  } else {
    x <- as.matrix(x)
  }
  x[is.infinite(x)] <- NA
  return(x)
}
.PreprocessObs <- function(x) {
  # check class of x (note that can have multiple classes)
  cl.x <- class(x)
  if ("integer" %in% cl.x) {
    x <- as.numeric(x)
    cl.x <- "numeric"
  }
  if(!(any(cl.x %in% c("matrix", "data.frame", "numeric","logical")))) {
    stop("Can only handle observations of class matrix, data.frame, or numeric.")
  }
  # collapse into a vector
  x <- drop(as.matrix(x))
  x[is.infinite(x)] <- NA
  return(x)
}
.RemoveNaMembers <- function(x) {
  na.cols <- apply(x, 2, function(z) all(is.na(z)))
  return(x[, !na.cols, drop=FALSE])
}



Preprocess <- function(ens=NULL, ens.ref=NULL, obs=NULL) {

  # logical vector indicating which of ens, ens.ref, obs was provided as an argument
  was.provided <- sapply(list(ens=ens, ens.ref=ens.ref, obs=obs), 
                  function(x) !(is.null(x)))

  ### transform the ensembles to matrices or leave as NA ###
  if (was.provided["ens"]) {
    ens <- .PreprocessEns(ens)
  }
  if (was.provided["ens.ref"]) {
    ens.ref <- .PreprocessEns(ens.ref)
  }

  ### transform observation to vector or leave as NA ### 
  if (was.provided["obs"]) {
    obs <- .PreprocessObs(obs)
  }
  
  # if more than one argument is provided, check for equal length of time
  # dimension and stop if lenghts are different
  if (sum(was.provided) > 1) {
    N.vec <- c(ens=NA, ens.ref=NA, obs=NA)
    if (was.provided["ens"]) N.vec["ens"] <- nrow(ens)
    if (was.provided["ens.ref"]) N.vec["ens.ref"] <- nrow(ens.ref)
    if (was.provided["obs"]) N.vec["obs"] <- length(obs)
       
    # the `as.matrix` below is required to avoid errors if ens == NA
    if(length(unique(N.vec[was.provided])) != 1) {
      stop("Inputs do not have equal time dimensions.")
    }
  }

  # erase ensemble members that are all NA
  if (was.provided["ens"]) {
    ens <- .RemoveNaMembers(ens)
  }
  if (was.provided["ens.ref"]) {
    ens.ref <- .RemoveNaMembers(ens.ref)
  }

  # return
  ret <- list(ens=ens, ens.ref=ens.ref, obs=obs)
  return(ret)

}


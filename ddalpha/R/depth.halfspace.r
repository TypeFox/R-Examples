################################################################################
# File:             depth.halfspace.r
# Created by:       Pavlo Mozharovskyi
# First published:  28.02.2013
# Last revised:     18.06.2015
# 
# Computation of the Tukey data depth.
################################################################################

.parse_HSD_pars <- function(exact, method){
  if(missing(exact) && missing(method))
    return(0)
  
  if(!missing(exact)){
    if(exact == F){
      if (missing(method))
        return(0)
      else if (method != 0 && method != "Sunif.1D")
        stop("Wrong combination of 'exact' and 'method' parameters.")
    }
    else{
      if (missing(method))
        return(1)           # default exact
      else
        if (!(method %in% 1:3 || method %in% c("recursive","plane","line")))
          stop("Wrong combination of 'exact' and 'method' parameters.")
    }
  }
  
  if (!(method %in% 0:3 || method %in% c("Sunif.1D","recursive","plane","line")))
    stop("Wrong parameter 'method'.")
      
  if (is.character(method))
    method =  switch (method, random = 0, recursive = 1, plane = 2, line = 3)
  
  return(method)
}

depth.halfspace <- function(x, data, exact, method, num.directions = 1000, seed = 0){
  if (!is.matrix(x) 
      && is.vector(x)){
    x <- matrix(x, nrow=1)
  }
  if (!(is.matrix(data) && is.numeric(data)
        || is.data.frame(data) && prod(sapply(data, is.numeric))) 
      || ncol(data) < 2){
    stop("Argument \"data\" should be a numeric matrix of at least 2-dimensional data")
  }
  if (!is.numeric(x)){
    stop("Argument \"x\" should be numeric")
  }
  if (ncol(x) != ncol(data)){
    stop("Dimensions of the arguments \"x\" and \"data\" should coincide")
  }
  if (ncol(data) + 1 > nrow(data)){ #?
    stop("To few data points")
  }

  points <- as.vector(t(data))
  objects <- as.vector(t(x))
  
  method = .parse_HSD_pars(exact, method)
  
  if (method == 0)
  if (!is.numeric(num.directions) 
      || is.na(num.directions) 
      || length(num.directions) != 1 
      || !.is.wholenumber(num.directions) 
      || !(num.directions > 1 && num.directions < 10000000)){
    numDirections <- 1000
    warning("Argument \"num.directions\" not specified correctly. 1000 is used as a default value")
  }else{
    numDirections <- num.directions
  }
  
  if (method == 0){
    c <- as.vector(nrow(data))
    k <- numDirections
    ds <- .C("HDepth", 
             as.double(points), 
             as.double(objects), 
             as.integer(nrow(x)), 
             as.integer(ncol(data)), 
             as.integer(c), 
             as.integer(1), 
             as.double(0), 
             as.double(0), 
             as.integer(k), 
             as.integer(1), # use the same directions and projections
             as.integer(seed),
             depths=double(nrow(x)))$depths
  } else 
    if (method %in% 1:3){
      ds <- .C("HDepthEx", 
               as.double(points), 
               as.double(objects), 
               as.integer(nrow(data)),
               as.integer(nrow(x)),  
               as.integer(ncol(data)), 
               as.integer(method), 
               depths=double(nrow(x)))$depths  
    }
  else 
    stop("wrong choise of the algorithm, method = ", method)
  
  return (ds)
}




depth.space.halfspace <- function(data, cardinalities, exact, method, num.directions = 1000, seed = 0){
  if (seed != 0) set.seed(seed)
  if (!(is.matrix(data) && is.numeric(data)
        || is.data.frame(data) && prod(sapply(data, is.numeric))) 
      || ncol(data) < 2){
    stop("Argument \"data\" should be a numeric matrix of at least 2-dimensional data")
  }
  if (!is.vector(cardinalities, mode = "numeric") 
      || is.na(min(cardinalities)) 
      || sum(.is.wholenumber(cardinalities)) != length(cardinalities) 
      || min(cardinalities) <= 0 
      || sum(cardinalities) != nrow(data)){
    stop("Argument \"cardinalities\" should be a vector of cardinalities of the classes in \"data\" ")
  }
  if (sum(cardinalities < ncol(data) + 1) != 0){
    stop("Not in all classes sufficiently enough objetcs")
  }
  if (!is.numeric(num.directions) 
      || is.na(num.directions) 
      || length(num.directions) != 1 
      || !.is.wholenumber(num.directions) 
      || !(num.directions > 1 && num.directions < 10000000)){
    numDirections <- 1000
    warning("Argument \"num.directions\" not specified correctly. 1000 is used as a default value")
  }else{
    numDirections <- num.directions
  }
  
  x <- as.vector(t(data))
  c <- as.vector(cardinalities)
  
  method = .parse_HSD_pars(exact, method)
  
  if (method == 0)
    if (!is.numeric(num.directions) 
        || is.na(num.directions) 
        || length(num.directions) != 1 
        || !.is.wholenumber(num.directions) 
        || !(num.directions > 1 && num.directions < 10000000)){
      numDirections <- 1000
      warning("Argument \"num.directions\" not specified correctly. 1000 is used as a default value")
    }else{
      numDirections <- num.directions
    }
  
  if (method == 0){
    k <- numDirections
    rez <- .C("HDSpace", 
              as.double(x), 
              as.integer(ncol(data)), 
              as.integer(c), 
              as.integer(length(cardinalities)), 
              as.integer(k), as.integer(1), 
              as.integer(seed),
              dspc=double(nrow(data)*length(cardinalities)), 
              dirs=double(k*ncol(data)), 
              prjs=double(k*nrow(data)))
    depth.space <- matrix(rez$dspc, nrow=nrow(data), ncol=length(cardinalities), byrow=TRUE)
  }else 
    if (method %in% 1:3){
      ds <- .C("HDepthSpaceEx", 
               as.double(x), 
               as.double(x), 
               as.integer(c), 
               as.integer(length(cardinalities)), 
               as.integer(nrow(data)),  
               as.integer(ncol(data)), 
               as.integer(method), 
               depths=double(nrow(data)*length(cardinalities)))$depths  
      
      depth.space <- matrix(ds, nrow=nrow(data), ncol=length(cardinalities), byrow=F)
    }
  else 
    stop("wrong choise of the algorithm, method = ", method)
  
  return (depth.space)
}


chooseDim <- function(array, dim, value, drop = FALSE) { 
  # Create list representing arguments supplied to [
  # bquote() creates an object corresponding to a missing argument
  dimnames <- attributes(array)$dimensions
  
  indices <- rep(list(bquote()), length(dim(array)))
  indices[[dim]] <- value
  
  if (dim(array)[dim] < max(value)) {
    stop('Chosen member exceeds the member range of the dataset.')
  }
  
  # Generate the call to [
  call <- as.call(c(
    list(as.name("["), quote(array)),
    indices,
    list(drop = drop)))
  # Print it, just to make it easier to see what's going on
  # Print(call)
  
  # Finally, evaluate it
  output <- eval(call)
  
  if (length(dim(output)) == length(dimnames)) {
    attributes(output)$dimensions <- dimnames
  } else if (length(dim(output)) < length(dimnames)){
    
    # In this case, one dimension is dropped, if value is a number 
    # and drop == T, this situation can appear. So the dropped dimemsion
    # should be the chosen dimension.
    i <- 1:length(dimnames)
    # get rid of the dropped dimensin
    i <- i[-dim]
    attributes(output)$dimensions <- dimnames[i]
  }
  
  return(output)
}


adjustDim <- function(data, ref = 'no') {
  # input data is an array
  # ref is the Data part of a hyfo file, used as reference
  # Further may be arranged into a seperate function
  # the input reference will be put at first, then the rest 
  if (is.null(data)) return(NULL)
  if (identical(ref, 'no')) {
    # Default
    refOrder <- c('lon', 'lat', 'time')
  } else if (is.character(ref)) {
    refOrder <- ref
  } else {
    # Get dimension from input
    refOrder <- attributes(ref)$dimensions
  }
  
  att <- attributes(data)$dimensions
  if (is.null(att)) stop('No dimnames in the input data attributes, please use loadNcdf to load data.')
  if (identical(att, refOrder)) return(data)
  
  dimIndex <- seq(1, length(att))
  dimIndex1 <- grepAndMatch(refOrder, att)# match can apply to simple cases
  
  
  # for array this works, or setdiff can be used here to find the nomatch element.
  dimIndex2 <- dimIndex[-dimIndex1]# choose nomatch
  
  
  data <- aperm(data, c(dimIndex1, dimIndex2))
  attributes(data)$dimensions <- att[c(dimIndex1, dimIndex2)]
  
  return(data)
}

# Belong to checkDimLength
calcuDim <- function(data, dim) {
  dimIndex <- grepAndMatch(dim, attributes(data)$dimensions)
  dimLength <- dim(data)[dimIndex]
  return(dimLength)
}

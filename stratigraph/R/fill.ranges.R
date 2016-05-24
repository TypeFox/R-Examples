"fill.ranges" <-
function(x, out = 'internals', ...){

if(!is.strat.column(x)){
  warning('argument to x is not of class strat.column')
  x <- as.matrix(x)
  if(!is.numeric(x)){
    stop('argument to x cannot be coerced to a numeric matrix')
  }
  y <- x
}else{
  y <- x$counts
}

nas <- sum(is.na(y))
if(nas > 0){
  warning(paste(nas, 'NAs converted to zero counts'))
}

if(sum(y < 0) > 0) stop('negative values in count matrix')

internals <- function(z){
  n <- length(z)
  if(all(z == 0)){
    warning('taxon with no counts > 0')
    return(rep(NA, n))
  }
  if(all(z > 0)) return(rep(FALSE, n))
  z.range.index <- range((1:n)[z > 0], na.rm = TRUE)
  in.z.range <- c(rep(FALSE, z.range.index[1] - 1),
                  rep(TRUE, z.range.index[2] - z.range.index[1] + 1),
                  rep(FALSE, n - z.range.index[2]))
  z[z == 0 & in.z.range] <- -1
  return(z < 0)
}

fill.vector <- function(z){
  n <- length(z)
  min.count <- min(z[z > 0], na.rm = TRUE)
  if(all(z == 0)){
    warning('taxon with no counts > 0')
    return(z)
  }
  if(all(z > 0)) return(z)
  z.range.index <- range((1:n)[z > 0], na.rm = TRUE)
  in.z.range <- c(rep(FALSE, z.range.index[1] - 1),
                  rep(TRUE, z.range.index[2] - z.range.index[1] + 1),
                  rep(FALSE, n - z.range.index[2]))
  z[z == 0 & in.z.range] <- min.count
  return(z)
}



if(out == 'internals'){
  y <- apply(y, 2, internals)
}else if(out == 'min'){
  y <- apply(y, 2, fill.vector)
}else if(out == 'pa'){
  y <- apply(y, 2, fill.vector)
  y <- y > 0
}else{
  stop('argument to out not understood')
}

if(is.strat.column(x)){
  x$counts <- y
  return(x)
}else{
  return(y)
}

} # End of function
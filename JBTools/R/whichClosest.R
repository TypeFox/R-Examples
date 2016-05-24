whichClosest=function(
  ##title<< Find closest matches in two vectors
   x         ##<< numeric vector: values that should be found in x.match
   ,x.match  ##<< numeric vector: values the values of x should be matched against
  , arr.ind = FALSE ##logical: whether to return matrix indices or not.
)
  ##description<< Function to find closest matches of vector A in vector B and return
  ##              the respective indices.
{
  match.ind <- array(0, dim=c(length(x), if(arr.ind){2} else {1}))
  
  if (!arr.ind) {
    for (i in 1:length(x)) 
      match.ind[i,] <- min(which(abs(x[i] - x.match) == min(abs(x[i] - x.match), na.rm = TRUE)))
  } else {
    for (i in 1:length(x)) 
      match.ind[i,] <- which(abs(x[i] - x.match) == min(abs(x[i] - x.match), na.rm = TRUE), arr.ind = TRUE)[1,]   
  }

  ##value<< integer vector: indices of the closest element of x.match to each entry in x
  drop(match.ind)
}

names.SpatialStreamNetwork <-
function(x)
{
  d <- x@obspoints@SSNPoints[[1]]@point.data
  no <- names(d)
  np <- length(x@predpoints@SSNPoints)
  namesList <- vector("list",np + 1)
  namesList[[1]] <- no
  names4List <- "Obs"
  if(np > 0) {
    for(i in 1:np) {
        names4List <- c(names4List, x@predpoints@ID[[i]])
        d <- x@predpoints@SSNPoints[[i]]@point.data
        namesList[[i+1]] <- names(d)
      }
    }
  names(namesList) <- names4List
  return(namesList)
}


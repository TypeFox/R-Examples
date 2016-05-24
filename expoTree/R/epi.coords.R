epi.coords <- function(epi) {
  root <- epi$id[epi$parent==-1]
  get.coords <- function(id,top.coord) {
    off <- epi$parent == id
    coords <- cbind(id,top.coord)
    if (any(off)) {
      off.ids <- epi$id[off]
      off.times <- epi$itimes[off]
      for (i in length(off.ids):1) {
        off.coord <- get.coords(off.ids[i],top.coord+1)
        top.coord <- off.coord[[1]]
        coords <- rbind(coords,off.coord[[2]])
      }
    }
    return(list(max(coords[,2]),coords))
  }
  get.coords(root,0)
}


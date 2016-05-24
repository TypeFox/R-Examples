"create.grid" <- function(grid.list, sort=TRUE){

  d <- length(grid.list)
  if (sort){
    for (j in 1:d){
      grid.list[[j]] <- sort(grid.list[[j]])
    }
  }
  grid.list <- grid.list[d:1]
  grid <- as.matrix(expand.grid(grid.list))
  grid <- grid[,d:1,drop=FALSE]

  return(grid)
}

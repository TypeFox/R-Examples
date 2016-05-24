comp.div <-
function(dat, group, hier, q = 1, sims = TRUE) {
  hierarchies <- list(NA, length = length(hier))
  for (j in 1:hier) {
    new.dat <- matrix(nrow = length(unique(group[, 
      j])), ncol = dim(dat)[2])
    new.group <- array(dim = c(length(unique(group[, 
      j])), dim(group)[2] - j))
    for (i in 1:length(unique(group[, j]))) {
      target.rows <- which(group[, j] == unique(group[, 
        j])[i])
      if (length(target.rows) > 1) {
        new.dat[i, ] <- colMeans(dat[target.rows, 
          ])
      }
      else {
        X <- matrix(as.numeric(dat[target.rows, 
          ]), nrow = 1, ncol = dim(new.dat)[2])
        new.dat[i, ] <- colMeans(X)
      }
      new.group[i, ] <- group[target.rows[1], 
        c((1 + j):(dim(group)[2]))]
    }
    res <- div.part(new.dat, new.group, q = q)
    hierarchies[[j]] <- list(new.dat, new.group, 
      res)
  }
  new.dat <- dat
  new.group <- group
  res <- div.part(new.dat, new.group, q = q)
  hierarchies[[j + 1]] <- list(new.dat, new.group, 
    res)
  tempH <- list(NA, length = length(hierarchies))
  tempH[[1]] <- hierarchies[[j + 1]]
  for (i in 1:(length(hierarchies) - 1)) {
    tempH[[i + 1]] <- hierarchies[[i]]
  }
  hierarchies <- tempH
  if (sims == TRUE) 
    return(hierarchies)
  else if (sims == FALSE) {
    r.hierarchies <- list(NA, length = length(hierarchies))
    for (i in 1:length(hierarchies)) {
      r.hierarchies[[i]] <- hierarchies[[i]][[3]]
    }
  }
  return(r.hierarchies)
}

ExpandDataset <-
function(dataset, impt.times) {
  missing.t.list <- tapply(dataset[, 2], dataset[, 1],
                           function(X) impt.times[impt.times>max(X)])
  missing.t.vect <- as.vector(unlist(missing.t.list))
  n.missing.t <- as.vector(unlist(lapply(missing.t.list, length)))
  dataset.add <- data.frame(matrix(NA, nrow = length(missing.t.vect), 
                                   ncol = ncol(dataset)))
  names(dataset.add) <- names(dataset)
  dataset.add$obs.type <- rep(3, dim(dataset.add)[1])
  dataset.add[, 1] <- rep(unique(dataset[, 1]), n.missing.t)
  dataset.add[, 2] <- missing.t.vect 
  newmat2 <- rbind(dataset, dataset.add)
  newmat2 <- newmat2[order(newmat2[, 1], newmat2[, 2]), ]
}

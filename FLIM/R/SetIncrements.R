SetIncrements <-
function(dataset, impt.times, responses) {
  dataset <- dataset[dataset[, 2] %in% impt.times, ]
  IncrFun <- function(vec) c(diff(vec), NA)
  addme.col <- as.data.frame(apply(as.matrix(dataset[, responses]), 2, IncrFun))
  names(addme.col) <- paste0(responses, ".inc")
  dat.exp <- cbind(dataset, addme.col)
  return(dat.exp)
}

is.alignable <- function(a, threshold = .5){
  d <- dist.dna(a, model = "raw", pairwise.deletion = TRUE, as.matrix = TRUE)
  ub <- list(names(d[, 1])[which(d[, 1] > threshold)],
             names(d[, 1])[which(d[, 1] <= threshold)])
  ub <- ub[order(sapply(ub, length), decreasing = TRUE)]
  length(ub[[2]]) == 0
}
introNAs <- function(data, percent = 25) {
  n <- round((nrow(data) * ncol(data)) * (percent / 100))
  inds <- as.matrix(expand.grid(1:nrow(data), 1:ncol(data)))
  inds <- matrix(inds[!is.na(data[inds])], ncol = 2)
  selected <- inds[sample(nrow(inds), n), ]
  data[selected] <- NA
  data
}
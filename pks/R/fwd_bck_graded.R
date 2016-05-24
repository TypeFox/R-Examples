## Check for forward- and backward-gradedness

is.forward.graded <- function(K) {
  sapply(colnames(K), function(q) {
           K.plus <- K
           K.plus[, q] <- 1
           all(as.pattern(K.plus) %in% as.pattern(K))
         })
}


is.backward.graded <- function(K) {
  sapply(colnames(K), function(q) {
           K.minus <- K
           K.minus[, q] <- 0
           all(as.pattern(K.minus) %in% as.pattern(K))
         })
}


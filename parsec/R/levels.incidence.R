levels.incidence <-
function(x) {
  n <- nrow(x)
  lev <- rep(NA, n)
  names(lev) <- rownames(x)
  i <- 1
  sub <- (1:n)[is.na(lev)]
  while(length(sub)>0) {
    tmp_z <- x[sub, sub]
    class(tmp_z) <- class(x)
    quali <- which(depths(tmp_z)==1)
    lev[sub[quali]] <- i
    i <- i + 1
    sub <- (1:n)[is.na(lev)]
  }
  lev
}

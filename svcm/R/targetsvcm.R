targetsvcm <- function(lambda) {
  
  predictsvcm(lambda, type = type)
  if (type == "TP") {
    GCV <- drop(as.matrix(n * crossprod(y - eta) / (n - ED)^2))
  } else if (type == "SEQ") {
    GCV <- n * sum((Y - eta)^2) / (n - ED)^2
  }
  assign("GCVtab", rbind(GCVtab, c(lambda, GCV)), env = svcm.env)
  return(GCV)
  
}

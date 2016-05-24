comp.test <- function(cl, spg){
  cs <- chisq.test(cl,spg,simulate.p.value=TRUE,B=10000)
  print(cs)
  invisible(cs)
}

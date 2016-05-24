delayt <- function(sec){
  start <- Sys.time()
  kk <- 0
  while (Sys.time() - start <= sec) {kk <- kk+1}
  kk
}

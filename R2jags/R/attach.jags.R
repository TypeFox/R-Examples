

attach.jags <- function(x, overwrite = NA){
  r2 <- attach.bugs(x$BUGSoutput, overwrite = overwrite)
  invisible(r2)
}

detach.jags <- function(){
  detach.bugs()
}

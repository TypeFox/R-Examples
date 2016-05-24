# Code is based on spdep::do_ldet written by Roger Bivand

lndet_imrw <- function(W){
  detval <- NULL
  rmin <- -0.9
  rmax <- 0.99
  env <- new.env(parent=globalenv())
  lw <- mat2listw(W)
  
  assign("n", nrow(W), envir=env)
  assign("listw", lw, envir=env)
  assign("similar", FALSE, envir=env)
  assign("family", "SAR", envir=env)
  
  mcdet_setup(env)
  detval1 <- seq(rmin, rmax, 0.001)
  detval2 <- sapply(detval1, do_ldet, env)
  detval  <- cbind(detval1, detval2)
  
  return( detval )
}

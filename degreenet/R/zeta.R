zeta <- function(v){
  # compute Reimann zeta function
  if (v[1] <= 1)
	stop("Invalid argument. Riemann Zeta function is only defined for v > 1.");
  return(.C("zetaC",v=as.double(v[1]),z=double(1), PACKAGE="degreenet")$z)
}

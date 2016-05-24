`rnchg` <- function(n, n1, n2, m1, psi, debug=0)
{
  .Call("rnchg", PACKAGE="RxCEcolInf",
	list("n"=as.integer(n),
	     "n1"=as.double(n1),
	     "n2"=as.double(n2),
	     "m1"=as.double(m1),
	     "psi"=as.double(psi)))
}

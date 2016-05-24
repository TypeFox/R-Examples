"gridCOP" <- function(cop=NULL, para=NULL, delta=0.05, transpose=TRUE, ...) {
   u <- seq(0,1, by=delta); v <- seq(0,1, by=delta)
   h <- sapply(u, function(au) {
        sapply(v, function(av) { return(cop(au,av, para=para, ...)) } ) })
   # The transpose is critical so that image() function will work correctly.
   if(transpose) h <- t(h)
   return(h)
}

real2canon <- function(R,B=c(0,1))
{
   storage.mode(R) <- "double"
   storage.mode(B) <- "integer"
   .Call("RealToCanonicalForR", R, B)
}

computeMLE <- function(R,B=c(0,1),max.inner=10,max.outer=1000,tol=1e-10)
{
   storage.mode(R) <- "double"
   storage.mode(B) <- "integer"
   storage.mode(max.inner) <- "integer"
   storage.mode(max.outer) <- "integer"
   storage.mode(tol) <- "double"
   .Call("ComputeMLEForR", R, B, max.inner, max.outer, tol)
}

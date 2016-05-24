canon2real <- function(Rcanon, R, B=c(0,1))
{
   storage.mode(Rcanon) <- "integer"
   storage.mode(R) <- "double"
   storage.mode(B) <- "integer"
   .Call("CanonicalToRealForR", Rcanon, R, B)
}

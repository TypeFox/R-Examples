reduc <- function(R,B=c(0,1),hm=FALSE,cm=FALSE)
{
   storage.mode(R) <- "double"
   storage.mode(B) <- "integer"
   storage.mode(hm) <- "integer"
   storage.mode(cm) <- "integer"
   .Call("ReductionStepForR", R, B, hm, cm)
}

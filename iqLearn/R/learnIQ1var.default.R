learnIQ1var.default <-
function (object, H1CVar=NULL, s1sInts=NULL, method="homo", ...){ 
  est <- iqQ1varEst (object, H1CVar, s1sInts, method, ...); 

  class (est) <- "learnIQ1var";
  est
}

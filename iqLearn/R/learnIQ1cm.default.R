learnIQ1cm.default <-
function (object, H1CMean, A1, s1cmInts=NULL,
  ...){  
  est <- iqQ1cmEst (object, H1CMean, A1, s1cmInts, ...); 

  class (est) <- "learnIQ1cm";
  est
}

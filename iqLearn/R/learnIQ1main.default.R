learnIQ1main.default <-
function (object, H1Main, A1, s1mainInts=NULL,
  ...){  
  est <- iqQ1MainEst (object, H1Main, A1, s1mainInts, ...); 

  class (est) <- "learnIQ1main";
  est
}

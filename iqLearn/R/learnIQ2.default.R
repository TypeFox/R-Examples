learnIQ2.default <-
function (H2, Y, A2, s2ints, ...){
  est <- iqQ2Est (H2, Y, A2, s2ints, ...); 

  class (est) <- "learnIQ2";
  est
}

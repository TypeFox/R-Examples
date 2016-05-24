qLearnS2.default <-
function (H2, Y, A2, s2ints=NULL, ...){ 
  est <- qLearnS2Est (H2, Y, A2, s2ints, ...);

  class (est) <- "qLearnS2";
  est
}

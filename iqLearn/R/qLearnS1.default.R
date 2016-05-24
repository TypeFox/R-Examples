qLearnS1.default <-
function (object, H1q, A1, s1ints=NULL, ...){ 
  est <- qLearnS1Est (object, H1q, A1, s1ints, ...);

  class (est) <- "qLearnS1";
  est
}

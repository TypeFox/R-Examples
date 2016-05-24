value.default <-
function (d1, d2, Y, A1, A2){
  est <- valueEst (d1, d2, Y, A1, A2);
  class (est) <- "value";

  est
}

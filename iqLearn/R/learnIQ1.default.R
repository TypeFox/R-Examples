learnIQ1.default <-
function (mainObj, cmObj, sigObj, dens="nonpar"){
  est <- learnIQ1Est (mainObj, cmObj, sigObj, dens);
  class (est) <- "learnIQ1";
  est
}

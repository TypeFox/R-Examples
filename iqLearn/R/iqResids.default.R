iqResids.default <-
function (object){
  resids = object$stdResids;
  class (resids) <- "iqResids";
  resids;
}

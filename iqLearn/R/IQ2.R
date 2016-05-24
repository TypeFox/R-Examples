IQ2 <-
function (object, h2){
  h2pos = c (1, h2, 1, h2[object$s2ints]);
  q2Pos = sum (h2pos*c (object$betaHat20, object$betaHat21));
  h2neg = c (1, h2, -1, -h2[object$s2ints]);
  q2Neg = sum (h2neg*c (object$betaHat20, object$betaHat21));
  q2opt = sign (q2Pos - q2Neg);
  return (list ("q2Pos"=q2Pos, "q2Neg"=q2Neg, "q2opt"=q2opt));
}

qLearnQ1 <-
function (object, h1q){

  if (length (object$s1ints) > 0){
    h1pos = c (1, h1q, 1, h1q[object$s1ints]);
    q1Pos = sum (h1pos*c (object$betaHat10, object$betaHat11));
    h1neg = c (1, h1q, -1, -h1q[object$s1ints]);
    q1Neg = sum (h1neg*c (object$betaHat10, object$betaHat11));
  }
  else{
    h1pos = c (1, h1q, 1);
    q1Pos = sum (h1pos*c (object$betaHat10, object$betaHat11));
    h1neg = c (1, h1q, -1);
    q1Neg = sum (h1neg*c (object$betaHat10, object$betaHat11));
  }
  
  q1opt = sign (q1Pos - q1Neg);
  return (list ("q1Pos"=q1Pos, "q1Neg"=q1Neg, "q1opt"=q1opt));
}

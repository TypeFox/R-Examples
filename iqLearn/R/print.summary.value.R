print.summary.value <-
function (x, ...){

  cat ("Estimated Optimal Value: \n");
  print (x$optVal)
  cat ("Estimated Value of (1, 1) Regime: \n");
  print (x$PosPos)
  cat ("Estimated Value of (1, -1) Regime: \n");
  print (x$PosNeg)
  cat ("Estimated Value of (-1, 1) Regime: \n");
  print (x$NegPos)
  cat ("Estimated Value of (-1, -1) Regime: \n");
  print (x$NegNeg)
}

print.summary.learnIQ1var <-
function (x, ...){

  if (x$varType){
    cat ("Constant standard deviation: \n");
    print (x$s1Reg)
  }
  else{
    cat ("Variance Model: \n");
    print (summary (x$s1Reg));
  }
}

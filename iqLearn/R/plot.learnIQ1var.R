plot.learnIQ1var <-
function (x, ...){
  if (!x$homo){
    par (mfrow=c(2,2))
    plot (x$s1VarFit, ...);
  }
  else{
    print (paste ("Variance assumed constant"));
  }
}

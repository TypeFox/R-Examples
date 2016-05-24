print.BAEssd <-
function(x,...){
  cat("\nBayesian Average Error Sample Size Determination\n")
  
  cat("\nCall: ")
  print(x$call)
  
  cat("\nSample Size: ",x$n)
  cat("\nTotal Average Error: ",attr(x$n,"TE"))
  
  if(attr(x$n,"TE")<=attr(x$n,"alpha")){
    cat("\n\nAcceptable sample size determined!\n\n")
  }
  else{
    cat("\n\nMaximum sample size attempted and was unacceptable. ",
        "Try a larger sample size.\n\n")
  }
}

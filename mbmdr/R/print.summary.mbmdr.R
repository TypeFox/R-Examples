print.summary.mbmdr <-
function(x,...){
  x <- as.data.frame(x) 
  if(x[1,"MIN.P"]>0.05){
  	cat("\n No significative interaction models found\n")
  	return(NULL)
  }
  else cat("\n Significative interaction models found are:\n\n")
  print(x)
}


computeRectangles<- function(series,alpha){
BPDriver <- .jnew("edu.dtulnu.stat.brpts.rect.BreakPointsMinMaxDriver") # create instance of BreakPointsMinMaxDriver class
out <- .jcall(BPDriver,  "[[D" , "computeRectangles",series,alpha,evalArray=TRUE) # invoke computeRectangles method
out <- sapply(out, .jevalArray)
out <- t(out)
return(out)
}

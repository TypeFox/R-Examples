computeBreakPoints<- function(series,alpha){
BPDriver <- .jnew("edu.dtulnu.stat.brpts.rect.BreakPointsMinMaxDriver") # create instance of BreakPointsMinMaxDriver class
out <- .jcall(BPDriver,  "[I" , "computeBreakPoints", series, alpha) # invoke computeBreakPoints method
return(out)
}

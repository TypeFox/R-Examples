print.bnp<-function(x,digits=max(3,getOption('digits')-3), ...){
   #Print formula
   cat('\nCall: ',deparse(x$call),'\n\n')
}

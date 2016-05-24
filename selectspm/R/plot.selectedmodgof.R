plot.selectedmodgof<- function(x,...){
	main=deparse(substitute(x))
	plot(x$envelopes,main=main,...)
}
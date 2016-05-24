"plot.crit.fun" <-
function(
	x,#an "opt.par" class object
	main=NULL,
	...	#aditional parameters to "plot.mat"
){
	if(is.null(main)) main <- deparse(substitute(x))
	plot.mat(x$M,clu=x$clu,IM=x$IM,main=main,...)
}


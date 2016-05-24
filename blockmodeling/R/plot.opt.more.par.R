"plot.opt.more.par" <-
function(
	x,#an "opt.par.mode" class object
	main=NULL,
	which=1,	#which (if there are more than one) of optimal solutions to plot
	...	#aditional parameters to "plot.mat"
){
	if(is.null(main)) main <- deparse(substitute(x))
	if(which>length(x$best)){
		warning("The selected (",which,") best solution does not exist!\nOnly ", length(x$best)," best solution(s) exist(s).\nThe first best solution will be ploted.\n")
		which<-1
	}
	plot.mat(x$M,clu=x$best[[which]]$clu,IM=x$best[[which]]$IM,main=main,...)
}



plot.Daim <- function(x, method=NULL, all.roc=FALSE, color="red", alpha=0.25,
					type="s", xlab="False positive rate", ylab="True positive rate", 
					main=NULL, add=FALSE, legend = FALSE, ...)
{		
	if(class(x)[2] != "cv"){
		if(is.null(method))
			method <- "0.632+"
		meth <- charmatch(method,c("0.632+","0.632","loob","sample"),nomatch = 0)
		if(meth == 0)
			stop("\nThe value of 'method' must be one of '0.632+', '0.632', 'loob' or 'sample'\n")
		if(is.null(main) && !add)
			main <- paste("Method -",method)
		if(meth == 4)
			all.roc <- TRUE
		if(!add){
			plot(-1, -1, xlim=c(0,1), ylim=c(0,1), xlab=xlab, ylab=ylab, main=main)
			grid()
		}
		if(all.roc || meth == 4){
			for(i in 1:length(x$sample.roc)){
				lines(x$sample.roc[[i]], col=rgb(0, 1, 0, 255*alpha, maxColorValue=255), ...)
			}
		}
		if(meth == 1){
			lines(1 - c(x$roc$spec632p, 1), c(x$roc$sens632p, 0), col=color, type=type, ...)
			if(legend){
				lines(c(0,1), c(0,1), col="black")
				points(1-x$spec632p, x$sens632p, col=1, pch=19)
				legend("bottomright", legend=
					   c(paste("cut-off:", formatC(x$method$cutoff, digits=max(3, getOption("digits") - 3))),
						 paste("FPR:", formatC(1-x$spec632p, digits=max(3, getOption("digits") - 3))),
						 paste("TPR:", formatC(x$sens632p, digits=max(3, getOption("digits") - 3))),
						 paste("AUC:", formatC(auc(x$roc$spec632p, x$roc$sens632p), 
											   digits=max(3, getOption("digits") - 3)))), inset=0.01)
			}
		}else if(meth == 2){
			lines(1 - c(x$roc$spec632, 1), c(x$roc$sens632, 0), col=color, type=type, ...)
			if(legend){
				lines(c(0,1), c(0,1), col="black")
				points(1-x$spec632, x$sens632, col=1, pch=19)
				legend("bottomright", legend=
					   c(paste("cut-off:", formatC(x$method$cutoff, digits=max(3, getOption("digits") - 3))),
						 paste("FPR:", formatC(1-x$spec632, digits=max(3, getOption("digits") - 3))),
						 paste("TPR:", formatC(x$sens632, digits=max(3, getOption("digits") - 3))),
						 paste("AUC:", formatC(auc(x$roc$spec632, x$roc$sens632), 
											   digits=max(3, getOption("digits") - 3)))), inset=0.01)
			}
		}else if(meth == 3){
			lines(1 - c(x$roc$specloob,1), c(x$roc$sensloob,0), col=color, type=type, ...)
			if(legend){
				lines(c(0,1), c(0,1), col="black")
				points(1-x$specloob, x$sensloob, col=1, pch=19)
				legend("bottomright", legend=
					   c(paste("cut-off:", formatC(x$method$cutoff, digits=max(3, getOption("digits") - 3))),
						 paste("FPR:", formatC(1-x$specloob, digits=max(3, getOption("digits") - 3))),
						 paste("TPR:", formatC(x$sensloob, digits=max(3, getOption("digits") - 3))),
						 paste("AUC:", formatC(auc(x$roc$specloob, x$roc$sensloob), 
											   digits=max(3, getOption("digits") - 3)))), inset=0.01)
			}
		}
	}
	else{
		if(is.null(method))
			method <- "cv"
		meth <- charmatch(method, c("cv", "sample"), nomatch = 0)
		if(meth == 0)
			stop("\nThe value of 'method' must be one of 'cv' or 'sample'\n")
		if(is.null(main) && !add)
			main <- paste("Method -",method)
		if(meth == 2)
			all.roc <- TRUE
		if(!add){
			plot(-1,-1, xlim=c(0,1), ylim=c(0,1), xlab=xlab, ylab=ylab,
				 main=main)
			grid()
		}
		if(all.roc || meth == 2){
			for(i in 1:length(x$sample.roc)){
				lines(x$sample.roc[[i]], col=rgb(0, 1, 0, 255*alpha, maxColorValue=255), ...)
			}
		}
		if(meth == 1)
			lines(1 - c(x$roc$specloob,1), c(x$roc$sensloob,0), col=color, type=type, ...)
		lines(c(0,1), c(0,1), col="black")
		points(1-x$specloob, x$sensloob, col=1, pch=19)
		legend("bottomright", legend=
			   c(paste("cut-off:", formatC(x$method$cutoff, digits=max(3, getOption("digits") - 3))),
				 paste("FPR:", formatC(1-x$specloob, digits=max(3, getOption("digits") - 3))),
				 paste("TPR:", formatC(x$sensloob, digits=max(3, getOption("digits") - 3))),
				 paste("AUC:", formatC(auc(x$roc$specloob, x$roc$sensloob), digits=max(3, getOption("digits") - 3)))), inset=0.01,...)
	}
}





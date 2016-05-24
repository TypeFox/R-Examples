`plot.raggr` <-
function(x, show.average = TRUE, show.legend = TRUE, colR="red", ...){	
	k <- length(x$top.list)

	if(is.null(x$num.iter))
		mat = matrix(c(1,2,1,2), 2)
	else
		mat = matrix(c(1, 3, 2, 3), 2)

	layout(mat)
	
	if(!is.null(x$num.iter)){
		plot(1:x$num.iter, x$summary[,1], col="red", type="l", xlab="Iteration", ylab="Scores",
			main="Minimum Path", bty="n")
		legend("topright", legend=c(paste("min =", round(x$optimal.value,3))), box.lty=0)
	}
	
	if(all(x$sample == min(x$sample)))
		hist(c(x$sample, x$sample[1]+.0001), x$sample.size, col="blue", border="darkblue", 
			xlab="Objective function scores", main=paste("Final Sample Distribution"),
			xlim=c(x$sample[1]-.5, x$sample[1]+.5))
	else
		hist(x$sample, x$sample.size, col="blue", border="darkblue", xlab="Objective function scores",
			main=paste("Final Sample Distribution"))

	xlabel <- paste("Optimal List:", paste(x$top.list, collapse=" "))
	if(nchar(xlabel)>75)
		xlabel="Optimal List"

	for(i in 1:nrow(x$lists)){
		mr <- match(x$top.list, x$lists[i,])
		mr[is.na(mr)] <- k+1
		if(i==1){
			plot(1:k, mr+runif(1,0,.1), ylim=c(0,k+2),
				xlab=xlabel, ylab="Ranks", main="Rank Aggregation", type="l", col="grey", xaxt="n", ...)
			axis(1,at=1:k, labels=x$top.list)}
		else
			lines(1:k, mr+runif(1,0,.1), col="grey", ...)
	}
	lines(1:k, 1:k, col=colR, ...)
	if(show.average){
		lines(1:k, apply(apply(x$lists, 1, function(z) {mr <- match(x$top.list, z)
										mr[is.na(mr)] <- k+1
										mr}),1,mean), col="black", ...) # 1 in the outside apply because inside apply
															  # transposes the result
		if(show.legend)
			legend("topleft", legend=c("Data", x$method, "Mean"), fill=c("grey",colR,"black"), horiz=TRUE)
	} else
		if(show.legend)
			legend("topleft", legend=c("Data", x$method), fill=c("grey",colR), horiz=TRUE)
   }
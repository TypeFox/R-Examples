plotfplsr <-
function(x, xlab1 = x$ypred$xname, ylab1 = "Basis function", xlab2 = "Time", ylab2 = "Coefficient", 
  mean.lab = "Mean", interaction.title = "Interaction")
{
	scores = x$T
	order  = dim(scores)[2]
	pred   = x$P
	resp   = x$Q
	par(mfrow = c(3, (order + 1)))
	plot(x$y1, x$meanX$y, type = "l", xlab = xlab1, ylab = mean.lab, main = "Predictor")
	for(i in 1:order)
	{
		plot(x$y1, pred[,i], type = "l", xlab = xlab1, ylab = paste(ylab1, i, sep = " "))
	}
	plot(x$y1, x$meanY$y, type = "l", xlab = xlab1, ylab = mean.lab, main = "Response")
	for(i in 1:order)
	{
		plot(x$y1, resp[,i], type = "l", xlab = xlab1, ylab = paste(ylab1, i, sep = " "))
	}
	plot(x$y1, resp[,1], type = "n", xlab = "", ylab = "", xaxt = "n", yaxt = "n", bty = "n")
	for(i in 1:order)
	{
		plot(x$x1, scores[,i], type = "l", xlab = xlab2, ylab = paste(ylab2, i, sep = " "))
	}
}



## Plotting in 1 or 2 dimensions
##
## May 8th, 2012

plot.GP <- function(x,M=1,range=c(0,1),resolution=50,colors=c('black','blue','red'),
	line_type=c(1,2),pch=20,cex=1,legends=FALSE,surf_check=FALSE,response=TRUE,...){
if (is.GP(x) == FALSE){
	stop("The object in question is not of class \"GP\" \n")
}
X = x$X;
Y = x$Y;
n = nrow(X);
d = ncol(X);

if (d>= 3){ stop("can not plot in higher than 2 dimensions.\n")}

if (d==1){
	# make the prediction values
	xvec = matrix(seq(from=range[1],to=range[2],length.out=resolution),ncol=1)
	GPprediction=predict.GP(x,xvec,M)
	Y_hat = GPprediction$Y_hat;
	MSE = GPprediction$MSE;

	# Finding a good range to plot over
	max_height = max(Y_hat, Y, Y_hat+2*sqrt(MSE),Y_hat-2*sqrt(MSE))
	min_height = min(Y_hat, Y, Y_hat+2*sqrt(MSE),Y_hat-2*sqrt(MSE))
	max_length = max(X,xvec)
	min_length = min(X,xvec)

	# Plotting the Prediction and errors
	leg.txt = c(expression("Model Prediction:     " ~ hat(y)(x)), 
		expression("Uncertanity Bounds: " ~ hat(y)(x) %+-% 2 %*% s(x)),
		"Design Points")
	matplot(X,Y,cex=cex[1],col=colors[1],pch=pch[1],
   		ylim=c(min_height,max_height),xlim=c(min_length,max_length),
		xlab = "x (Input Variable)", ylab = "Model Prediction")
	lines(xvec,Y_hat,col=colors[2],lty=line_type[1])			#Predicted Function
	lines(xvec,Y_hat-2*sqrt(MSE),col=colors[3],lty=line_type[2])	#Errors
	lines(xvec,Y_hat+2*sqrt(MSE),col=colors[3],lty=line_type[2])	#Errors
	if(legends==TRUE){
		legend(min_length, max_height, leg.txt, col = c(colors[2], colors[3], colors[1]), 
			lty = c(line_type[1],line_type[2], -1), pch = c(-1, -1, pch[1]), 
			pt.cex = cex[1])
	}
}

if (d==2){
	# Making the predictions
	xvector = seq(from=range[1],to=range[2],length.out=resolution)
	xvec = expand.grid(x = xvector, y=xvector)
	xvec = as.matrix(xvec)
	GPprediction = predict.GP(x,xvec,M)
	Y_hat = GPprediction$Y_hat;
	MSE = GPprediction$MSE;
	dim(Y_hat) = c(length(xvector),length(xvector))
	dim(MSE) = c(length(xvector),length(xvector))

	if (surf_check == TRUE){
		## Wireframe Plots
		if (response == TRUE){
			h1 = wireframe(Y_hat,scales = list(arrows=FALSE),row.values = xvector,
				column.values = xvector,
				xlab = expression(X[1]),ylab = expression(X[2]),
				zlab = list('Model Prediction', rot = 90),...)
			print(h1)
		} else {
			h2 = wireframe(MSE,scales = list(arrows=FALSE),row.values = xvector,
				column.values = xvector,
				xlab = expression(X[1]),ylab = expression(X[2]),
				zlab = list('MSE', rot = 90),...)
			print(h2)
		}

	} else {
		## Contour Plotting
		if (response == TRUE){
			h1 = levelplot(Y_hat, xlab = expression(X[1]),ylab = expression(X[2]), 
				row.values = xvector, column.values = xvector,
				xlim = range, ylim = range, ...)
			print(h1)
		} else {
			h2 = levelplot(MSE, xlab = expression(X[1]),ylab = expression(X[2]), 
				row.values = xvector, column.values = xvector, 
				xlim = range, ylim = range, ...)
			print(h2)
		}
	}
}
}
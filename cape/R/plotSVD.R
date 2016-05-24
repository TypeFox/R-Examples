plotSVD <-
function(data.obj, orientation = c("vertical", "horizontal"), neg.col = "blue", pos.col =  "brown", light.dark = "dark"){
	
	
	#test to see if there is a v in the orientation vector
	orient.test <- grep("v", orientation) 

	if(is.null(data.obj$right.singular.vectors)){
		stop("get.eigentraits() must be run before plotSVD.")
		}

	#The contribution of each eigentrait to each phenotype
	#is encoded in the right singular vectors (eigenvectors)	
	eigen.weights <- t(data.obj$right.singular.vectors)[1:dim(data.obj$ET)[2],]
	colnames(eigen.weights) <- colnames(data.obj$pheno)
	rownames(eigen.weights) <- colnames(data.obj$ET)
	
	#use the singular values to calculate the variance
	#captured by each mode
	singular.vals <- data.obj$singular.values[1:dim(data.obj$ET)[2]]
	
	if(is.null(singular.vals)){
		stop("get.eigentraits() must be run before plotSVD()")
		}
		
	eigen.vals <- singular.vals^2
	var.accounted <- eigen.vals/sum(eigen.vals)


	#if the orientation is vertical, use the matrix as is.
	#otherwise, rotate the matrix, so the ETs are plotted 
	#in each row
	if(length(orient.test) > 0){
		rotated.eigen <- eigen.weights
		}else{
			rotated.eigen <- rotate.mat(eigen.weights)
			}


	min.weight <- min(eigen.weights)
	max.weight <- max(eigen.weights)
	
	pos.cols <- get.col(pos.col, light.dark)
	neg.cols <- get.col(neg.col, light.dark)[3:1]

	mypal.pos <- colorRampPalette(pos.cols)
	mypal.neg <- colorRampPalette(neg.cols)


	ColorLevels <- seq(min.weight, max.weight, length=256)
	ColorRamp <- c(mypal.neg(length(which(ColorLevels < 0))), mypal.pos(length(which(ColorLevels >= 0))))
	

	barplot.axis.cex <- 1.7; barplot.labels.cex = 2; barplot.title.cex = 1.7

	#plot the horizontal configuration
	if(length(orient.test) == 0){
		layout.mat <- matrix(c(0,1,1,0,3,2,5,7,0,4,6,0), nrow = 3, ncol = 4, byrow = TRUE)
		layout(layout.mat, heights = c(0.15, 0.7, 0.15), widths = c(0.1, 0.35, 0.45, 0.1))


		#1) plot the title in its own box
		par(mar = c(0,0,0,0))
		plot.new()
		plot.window(xlim = c(0, 1), ylim = c(0, 1))
		text(0.5, 0.5, expression(bold("Eigentrait Contributions to Phenotypes")), cex = 2)

		#2) plot the weight matrix
		par(mar = c(0, 0, 2, 0))
		# myImagePlot(rotated.eigen)
		image(x = 1:length(rotated.eigen[,1]), y = 1:length(rotated.eigen[1,]), rotated.eigen, col = ColorRamp, axes = FALSE, xlab = "", ylab = "")
		
		#3) fit in the text for the y axis of the weight matrix
		par(mar = c(0,0,0,0))
		plot.new()
		plot.window(xlim = c(0, 1), ylim = c(0, length(eigen.weights[,1])))
		text(x = 0.5, y = seq(0.5, (length(eigen.weights[,1])-0.5), 1), rev(rownames(eigen.weights)), cex = 1.7) #the labels for the y axis (phenotypes)
		
		#4) fit in the text for the x axis of the weight matrix
		par(mar = c(0,0,0,0))
		plot.new()
		plot.window(ylim = c(0, 1), xlim = c(0, length(eigen.weights[1,])))
		text(y = 0.5, x = seq(0.5, (length(eigen.weights[1,])-0.5), 1), colnames(eigen.weights), cex = 1.7, srt = 90) #the labels for the y axis (phenotypes)

		
		#5) plot the barplot
		par(mar = c(0, 0, 2, 2))
		barplot(rev(var.accounted), names = "", horiz = TRUE, cex.axis = barplot.axis.cex, cex.lab = barplot.labels.cex, xlab = "")
		
		#6) add the axis label for the barplot
		par(mar = c(0,0,0,0))
		plot.new()
		plot.window(xlim = c(0, 1), ylim = c(0, 1))
		text(x = 0.4, y = 0.5, labels = "Percent Total Variance", cex = barplot.title.cex) #the labels for the y axis (phenotypes)
		
		
		#7) add the color ramp
		par(mar = c(0,2,2,2))
		image(1, ColorLevels, matrix(data=ColorLevels, ncol=length(ColorLevels),nrow=1), col=ColorRamp, xlab="",ylab="",xaxt="n")		

		# dev.off()
		}

	#If the orientation is set to vertical, there are lots of little parameters to adjust...
	if(length(orient.test) > 0){
		
		layout.mat <- matrix(c(0,1,0,3,2,0,5,4,7,0,6,0), nrow = 4, ncol = 3, byrow = TRUE)
		layout(layout.mat, widths = c(0.2, 0.7, 0.1), heights = c(0.1, 0.35, 0.45, 0.1))

		#1) plot the title in its own box
		par(mar = c(0,0,0,0))
		plot.new()
		plot.window(xlim = c(0, 1), ylim = c(0, 1))
		text(0.5, 0.5, expression(bold("Eigentrait Contributions to Phenotypes")), cex = 2)
		
		#2) plot the barplot
		par(mar = c(0, 0, 4, 2))
		mids <- barplot(var.accounted, plot = FALSE)
		barplot(var.accounted, names = "", horiz = FALSE, cex.axis = barplot.axis.cex, cex.lab = barplot.labels.cex, ylab = "Percent Total Variance")
		
		#3) add the axis label for the barplot
		par(mar = c(0,0,0,0))
		plot.new()
		plot.window(xlim = c(0, 1), ylim = c(0, 1))
		text(x = 0.65, y = 0.4, labels = "Percent Total Variance", cex = barplot.title.cex, srt = 90) #the labels for the y axis (phenotypes)
		
		
		#3) plot the weight matrix
		par(mar = c(0, 0, 0, 2))
		image(x = 1:length(eigen.weights[,1]), y = 1:length(eigen.weights[1,]), rotated.eigen, col = ColorRamp, axes = FALSE, xlab = "", ylab = "")

		#4) fit in the text for the y axis of the weight matrix
		par(mar = c(0,0,0,0))
		plot.new()
		plot.window(xlim = c(0, 1), ylim = c(0, length(eigen.weights[1,])))
		text(x = 0.5, y = seq(0.5, (length(eigen.weights[1,])-0.5), 1), colnames(eigen.weights), cex = 1.7) #the labels for the y axis (phenotypes)

		#5) fit in the text for the x axis of the weight matrix
		par(mar = c(0,0,0,0))
		plot.new()
		plot.window(xlim = c(0, length(eigen.weights[,1])), ylim = c(0, 1))
		text(x = seq(0.5, (length(eigen.weights[,1])-0.5), 1), y = 0.5, rownames(eigen.weights), cex = 2) #the labels for the x axis (ETs)
	
		#6) add the color ramp
		par(mar = c(0,2,0,2))
		image(1, ColorLevels, matrix(data=ColorLevels, ncol=length(ColorLevels),nrow=1), col=ColorRamp, xlab="",ylab="",xaxt="n", cex.axis = 2)		

		names(var.accounted) <- colnames(data.obj$ET)
		invisible(var.accounted)

		# dev.off()			
		
		}


}

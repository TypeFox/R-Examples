result.pic.scatterplot <- function(values, picPath=getwd(), 
					picName="beech-budburst", 
					createFile=TRUE){
	x <- values$doy.observed
	y <- values$doy.model

	check <- which((x > -9999) & (!is.na(x)) & (y > -9999) & (!is.na(y)))

	x <- x[check]
	y <- y[check]

	min.val <- min(unique(x,y), na.rm=TRUE)
	max.val <- max(unique(x,y), na.rm=TRUE)

	lin.model <- lm(y~x)
	slm <- summary(lin.model)
	test <- try(pvalue <- as.numeric(pf(slm$fstatistic[1], slm$fstatistic[2], 
			slm$fstatistic[3], lower.tail=FALSE)),silent=TRUE)
	if (inherits(test, "try-error")){
		subtitle.text <- "Linear model calculation failed"
		lm.success <- FALSE
	} else {
		subtitle.text <- paste("Linear model: ",
				round(as.numeric(lin.model$coefficients[2]),3),"x", 
				ifelse(as.numeric(lin.model$coefficients[1]) >= 0, " + ", " - "), 
				round(abs(as.numeric(lin.model$coefficients[1])),3), 
				"; p-value of F-test: ",
				ifelse(pvalue==0, " < 1e-15",pvalue),sep="")
		lm.success <- TRUE
	}

	if (createFile){ png(filename=paste(picPath,"/",picName,"-scatterplot.png",sep=""), 
				width=1000, height=1000, res=100) }
	plot(x=x, y=y, xlab="Observed", ylab="Modelled",
		xlim=c(min.val,max.val), ylim=c(min.val,max.val), font.axis=2, 
		font.lab=2, lwd=2, main="Scatterplot", 
		font.sub=2,cex.sub=1.3, col.sub="blue", sub=subtitle.text)
	axis(2, labels=FALSE, at=seq(from=0, to=360, by=10), tck=0.02)
	axis(2, labels=FALSE, at=seq(from=0, to=366, by=1), tck=0.01)
	axis(1, labels=FALSE, at=seq(from=0, to=360, by=10), tck=0.02)
	axis(1, labels=FALSE, at=seq(from=0, to=366, by=1), tck=0.01)
	
	if (lm.success) {
		lines(min.val:max.val, min.val:max.val, col="green", lwd=2)
		abline(lin.model, lwd=2, col="orange")
		legend(x="topleft", bty="n", legend=c("Linear model fit", 
			"Line indicating ideal fit"), col=c("orange", "green"), 
			pch=c(-1,-1), lty=c(1,1), lwd=c(2,2))
	}
	if (createFile){ dev.off() }

	return(TRUE)
}

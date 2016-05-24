
ROC.graphic.ktsp <- function(roc, m = 1, boxplot = TRUE, AUC = TRUE, auc.x = 0.55, auc.y = 0.25, box.col = "yellow", line.col = "red", multiple.col = c("red", "orange", "blue"), maintitle = NULL, mtext = NULL, undertitle = NULL, graphic = TRUE){

	mult.cutoff <- roc$mult.cutoff
	spec_list <- roc$spec
	sens_list <- roc$sens
	n <- roc$n
	
	if(is.null(mult.cutoff)){mult.cutoff <- FALSE}
	if(mult.cutoff){
		spec_25 <- list(length=n)
		sens_25 <- list(length=n)

		for(i in 1:n){
			spec_25[[i]] <- unique(spec_list[[1]][i,])
			sens_25[[i]] <- c()
			mean <- c()
			for(j in 1:length(spec_25[[i]])){
				mean[j] <- mean(sens_list[[1]][i,][spec_list[[1]][i,]==spec_25[[i]][j]], na.rm=TRUE)
			}
			sens_25[[i]] <- mean
		}

		spec_25_1 <- 1-unlist(spec_25)

		spec_25_1_order <- spec_25_1[order(spec_25_1)]
		sens_25_order <- unlist(sens_25)[order(spec_25_1)]

		uniquespec_25_1 <- unique(spec_25_1_order)

		box <- list()
		for(i in 1:length(uniquespec_25_1)){
			box[[i]]<- sens_25_order[which(spec_25_1_order==uniquespec_25_1[i])]
		}

		spec_1 <-c()
		sens <- c()
		median1 <- c()
		at1 <- c()

		for(i in 1:(length(box)/m)){
			sens <- c(sens, box[[m*i]])
			spec_1 <- c(spec_1, rep(uniquespec_25_1[m*i],length(box[[m*i]])))
			median1[i] <- median(box[[m*i]])
			at1[i] <- uniquespec_25_1[m*i]
		}

		if(graphic){
			layout(rbind(c(1,1,2,2),c(3,3,4,4)))
			par(oma=c(0,0,3,0))
			par(mar=c(4,4,8.3,2))
	
			if(boxplot){
				boxplot(sens ~ spec_1, xlab ="1-Specificity", ylab="Sensitivity", main="ROC curve of the k-TSP with M=0.25", col=box.col, at=at1, xlim=c(0,1), ylim=c(0,1), boxwex = 0.03, axes=FALSE)

				axis(side = 1, at = c(seq(0,1,0.2)), labels = c(seq(0,1,0.2)))
				axis(side = 2, at = c(seq(0,1,0.2)), labels = c(seq(0,1,0.2)))
				abline(0,1, col="gray60")
	
				points(at1,median1, col=line.col, pch=20)
				lines(at1, median1, col=line.col, lwd=2)
			}
			else{
				plot(spec_1, sens, xlim=c(0,1), ylim=c(0,1), xlab ="1-Specificity", ylab="Sensitivity", main="ROC curve of the k-TSP with M=0.25", axes=FALSE)
				axis(side = 1, at = c(seq(0,1,0.2)), labels = c(seq(0,1,0.2)))
				axis(side = 2, at = c(seq(0,1,0.2)), labels = c(seq(0,1,0.2)))
				abline(0,1, col="gray60")
			}
		}

		if(AUC){
			auc1 <- round(AUC.calc(1-at1,median1), digits=2)
			if(graphic){
				legend(auc.x,auc.y,bquote(median  ==.(auc1)), title="AUC:", title.col="blue", text.col=c("red"), cex = 0.8, box.col="orange")
			}

		}
		spec_5 <- list(length=n)
		sens_5 <- list(length=n)

		for(i in 1:n){
			spec_5[[i]] <- unique(spec_list[[2]][i,])
			sens_5[[i]] <- c()
			mean <- c()
			for(j in 1:length(spec_5[[i]])){
				mean[j] <- mean(sens_list[[2]][i,][spec_list[[2]][i,]==spec_5[[i]][j]], na.rm=TRUE)
			}
			sens_5[[i]] <- mean
		}	

		spec_5_1 <- 1-unlist(spec_5)
		spec_5_1_order <- spec_5_1[order(spec_5_1)]
		sens_5_order <- unlist(sens_5)[order(spec_5_1)]

		uniquespec_5_1 <- unique(spec_5_1_order)
		box <- list()
		for(i in 1:length(uniquespec_5_1)){
			box[[i]]<- sens_5_order[which(spec_5_1_order==uniquespec_5_1[i])]
		}

		spec_15 <-c()
		sens5 <- c()
		median5 <- c()
		at5 <- c()

		for(i in 1:(length(box)/m)){
			sens5 <- c(sens5, box[[m*i]])
			spec_15 <- c(spec_15, rep(uniquespec_5_1[m*i],length(box[[m*i]])))
			median5[i] <- median(box[[m*i]])
			at5[i] <- uniquespec_5_1[m*i]
		}

		if(graphic){
			if(boxplot){
				boxplot(sens5 ~ spec_15, xlab ="1-Specificity", ylab="Sensitivity", main="ROC curve of the k-TSP with M=0.5", col=box.col, at=at5, xlim=c(0,1), ylim=c(0,1), boxwex = 0.03, axes=FALSE)
				axis(side = 1, at = c(seq(0,1,0.2)), labels = c(seq(0,1,0.2)))
				axis(side = 2, at = c(seq(0,1,0.2)), labels = c(seq(0,1,0.2)))
				abline(0,1, col="gray60")
	
				points(at5,median5, col=line.col, pch=20)
				lines(at5, median5, col=line.col, lwd=2)
			}

			else{
				plot(spec_15, sens5, xlim=c(0,1), ylim=c(0,1), xlab ="1-Specificity", ylab="Sensitivity", main="ROC curve of the k-TSP with M=0.5", axes=FALSE)
				axis(side = 1, at = c(seq(0,1,0.2)), labels = c(seq(0,1,0.2)))
				axis(side = 2, at = c(seq(0,1,0.2)), labels = c(seq(0,1,0.2)))
				abline(0,1, col="gray60")
			}
		}
		if(AUC){
			auc5 <- round(AUC.calc(1-at5,median5), digits=2)
			if(graphic){
				legend(auc.x,auc.y,bquote(median  ==.(auc5)), title="AUC:", title.col="blue", text.col=c("red"), cex = 0.8, box.col="orange")
			}

		}
		spec_7 <- list(length=n)
		sens_7 <- list(length=n)

		for(i in 1:n){
			spec_7[[i]] <- unique(spec_list[[3]][i,])
			sens_7[[i]] <- c()
			mean <- c()
			for(j in 1:length(spec_7[[i]])){
				mean[j] <- mean(sens_list[[3]][i,][spec_list[[3]][i,]==spec_7[[i]][j]], na.rm=TRUE)
			}
			sens_7[[i]] <- mean
		}

		spec_7_1 <- 1-unlist(spec_7)

		spec_7_1_order <- spec_7_1[order(spec_7_1)]
		sens_7_order <- unlist(sens_7)[order(spec_7_1)]

		uniquespec_7_1 <- unique(spec_7_1_order)
		box <- list()
		for(i in 1:length(uniquespec_7_1)){
			box[[i]]<- sens_7_order[which(spec_7_1_order==uniquespec_7_1[i])]
		}

		spec_17 <-c()
		sens7 <- c()
		median7 <- c()
		at7 <- c()

		for(i in 1:(length(box)/m)){
			sens7 <- c(sens7, box[[m*i]])
			spec_17 <- c(spec_17, rep(uniquespec_7_1[m*i],length(box[[m*i]])))
			median7[i] <- median(box[[m*i]])
			at7[i] <- uniquespec_7_1[m*i]
		}

		if(graphic){
			if(boxplot){
				boxplot(sens7 ~ spec_17, xlab ="1-Specificity", ylab="Sensitivity", main="ROC curve of the k-TSP with M=0.75", col=box.col, at=at7, xlim=c(0,1), ylim=c(0,1), boxwex = 0.03, axes=FALSE)
				axis(side = 1, at = c(seq(0,1,0.2)), labels = c(seq(0,1,0.2)))
				axis(side = 2, at = c(seq(0,1,0.2)), labels = c(seq(0,1,0.2)))
				abline(0,1, col="gray60")

				points(at7, median7, col=line.col, pch=20)
				lines(at7, median7, col=line.col, lwd=2)
			}

			else{
				plot(spec_17, sens7, xlim=c(0,1), ylim=c(0,1), xlab ="1-Specificity", ylab="Sensitivity", main="ROC curve of the k-TSP with M=0.75", axes=FALSE)
				axis(side = 1, at = c(seq(0,1,0.2)), labels = c(seq(0,1,0.2)))
				axis(side = 2, at = c(seq(0,1,0.2)), labels = c(seq(0,1,0.2)))
				abline(0,1, col="gray60")
			}
		}

		if(AUC){
			auc7 <- round(AUC.calc(1-at7,median7), digits=2)
			if(graphic){
				legend(auc.x,auc.y,bquote(median  ==.(auc7)), title="AUC:", title.col="blue", text.col=c("red"), cex = 0.8, box.col="orange")
			}
		}
		if(graphic){

			plot(at1,median1, pch=20, xlab="1-Specificity", xlim=c(0,1), ylim=c(0,1), ylab="Median of the sensitivity", main="Summary of the ROC curves
for the different thresholds")
   
			points(at1,median1, col=multiple.col[1], pch=20)
			lines(at1, median1, col=multiple.col[1])

			points(at5,median5, col=multiple.col[2], pch=20)
			lines(at5, median5, col=multiple.col[2])

			points(at7,median7, col=multiple.col[3], pch=20)
			lines(at7, median7, col=multiple.col[3])

			abline(0, 1, col="gray60")
			legend(0.75, 0.3, c("M=0.25","M=0.5", "M=0.75"), cex=0.8, 
			col=multiple.col, pch=20, lty=1)

			if(!is.null(maintitle)){
				title(maintitle, line=0, outer=TRUE, cex.main=2)
			}
			if(!is.null(mtext)){
				mtext(mtext, side=3, line=-2.5, outer=TRUE, cex.main=1.5)
			}

		}
	}

	if(mult.cutoff==FALSE){
		spec <- list(length=n)
		sens <- list(length=n)

		for(i in 1:n){
			spec[[i]] <- unique(spec_list[i,])
			sens[[i]] <- c()
			mean <- c()
			for(j in 1:length(spec[[i]])){
				mean[j] <- mean(sens_list[i,][spec_list[i,]==spec[[i]][j]], na.rm=TRUE)
			}
			sens[[i]] <- mean
		}

		spec_1 <- 1-unlist(spec)

		spec_1_order <- spec_1[order(spec_1)]
		sens_order <- unlist(sens)[order(spec_1)]

		uniquespec_1 <- unique(spec_1_order)
		box <- list()
		for(i in 1:length(uniquespec_1)){
			box[[i]]<- sens_order[which(spec_1_order==uniquespec_1[i])]
		}

		spec_1 <-c()
		sens <- c()
		median <- c()
		at <- c()

		for(i in 1:(length(box)/m)){
			sens <- c(sens, box[[m*i]])
			spec_1 <- c(spec_1, rep(uniquespec_1[m*i],length(box[[m*i]])))
			median[i] <- median(box[[m*i]])
			at[i] <- uniquespec_1[m*i]
		}

		if(graphic){
			if(boxplot){
				if(is.null(undertitle)){
					boxplot(sens ~ spec_1, xlab ="1-Specificity", ylab="Sensitivity", main="", col=box.col, at=at, xlim=c(0,1), ylim=c(0,1), boxwex = 0.03, axes=FALSE)
				}
				else{
					boxplot(sens ~ spec_1, xlab ="1-Specificity", ylab="Sensitivity", main=undertitle, col=box.col, at=at, xlim=c(0,1), ylim=c(0,1), boxwex = 0.03, axes=FALSE,line=1)
				}
				axis(side = 1, at = c(seq(0,1,0.2)), labels = c(seq(0,1,0.2)))
				axis(side = 2, at = c(seq(0,1,0.2)), labels = c(seq(0,1,0.2)))
				abline(0,1, col="gray60")
	
				points(at,median, col=line.col, pch=20)
				lines(at, median, col=line.col, lwd=2)
			}
			else{
				if(is.null(undertitle)){
					plot(spec_1, sens, xlim=c(0,1), ylim=c(0,1), xlab ="1-Specificity", ylab="Sensitivity", main="", axes=FALSE)
				}
				else{
					plot(spec_1, sens, xlim=c(0,1), ylim=c(0,1), xlab ="1-Specificity", ylab="Sensitivity", main=undertitle, axes=FALSE)
				}
				axis(side = 1, at = c(seq(0,1,0.2)), labels = c(seq(0,1,0.2)))
				axis(side = 2, at = c(seq(0,1,0.2)), labels = c(seq(0,1,0.2)))
				abline(0,1, col="gray60")
			}

			if(AUC){
				auc <- round(AUC.calc(1-at,median[order(median)]), digits=2)
				legend(auc.x,auc.y,bquote(median  ==.(auc)), title="AUC:", title.col="blue", text.col=c("red"), cex = 0.8, 	box.col="orange")
			}
			if(!is.null(maintitle)){
				title(maintitle, line=-2, outer=TRUE, cex.main=2)
			}
			if(!is.null(mtext)){
				mtext(mtext, side=3, line=-3.5, outer=TRUE, cex.main=1.5)
			}
		}
		if(AUC){
			auc <- round(AUC.calc(1-at,median[order(median)]), digits=2)
			if(graphic){legend(auc.x,auc.y,bquote(median  ==.(auc)), title="AUC:", title.col="blue", text.col=c("red"), cex = 0.8, 	box.col="orange")}
		}
	
	}

	if(mult.cutoff){
		at<-list(at1, at5, at7)
		median <- list(median1, median5, median7)
		if(AUC){
			auc <- list(auc1, auc5, auc7)
		}
	}	
	ROC.graphic <- list(at = at, median = median, auc = auc)
	class(ROC.graphic) <- "ROC.graphic"
	return(ROC.graphic)
}


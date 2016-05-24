
`summary.HTSCluster` <-
function (object, ...) 
{
	x <- object
    	if (class(x) != "HTSCluster") {
        	stop(paste(sQuote("x"), sep = ""), " must be of class ", 
            paste(dQuote("HTSCluster"), sep = ""), sep = "")
    	}

	probaPost <- x$probaPost
	labels <- x$labels
	lambda <- x$lambda
	pi <- x$pi
	g <- length(pi)

	map <- apply(probaPost, 1, max)
	length(which(map > 0.9))/length(map)

	cat("*************************************************\n")
	cat("Number of clusters = ", g, "\n", sep = "")
	if(is.na(x$model.selection) == FALSE) {
		cat("Model selection via ", x$model.selection, "\n", sep = "")
	}
	cat("*************************************************\n")
	tab <- table(labels)
	names(tab) <- paste("Cluster", names(tab))
	cat("Cluster sizes:\n"); print(tab); cat("\n")
	cat("Number of observations with MAP > 0.90 (% of total):\n")
	cat(length(which(map > 0.9)), " (", round(length(which(map > 0.9))/length(map)*100,2),
		"%)\n\n", sep = "")
	cat("Number of observations with MAP > 0.90 per cluster (% of total per cluster):\n"); 

	tab2 <- matrix(NA, nrow = 2, ncol = g)
	colnames(tab2) <- paste("Cluster", 1:g); rownames(tab2) <- rep("", 2)
	for(i in 1:g) {
		if(sum(labels == i) > 1) {
			map.clust <- apply(matrix(probaPost[labels == i,], ncol=g), 1, max)
			tab2[1,i] <- length(which(map.clust > 0.9))
			tab2[2,i] <- paste("(", round(100*length(which(map.clust > 0.9))/length(map.clust),2),
				"%)", sep = "")
		}
		if(sum(labels == i) == 1) {
			map.clust <- max(probaPost[labels == i,])
			tab2[1,i] <- length(which(map.clust > 0.9))
			tab2[2,i] <- paste("(", round(100*length(which(map.clust > 0.9))/length(map.clust),2),
				"%)", sep = "")
		}
		if(sum(labels == i) == 0) {
			tab2[1,i] <- "---"
			tab2[2,i] <- "---"
		}
	}
	print(tab2, quote = FALSE); cat("\n")

	cat("Lambda:\n"); print(round(lambda,2)); cat("\n")
	cat("Pi:\n"); print(round(pi,2)); cat("\n")
}




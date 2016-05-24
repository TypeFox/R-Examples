#################################################################################
#' CDROM: Classification of Duplicate gene RetentiOn Mechanisms
#'
#' Classification is based on the recently developed phylogenetic approach by Assis and Bachtrog (2013). The method classifies the evolutionary mechanisms
#' retaining pairs of duplicate genes (conservation, neofunctionalization, subfunctionalization, or specialization) by comparing gene expression profiles of 
#' duplicate genes in one species to those of their single-copy ancestral genes in a sister species. 
#'
#' @param dupFile a tab-separated file containing duplicate gene pairs and their orthologs (three gene IDs per line).
#' @param singleFile a tab-separated file containing orthologous single-copy genes (two gene IDs per line).
#' @param exprFile1 a tab-separated file containing gene expression levels for all genes in species 1 (one gene ID and its expression levels per line).
#' @param exprFile2 a tab-separated file containing gene expression levels for all genes in species 2 (one gene ID and its expression levels per line).
#' @param out the prefix to be used in the names of the three output files. Defaults to 'out'.
#' @param PC a logical value indicating whether parent and child copies are separated in dupFile. Defaults to FALSE.
#' @param Ediv the divergence cutoff to be used in classifications. Defaults to semi-interquartile range of the median (SIQR).
#' @param useAbsExpr a logical value indicating whether absolute or relative expression values are used for Euclidean distance calculations. Defaults to FALSE.
#' @param head1 a logical value indicating whether exprFile1 contains the names of its variables as the first line. Defaults to TRUE.
#' @param head2 a logical value indicating whether exprFile2 contains the names of its variables as the first line. Defaults to TRUE.
#' @param head3 a logical value indicating whether dupFile contains the names of its variables as the first line. Defaults to TRUE.
#' @param head4 a logical value indicating whether singleFile contains the names of its variables as the first line. Defaults to TRUE.
#' @param legend a keyword indicating the position of the legend in the output plot. Options are 'topright' and 'topleft'. Defaults to 'topleft'.
#' @return A table with the classifications of all duplicate gene pairs.
#' @return A table with counts of classifications with five Ediv values.
#' @return A plot showing distributions of all Euclidean distances and the position of Ediv.
#' @author Brent Perry <brp5173@psu.edu>, Raquel Assis <rassis@psu.edu>
#' @importFrom grDevices dev.off png
#' @importFrom graphics abline lines par plot title
#' @importFrom stats IQR density median quantile sd
#' @importFrom utils read.table write.table
#' @examples CDROM(dupFile=system.file("extdata","human_chicken_dups.txt",package="CDROM"),
#' singleFile=system.file("extdata","human_chicken_singles.txt",package="CDROM"),
#' exprFile1=system.file("extdata","human_expr.txt",package="CDROM"),
#' exprFile2=system.file("extdata","chicken_expr.txt",package="CDROM"))
#' @export
#' 
#'
#' @keywords gene duplication, neofunctionalization, subfunctionalization
#' @references 
#' [1] Assis R and Bachtrog D. Neofunctionalization of young duplicate genes in Drosophila. Proc. Natl. Acad. Sci. USA. 110: 17409-17414 (2013). 
#################################################################################

CDROM <- function(dupFile, singleFile, exprFile1, exprFile2, out = "out", Ediv,  
	PC = FALSE, useAbsExpr = FALSE, head1 = TRUE, head2 = TRUE, head3 = TRUE, head4 = TRUE, legend = "topleft") {
	

	## General checks are made
    
	if (missing(exprFile1)) 
		stop("'exprFile1' is missing")
	if (missing(exprFile2)) 
		stop("'exprFile2' is missing")
	if (missing(dupFile)) 
		stop("'dupFile' is missing")
	if (missing(singleFile)) 
		stop("'singleFile' is missing")
	if (missing(Ediv)) 
		cat("Note: 'SIQR' will be used as Ediv\n")
	if(! any(legend == c("topleft", "topright"))) {
		legend <- "topleft"
		cat("Legend will be set to default position (topleft)\n")
	}
    

	## Input files are read

	if (head1 == TRUE) {
		expr1 <- read.table(exprFile1, row.names = 1, header = TRUE)
	} else {
		expr1 <- read.table(exprFile1, row.names = 1, header = FALSE)
	}
    
	if (head2 == TRUE) {
		expr2 <- read.table(exprFile2, row.names = 1, header = TRUE)
	} else {
		expr2 <- read.table(exprFile2, row.names = 1, header = FALSE)
	}

	if (head3 == TRUE) {
		dups <- read.table(dupFile, header = TRUE)
	} else {
		dups <- read.table(dupFile, header = FALSE)
	}
    
	if (head4 == TRUE) {
		singles <- read.table(singleFile, header = TRUE)
	} else {
		singles <- read.table(singleFile, header = FALSE)
	}

    
	## Expression data are obtained from expression files

	colnames(expr2) <- colnames(expr1)
	expr1$zeros <- 0
	expr2$zeros <- 0
	expr <- rbind(expr1[! (row.names(expr1) %in% row.names(expr2)), ], expr2)

	P <- dups[[1]]
	C <- dups[[2]]
	A <- dups[[3]]
	S1 <- singles[[1]]
	S2 <- singles[[2]]
  
	getP <- expr[row.names(expr) %in% P, ]
	getC <- expr[row.names(expr) %in% C, ]
	getA <- expr[row.names(expr) %in% A, ]
	getS1 <- expr[row.names(expr) %in% S1, ]
	getS2 <- expr[row.names(expr) %in% S2, ]
    
	exprP <- data.matrix(getP[match(P, rownames(getP)), ])
	exprC <- data.matrix(getC[match(C, rownames(getC)), ])
	exprA <- data.matrix(getA[match(A, rownames(getA)), ])
	exprS1 <- data.matrix(getS1[match(S1, rownames(getS1)), ])
	exprS2 <- data.matrix(getS2[match(S2, rownames(getS2)), ])
	exprPC <- exprP + exprC
 
  if (useAbsExpr == FALSE) {
	 
	  ## Relative expression values are calculated 
      
	  relP <- exprP / rowSums(exprP)
	  relC <- exprC / rowSums(exprC)
	  relA <- exprA / rowSums(exprA)
	  relPC <- exprPC / rowSums(exprPC)
	  relS1 <- exprS1 / rowSums(exprS1)
	  relS2 <- exprS2 / rowSums(exprS2)
    
	  ## Euclidean distances are calculated
    
	  eucPA <- (rowSums((relP - relA) ^ 2)) ^ (1/2)
	  eucCA <- (rowSums((relC - relA) ^ 2)) ^ (1/2)
	  eucPCA <- (rowSums((relPC - relA) ^ 2)) ^ (1/2)
	  eucS1S2 <- (rowSums((relS1 - relS2) ^ 2)) ^ (1/2)
  }
  
	if (useAbsExpr == TRUE) {
	  
	  ## Euclidean distances are calculated
	  
	  eucPA <- (rowSums((exprP - exprA) ^ 2)) ^ (1/2)
	  eucCA <- (rowSums((exprC - exprA) ^ 2)) ^ (1/2)
	  eucPCA <- (rowSums((exprPC - exprA) ^ 2)) ^ (1/2)
	  eucS1S2 <- (rowSums((exprS1 - exprS2) ^ 2)) ^ (1/2)
	}

	## Ediv is calculated 
    
	if (missing(Ediv)) {
		SIQR <- (IQR(eucS1S2, na.rm = TRUE) / 2)
		Ediv <- median(eucS1S2, na.rm = TRUE) + SIQR
	} else {
		Ediv
		if (Ediv < 0) 
		stop("'Ediv' must be greater than or equal to zero\n")
	}

	## Classifications are made

	eucDists <- data.frame((eucPA), (eucCA), (eucPCA))
	eucDists <- replace(eucDists, is.na(eucDists), "NA")
    
	if (PC == FALSE) {
		eucDists <- within(eucDists, Classification <- "NA")
		eucDists[(eucDists$X.eucPA. <= Ediv & eucDists$X.eucCA. <= 
			Ediv), "Classification"] <- "Conservation"
		eucDists[(eucDists$X.eucPA. > Ediv & eucDists$X.eucCA. <= 
			Ediv), "Classification"] <- "Neofunctionalization(Dup1)"
		eucDists[(eucDists$X.eucPA. <= Ediv & eucDists$X.eucCA. > 
			Ediv), "Classification"] <- "Neofunctionalization(Dup2)"
		eucDists[(eucDists$X.eucPA. > Ediv & eucDists$X.eucCA. > 
			Ediv & eucDists$X.eucPCA. <= Ediv), "Classification"] <- "Subfunctionalization"
		eucDists[(eucDists$X.eucPA. > Ediv & eucDists$X.eucCA. > 
			Ediv & eucDists$X.eucPCA. > Ediv), "Classification"] <- "Specialization"
		eucDists[(eucDists$X.eucPA. == "NA" | eucDists$X.eucCA. == "NA" | 
			eucDists$X.eucPCA. == "NA"), "Classification"] <- "NA"


		## Densities of Euclidian distances are plotted 
    
		densS1S2 <- density(eucS1S2, na.rm = TRUE)
		Ymax_S1S2 <- max(densS1S2$y)
		Xmax_S1S2 <- max(densS1S2$x)
		Xmin_S1S2 <- min(densS1S2$x)
        
		densBoth <- density(c(eucCA, eucPA), na.rm = TRUE)
		Ymax_Both <- max(densBoth$y)
		Xmax_Both <- max(densBoth$x)
		Xmin_Both <- min(densBoth$x)

		densPCA <- density(eucPCA, na.rm = TRUE)
		Ymax_PCA <- max(densPCA$y)
		Xmax_PCA <- max(densPCA$x)
		Xmin_PCA <- min(densPCA$x)

		Ymax <- max(Ymax_S1S2, Ymax_Both, Ymax_PCA)
		Xmax <- max(Xmax_S1S2, Xmax_Both, Xmax_PCA)
		Xmin <- min(Xmin_S1S2, Xmin_Both, Xmin_PCA)

		Yrange <- (Ymax - 0)
		Xrange <- (Xmax - Xmin)

		if (useAbsExpr == FALSE) {
			upperY <- (Ymax + (0.05 * Yrange))
			lowerX <- (Xmin - (0.02 * Xrange))
			upperX <- (Xmax + (0.02 * Xrange))
		}

		if (useAbsExpr == TRUE) {
			upperY <- (Ymax + (0.05 * Yrange))			
			upperX <- (20 * IQR(eucS1S2, na.rm = TRUE))
			lowerX <- -IQR(eucS1S2, na.rm = TRUE)
		}

	    png(filename = paste0(out, ".png"), width = 700, height = 700)
    	par(mar = c(4, 5, 1, 1), cex.axis = 1.4, cex.lab = 2)
		plot(densS1S2, col = "black", main = "", xlab = "Euclidean Distance", ylab = "", 
			lwd = 6, mgp = c(3, 1, 0), yaxt = "n", ylim = c(0, upperY), 
			xlim = c(lowerX, upperX), xaxs = "i", yaxs = "i")
		lines(densBoth, col = "green3", lwd = 6)
		lines(densPCA, col = "purple3", lwd = 6)
		abline(v = Ediv, col = "grey50", lty = 2, lwd = 4)

		if (legend == "topright") {
			legend("topright", c(expression(italic("E")["S1,S2"]), expression(italic("E")["D1,A"]+italic("E")["D2,A"]),
				expression(italic("E")["D1+D2,A"]), paste("Ediv = ",round(Ediv, 4))), 
				fill = c("black", "green3", "purple3", "grey50"), cex = 2, pt.cex = 1.1)
		} else {
			legend("topleft", c(expression(italic("E")["S1,S2"]), expression(italic("E")["D1,A"]+italic("E")["D2,A"]),
				expression(italic("E")["D1+D2,A"]), paste("Ediv = ",round(Ediv, 4))), 
				fill = c("black", "green3", "purple3", "grey50"), cex = 2, pt.cex = 1.1)
		}
    
		title(ylab = "Density", mgp = c(1.6, 1, 0), cex.lab = 2.2)
		dev.off()
    

	## Output file1 is written
    
		classes <- data.frame(P, C, A, eucDists)
		rownames(classes) <- NULL
		colnames(classes) <- c("Dup1", "Dup2", "Ancestor", "E_D1,A", "E_D2,A", 
			"E_D1+D2,A", "Classification")
		sink(paste0(out, "1.txt"))
		write.table(classes, quote = FALSE, sep = "\t", row.names = FALSE)
		sink()


		## Five E_div values are calculated

		meanSD <- mean(eucS1S2, na.rm = TRUE) + sd(eucS1S2, na.rm = TRUE)
		mean2SD <- mean(eucS1S2, na.rm = TRUE) + 2*sd(eucS1S2, na.rm = TRUE)
		medSIQR <- median(eucS1S2, na.rm = TRUE) + (IQR(eucS1S2, na.rm = TRUE) / 2)
		medIQR <- median(eucS1S2, na.rm = TRUE) + IQR(eucS1S2, na.rm = TRUE)
		quant75 <- quantile(eucS1S2, 0.75, na.rm = TRUE)


		## Output2 file with five E_div values is generated

		Edivs <- (c(meanSD,mean2SD,medSIQR,medIQR,quant75))
		classify <- matrix(nrow = 5, ncol = 5)
		eucDists <- data.frame((eucPA), (eucCA), (eucPCA))

		for ( i in 1:5){

			Ediv <- Edivs[i]
			eucDists <- replace(eucDists, is.na(eucDists), "NA")
		
			eucDists <- within(eucDists, Classification <- "NA")
			eucDists[(eucDists$X.eucPA. <= Ediv & eucDists$X.eucCA. <= 
				Ediv), "Classification"] <- "Conservation"
			eucDists[(eucDists$X.eucPA. > Ediv & eucDists$X.eucCA. <= 
				Ediv), "Classification"] <- "Neofunctionalization(Dup1)"
			eucDists[(eucDists$X.eucPA. <= Ediv & eucDists$X.eucCA. > 
				Ediv), "Classification"] <- "Neofunctionalization(Dup2)"
			eucDists[(eucDists$X.eucPA. > Ediv & eucDists$X.eucCA. > 
				Ediv & eucDists$X.eucPCA. <= Ediv), "Classification"] <- "Subfunctionalization"
			eucDists[(eucDists$X.eucPA. > Ediv & eucDists$X.eucCA. > 
				Ediv & eucDists$X.eucPCA. > Ediv), "Classification"] <- "Specialization"
			eucDists[(eucDists$X.eucPA. == "NA" | eucDists$X.eucCA. == "NA" | 
				eucDists$X.eucPCA. == "NA"), "Classification"] <- "NA"

			counts <- as.data.frame(table(eucDists[[4]]))
			countsOrdered <- counts[match(c("Conservation","Neofunctionalization(Dup1)",
				"Neofunctionalization(Dup2)","Subfunctionalization","Specialization"), counts[[1]]),]
			classify[i,] <- countsOrdered[[2]] 
		}

		classDF <- as.data.frame(classify)
		classDF[6] <- c("meanSD","mean2SD","medSIQR","medIQR","quant75")
		names(classDF)[6] <- "Ediv"

		idDF <- data.frame(Edivs)
		idDF[1] <- c("meanSD","mean2SD","medSIQR","medIQR","quant75")
		idDF[2] <- c(meanSD,mean2SD,medSIQR,medIQR,quant75)

		names(idDF)[1] <- "Ediv"

		total <- merge(idDF,classDF,by="Ediv")
		colnames(total) <- c("E_div","Value","Conservation","Neofunctionalization(Dup1)",
			"Neofunctionalization(Dup2)","Subfunctionalization","Specialization")
		total[is.na(total)] <- 0
		totalOrdered <- total[match(c("meanSD","mean2SD","medSIQR","medIQR","quant75"), total[[1]]),]

		sink(paste0(out, "2.txt"))
		write.table(totalOrdered, quote = FALSE, sep = "\t", row.names = FALSE)
		sink()

	}

	if (PC == TRUE) {
		eucDists <- within(eucDists, Classification <- "NA")
		eucDists[(eucDists$X.eucPA. <= Ediv & eucDists$X.eucCA. <= 
			Ediv), "Classification"] <- "Conservation"
		eucDists[(eucDists$X.eucPA. > Ediv & eucDists$X.eucCA. <= 
			Ediv), "Classification"] <- "Neofunctionalization(Parent)"
		eucDists[(eucDists$X.eucPA. <= Ediv & eucDists$X.eucCA. > 
			Ediv), "Classification"] <- "Neofunctionalization(Child)"
		eucDists[(eucDists$X.eucPA. > Ediv & eucDists$X.eucCA. > 
			Ediv & eucDists$X.eucPCA. <= Ediv), "Classification"] <- "Subfunctionalization"
		eucDists[(eucDists$X.eucPA. > Ediv & eucDists$X.eucCA. > 
			Ediv & eucDists$X.eucPCA. > Ediv), "Classification"] <- "Specialization"
		eucDists[(eucDists$X.eucPA. == "NA" | eucDists$X.eucCA. == "NA" | 
			eucDists$X.eucPCA. == "NA"), "Classification"] <- "NA"


		## Densities of Euclidian distances are plotted 
    
		densS1S2 <- density(eucS1S2, na.rm = TRUE)
		Ymax_S1S2 <- max(densS1S2$y)
		Xmax_S1S2 <- max(densS1S2$x)
		Xmin_S1S2 <- min(densS1S2$x)

		densCA <- density(eucCA, na.rm = TRUE)
		Ymax_CA <- max(densCA$y)
		Xmax_CA <- max(densCA$x)
		Xmin_CA <- min(densCA$x)

		densPA <- density(eucPA, na.rm = TRUE)
		Ymax_PA <- max(densPA$y)
		Xmax_PA <- max(densPA$x)
		Xmin_PA <- min(densPA$x)

		densPCA <- density(eucPCA, na.rm = TRUE)
		Ymax_PCA <- max(densPCA$y)
		Xmax_PCA <- max(densPCA$x)
		Xmin_PCA <- min(densPCA$x)

		Ymax <- max(Ymax_S1S2, Ymax_CA, Ymax_PA, Ymax_PCA)
		Xmax <- max(Xmax_S1S2, Xmax_CA, Xmax_PA, Xmax_PCA)
		Xmin <- min(Xmin_S1S2, Xmin_CA, Xmin_PA, Xmin_PCA)

		Yrange <- (Ymax - 0)
		Xrange <- (Xmax - Xmin)
		
		if (useAbsExpr == FALSE) {
			upperY <- (Ymax + (0.05 * Yrange))
			lowerX <- (Xmin - (0.02 * Xrange))
			upperX <- (Xmax + (0.02 * Xrange))
		}

		if (useAbsExpr == TRUE) {
			upperY <- (Ymax + (0.05 * Yrange))			
			upperX <- (20 * IQR(eucS1S2, na.rm = TRUE))
			lowerX <- -IQR(eucS1S2, na.rm = TRUE)
		}

    	png(filename = paste0(out, ".png"), width = 700, height = 700)
    	par(mar = c(4, 5, 1, 1), cex.axis = 1.4, cex.lab = 2)
		plot(densS1S2, col = "black", main = "", xlab = "Euclidean Distance", ylab = "", 
			lwd = 6, mgp = c(3, 1, 0), yaxt = "n", ylim = c(0, upperY), 
			xlim = c(lowerX, upperX), xaxs = "i", yaxs = "i")
		lines(densCA, col = "red", lwd = 6)
		lines(densPA, col = "blue3", lwd = 6)
		lines(densPCA, col = "purple3", lwd = 6)
		abline(v = Ediv, col = "grey50", lty = 2, lwd = 4)

		if (legend == "topright") {
			legend("topright", c(expression(italic("E")["S1,S2"]), expression(italic("E")["C,A"]), 
				expression(italic("E")["P,A"]), expression(italic("E")["P+C,A"]), paste("Ediv = ",round(Ediv, 4))), 
				fill = c("black", "red", "blue3", "purple3", "grey50"), cex = 2, pt.cex = 1.1)
		} else {
		legend("topleft", c(expression(italic("E")["S1,S2"]), expression(italic("E")["C,A"]), 
				expression(italic("E")["P,A"]), expression(italic("E")["P+C,A"]), paste("Ediv = ",round(Ediv, 4))), 
				fill = c("black", "red", "blue3", "purple3", "grey50"), cex = 2, pt.cex = 1.1)
		}
    
		title(ylab = "Density", mgp = c(1.6, 1, 0), cex.lab = 2.2)
		dev.off()
    

		## Output file1 is written
    
		classes <- data.frame(P, C, A, eucDists)
		rownames(classes) <- NULL
		colnames(classes) <- c("Parent", "Child", "Ancestor", "E_P,A", "E_C,A", 
			"E_P+C,A", "Classification")
		sink(paste0(out, "1.txt"))
		write.table(classes, quote = FALSE, sep = "\t", row.names = FALSE)
		sink()


		## Five E_div values are calculated

		meanSD <- mean(eucS1S2, na.rm = TRUE) + sd(eucS1S2, na.rm = TRUE)
		mean2SD <- mean(eucS1S2, na.rm = TRUE) + 2*sd(eucS1S2, na.rm = TRUE)
		medSIQR <- median(eucS1S2, na.rm = TRUE) + (IQR(eucS1S2, na.rm = TRUE) / 2)
		medIQR <- median(eucS1S2, na.rm = TRUE) + IQR(eucS1S2, na.rm = TRUE)
		quant75 <- quantile(eucS1S2, 0.75, na.rm = TRUE)


		## Output2 file with five E_div values is generated

		Edivs <- (c(meanSD,mean2SD,medSIQR,medIQR,quant75))
		classify <- matrix(nrow = 5, ncol = 5)
		eucDists <- data.frame((eucPA), (eucCA), (eucPCA))

		for ( i in 1:5){

			Ediv <- Edivs[i]
			eucDists <- replace(eucDists, is.na(eucDists), "NA")
		
			eucDists <- within(eucDists, Classification <- "NA")
			eucDists[(eucDists$X.eucPA. <= Ediv & eucDists$X.eucCA. <= 
				Ediv), "Classification"] <- "Conservation"
			eucDists[(eucDists$X.eucPA. > Ediv & eucDists$X.eucCA. <= 
				Ediv), "Classification"] <- "Neofunctionalization(Parent)"
			eucDists[(eucDists$X.eucPA. <= Ediv & eucDists$X.eucCA. > 
				Ediv), "Classification"] <- "Neofunctionalization(Child)"
			eucDists[(eucDists$X.eucPA. > Ediv & eucDists$X.eucCA. > 
				Ediv & eucDists$X.eucPCA. <= Ediv), "Classification"] <- "Subfunctionalization"
			eucDists[(eucDists$X.eucPA. > Ediv & eucDists$X.eucCA. > 
				Ediv & eucDists$X.eucPCA. > Ediv), "Classification"] <- "Specialization"
			eucDists[(eucDists$X.eucPA. == "NA" | eucDists$X.eucCA. == "NA" | 
				eucDists$X.eucPCA. == "NA"), "Classification"] <- "NA"

			counts <- as.data.frame(table(eucDists[[4]]))
			countsOrdered <- counts[match(c("Conservation","Neofunctionalization(Parent)",
				"Neofunctionalization(Child)","Subfunctionalization","Specialization"), counts[[1]]),]
			classify[i,] <- countsOrdered[[2]] 
		}

		classDF <- as.data.frame(classify)
		classDF[6] <- c("meanSD","mean2SD","medSIQR","medIQR","quant75")
		names(classDF)[6] <- "Ediv"

		idDF <- data.frame(Edivs)
		idDF[1] <- c("meanSD","mean2SD","medSIQR","medIQR","quant75")
		idDF[2] <- c(meanSD,mean2SD,medSIQR,medIQR,quant75)

		names(idDF)[1] <- "Ediv"

		total <- merge(idDF,classDF,by="Ediv")
		colnames(total) <- c("E_div","Value","Conservation","Neofunctionalization(Parent)",
			"Neofunctionalization(Child)","Subfunctionalization","Specialization")
		total[is.na(total)] <- 0
		totalOrdered <- total[match(c("meanSD","mean2SD","medSIQR","medIQR","quant75"), total[[1]]),]

		sink(paste0(out, "2.txt"))
		write.table(totalOrdered, quote = FALSE, sep = "\t", row.names = FALSE)
		sink()
	}
} 

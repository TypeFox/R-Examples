`LSD.test` <-
		function (y, trt, DFerror, MSerror, alpha = 0.05, p.adj = c("none",
						"holm", "hochberg", "bonferroni", "BH", "BY", "fdr"), group = TRUE,
				main = NULL,console=FALSE)
{
	p.adj <- match.arg(p.adj)
	clase <- c("aov", "lm")
	name.y <- paste(deparse(substitute(y)))
	name.t <- paste(deparse(substitute(trt)))
	if(is.null(main))main<-paste(name.y,"~", name.t)
	if ("aov" %in% class(y) | "lm" %in% class(y)) {
		if(is.null(main))main<-y$call
		A <- y$model
		DFerror <- df.residual(y)
		MSerror <- deviance(y)/DFerror
		y <- A[, 1]
		ipch <- pmatch(trt, names(A))
		nipch<- length(ipch)
		for(i in 1:nipch){
			if (is.na(ipch[i]))
				return(if(console)cat("Name: ", trt, "\n", names(A)[-1], "\n"))
		}
		name.t<- names(A)[ipch][1]
		trt <- A[, ipch]
		if (nipch > 1){
			trt <- A[, ipch[1]]
			for(i in 2:nipch){
				name.t <- paste(name.t,names(A)[ipch][i],sep=":")
				trt <- paste(trt,A[,ipch[i]],sep=":")
			}}
		name.y <- names(A)[1]
	}
	junto <- subset(data.frame(y, trt), is.na(y) == FALSE)
	Mean<-mean(junto[,1])
	CV<-sqrt(MSerror)*100/Mean
	means <- tapply.stat(junto[, 1], junto[, 2], stat = "mean")
	sds <- tapply.stat(junto[, 1], junto[, 2], stat = "sd")
	nn <- tapply.stat(junto[, 1], junto[, 2], stat = "length")
	mi<-tapply.stat(junto[,1],junto[,2],stat="min") # change
	ma<-tapply.stat(junto[,1],junto[,2],stat="max") # change
	std.err <- sqrt(MSerror)/sqrt(nn[, 2]) # change sds[,2]
	Tprob <- qt(1 - alpha/2, DFerror)
	LCL <- means[, 2] - Tprob * std.err
	UCL <- means[, 2] + Tprob * std.err
	means <- data.frame(means, std=sds[,2], r = nn[, 2],
	LCL, UCL,Min=mi[,2],Max=ma[,2])
	names(means)[1:2] <- c(name.t, name.y)
	ntr <- nrow(means)
	nk <- choose(ntr, 2)
	if (p.adj != "none") {
		a <- 1e-06
		b <- 1
		for (i in 1:100) {
			x <- (b + a)/2
			xr <- rep(x, nk)
			d <- p.adjust(xr, p.adj)[1] - alpha
			ar <- rep(a, nk)
			fa <- p.adjust(ar, p.adj)[1] - alpha
			if (d * fa < 0)
				b <- x
			if (d * fa > 0)
				a <- x
		}
		Tprob <- qt(1 - x/2, DFerror)
	}
	nr <- unique(nn[, 2])
	if(console){
	cat("\nStudy:", main)
	if(console)cat("\n\nLSD t Test for", name.y, "\n")
	if (p.adj != "none")cat("P value adjustment method:", p.adj, "\n")
	cat("\nMean Square Error: ", MSerror, "\n\n")
	cat(paste(name.t, ",", sep = ""), " means and individual (",
			(1 - alpha) * 100, "%) CI\n\n")
	print(data.frame(row.names = means[, 1], means[, -1]))
	cat("\nalpha:", alpha, "; Df Error:", DFerror)
	cat("\nCritical Value of t:", Tprob, "\n")
	}
	if (group) {
		if (length(nr) == 1) {
			LSD <- Tprob * sqrt(2 * MSerror/nr)
			if(console)cat("\nLeast Significant Difference", LSD)
    statistics<-data.frame(Mean=Mean,CV=CV,MSerror=MSerror,LSD=LSD)
		}
		else {
		if(console)cat("\nMinimum difference changes for each comparison\n")
		#	nr1 <- 1/mean(1/nn[, 2])
		#	LSD <- Tprob * sqrt(2 * MSerror/nr1)
		#	cat("\nLeast Significant Difference", LSD)
		#	cat("\nHarmonic Mean of Cell Sizes ", nr1)
    statistics<-data.frame(Mean=Mean,CV=CV,MSerror=MSerror)
		}
		if(console){cat("\nMeans with the same letter are not significantly different.")
		cat("\n\nGroups, Treatments and means\n")}
		groups <- order.group(means[, 1], means[, 2], means[,
						4], MSerror, Tprob, means[, 3],alpha=alpha,console=console)
		w <- order(means[, 2], decreasing = TRUE)
		groups <- data.frame(groups[,1:3])
		comparison=NULL
	}
	if (!group) {
		LSD=" "
		comb <- utils::combn(ntr, 2)
		nn <- ncol(comb)
		dif <- rep(0, nn)
		pvalue <- dif
		sdtdif <- dif
		sig <- rep(" ", nn)
		for (k in 1:nn) {
			i <- comb[1, k]
			j <- comb[2, k]
#            if (means[i, 2] < means[j, 2]) {
#                comb[1, k] <- j
#                comb[2, k] <- i
#            }
			dif[k] <-means[i, 2] - means[j, 2]
			sdtdif[k] <- sqrt(MSerror * (1/means[i, 4] + 1/means[j,
										4]))
			pvalue[k] <- 2 * (1 - pt(abs(dif[k])/sdtdif[k], DFerror))
		}
		if (p.adj != "none")
			pvalue <- p.adjust(pvalue, p.adj)
		pvalue <- round(pvalue,4)
		LCL1 <- dif - Tprob * sdtdif
		UCL1 <- dif + Tprob * sdtdif
		for (k in 1:nn) {
			if (pvalue[k] <= 0.001)
				sig[k] <- "***"
			else if (pvalue[k] <= 0.01)
				sig[k] <- "**"
			else if (pvalue[k] <= 0.05)
				sig[k] <- "*"
			else if (pvalue[k] <= 0.1)
				sig[k] <- "."
		}
		tr.i <- means[comb[1, ], 1]
		tr.j <- means[comb[2, ], 1]
		comparison <- data.frame(Difference = dif, pvalue = pvalue,
				"sig."=sig, LCL = LCL1, UCL = UCL1)
		rownames(comparison) <- paste(tr.i, tr.j, sep = " - ")
		if(console){cat("\nComparison between treatments means\n\n")
		print(comparison)}
		groups <- NULL
		statistics<-data.frame(Mean=Mean,CV=CV,MSerror=MSerror)
	}
		test<-"Fisher-LSD"
		if(p.adj!="none")test<-p.adj
		parameters<-data.frame(Df=DFerror,ntr = ntr, t.value=Tprob,alpha=alpha,test=test,name.t=name.t)
		if(p.adj!="none") names(parameters)[3]<-p.adj
		rownames(parameters)<-" "
		rownames(statistics)<-" "
		rownames(means)<-means[,1]
		means<-means[,-1]
	output<-list(statistics=statistics,parameters=parameters, 
	means=means,comparison=comparison,groups=groups)
	invisible(output)
}

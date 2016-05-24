rowEpistatic <- function(mat.snp, cl, genes=NULL, warnError=TRUE){
	if(!is.matrix(mat.snp))
		stop("mat.snp must be a matrix.")
	colEpistatic(t(mat.snp), cl, genes=genes, warnError=warnError)
}

colEpistatic <- function(mat.snp, cl, genes=NULL, warnError=TRUE){
	if(!is.matrix(mat.snp))
		stop("mat.snp has to be a matrix.")
	if(any(!mat.snp %in% c(0:2, NA)))
		stop("The values in mat.snp must be 0, 1, and 2.")
	if(any(is.na(cl)))
		stop("No missing values allowed in cl.")
	if(any(!cl %in% (0:1)))
		stop("cl must consist of 0's (coding for controls) and 1's (coding for cases).")
	if(nrow(mat.snp) != length(cl))
		stop("The number of rows in mat.snp (representing subjects) must be equal to\n",
			"to the length of cl.")
	x <- mat.snp - 1
	z <- (x==0) - 0.5
	n.snp <- ncol(mat.snp)
	if(is.null(genes))
		combs <- allCombs(n.snp)
	else{
		if(!is.character(genes))
			stop("genes must be a vector of character strings.")
		if(length(genes) != n.snp)
			stop("The length of genes must be equal to the number of SNPs,\n",
				"i.e. the number of columns in mat.snp.")
			ids.genes <- as.numeric(as.factor(genes))
			combs <- allBetweenCombs(ids.genes)
	}
	n.combs <- nrow(combs)
	ll.main <- ll.full <- numeric(n.combs)
	vec.error <- vector("list", n.combs)
	if(warnError){
		wa <- options()$warn
		on.exit(options(warn=wa))
		options(warn = 2)
	}
	for(i in 1:n.combs){
		x1 <- x[,combs[i,1]]
		z1 <- z[,combs[i,1]]
		x2 <- x[,combs[i,2]]
		z2 <- z[,combs[i,2]]
		woIA <- try(glm(cl ~ x1 + z1 + x2 + z2, family="binomial"), silent=TRUE)
		if(is(woIA, "try-error")){
			ll.main[i] <- NA
			vec.error[[i]] <- woIA
			cat("NOTE: A warning has occurred when testing the SNPs in columns ", combs[i,1],
				" and ", combs[i,2], " of mat.snp.\n", "Therefore, the likelihood and ",
				"the test statistic is set to NA.\n", "(For this warning, see vec.error[[",
				 i, "]] in the output of colCordell.)\n\n", sep="")
		}
		else
			ll.main[i] <- woIA$deviance/-2
		full <- try(glm(cl ~ x1 + z1 + x2 + z2 + x1*x2 + x1*z2 + z1*x2 + z1*z2,
			family="binomial"), silent=TRUE)
		if(is(full, "try-error")){
			ll.full[i] <- NA
			vec.error[[i]] <- full
			cat("NOTE: A warning has occurred when testing the SNPs in columns ", combs[i,1],
				" and ", combs[i,2], " of mat.snp.\n", "Therefore, the likelihood and ",
				"the test statistic are set to NA.\n", "(For this warning, see vec.error[[",
				 i, "]] in the output of colCordell.)\n\n", sep="")
		}
		else
			ll.full[i] <- full$deviance/-2
	}
	if(warnError)
		options(warn=wa)
	stat <- -2 * (ll.main - ll.full)
	pval <- pchisq(stat, 4, lower.tail=FALSE)
	snpnames <- colnames(mat.snp)
	if(is.null(snpnames))
		snpnames <- paste("SNP", 1:n.snp, sep="")
	names(ll.full) <- names(stat) <- names(pval) <- names(ll.main) <- 
		names(vec.error) <- paste(snpnames[combs[,1]], snpnames[combs[,2]], sep=" : ")
	if(!is.null(genes))
		genecombs <- paste(genes[combs[,1]], genes[combs[,2]], sep=" : ")
	else
		genecombs <- NULL
	out <- list(ll.main=ll.main, ll.full=ll.full, stat=stat, pval=pval, genes=genecombs, vec.error=vec.error)
	class(out) <- "colEpi"
	out
}

print.colEpi <- function(x, top=5, digits=4, ...){
	pval <- format.pval(x$pval, digits=digits)
	out <- data.frame("LL (with IAs)" = x$ll.full, "LL (w/o IAs)" =x$ll.main,
		Statistic=x$stat, "P-Value"=pval, check.names=FALSE, stringsAsFactors=FALSE)
	if(!is.null(x$genes))
		out <- data.frame(out, Genes=x$genes, check.names=FALSE, stringsAsFactors=FALSE)
	cat("      Cordell's Likelihood Ratio Test for Epistatic Interactions", "\n\n")
	if(top>0 && length(x$ll.main)>top){
		ord <- order(x$pval)[1:top]
		out <- out[ord,]
		cat("Top", top, "SNP Interactions:\n")
	}
	print(format(out,digits=digits))
}

allCombs <- function (m, upperID = FALSE){
	vec1 <- rep.int(1:(m - 1), (m - 1):1)
    	txt <- paste("c(", paste(2:m, "m", sep = ":", collapse = ","), ")", sep = "")
    	vec2 <- eval(parse(text = txt))
    	if (upperID) 
        	return(vec1 + m * (vec2 - 1))
    	cbind(vec1, vec2)
}

allBetweenCombs <- function (gene){
	ids <- 1:length(gene)
    	mat.comb <- NULL
    	for (i in 1:(max(gene) - 1)){
        	tmp1 <- ids[gene == i]
        	tmp2 <- ids[gene > i]
        	tmpmat <- cbind(rep(tmp1, e = length(tmp2)), rep.int(tmp2, length(tmp1)))
        	mat.comb <- rbind(mat.comb, tmpmat)
    	}
    	mat.comb
}


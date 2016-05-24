# converts results of hwx.test into a list of hwtest objects
#' Convert results of \code{\link{hwx.test}} to a single list of \code{hwtest} objects.
#' 
#' There are two main uses of \code{listify}. You can simplify a complex result from \code{\link{hwx.test}} containing multiple populations and multiple loci into a simple list of \code{hwtest} objects. At the same time, you have a chance to change the parameters \code{detail} and \code{statName}. Useful to get output from a test.
#' 
#' @param hwlist the results of a call to \code{\link{hwx.test}}. It can be an \code{hwtest} object, a list of them or a list of lists of them.
#' @param detail Used only if you wish to reset the \code{detail} of each object.
#' @param statName Used only if you want to rest the \code{statName} of each object
#' 
#' 
#' @return a list of \code{hwtest} objects, possibly with their \code{detail} and \code{statName} parameters reset
#' 
#' @examples
#' data(HWcases)
#' outcome <- hwx.test(HWcases, detail=4, statName="LLR")
#' listify(outcome, detail=1, statName="U")

#' @export
listify <- function(hwlist, detail = NA, statName = NA) {
	# If it's a plain old hwtest, make it a list
	if (class(hwlist) == "hwtest") 
		hwlist <- list(result = hwlist)
	# It should now be a list
	if (class(hwlist) != "list") 
		stop("Format must be list of hwtest")
	# If it is a list of lists, flatten it by 1
	if (all(lapply(hwlist, class) == "list", na.rm = T)) 
		hwlist <- unlist(hwlist, recursive = F)
	# All elements should now be hwtest
	if (any(lapply(hwlist, class) != "hwtest", na.rm = T)) 
		stop("Format must be a list of hwtest")
	if (!is.na(detail)) 
		hwlist <- lapply(hwlist, function(z) {
			z$detail <- detail
			z
		})
	if (!is.na(statName)) {
		statNames <- c("LLR", "Prob", "U", "Chisq");
		statID <- which(statNames==statName);
		hwlist <- lapply(hwlist, function(z) {
			z$statName <- statName
			z$p.value <- z$Pvalues[[statID]]
			z
		})
	}
	isnum <- sapply(hwlist, function(x) is.character(x$method))
	hwlist <- hwlist[isnum]
	class(hwlist) <- "hwlist"
	hwlist
}


#' Construct a data frame from \code{\link{hwx.test}} output
#' 
#' If the \code{\link{hwx.test}} output has multiple populations and/or multiple loci, use this function to make a data frame to display the results in tabular form. 

#' 
#' @param hwlist The output from a call to \code{\link{hwx.test}}
#' @param statName gives you the option of changing which statistic's P value is reported
#' @param showN whether to show a column of sample size (number of diploids in the sample)
#' @param showk whether to show the number of alleles
#' @param showMethod whether to show whether the exact or Monte Carlo method was used
#' @param showSE whether to include the standard error for those tests which used the Monte Carlo method
#' @param showTables whether to show the total number of tables examined when full enumeration (exact) method is used
#' @param showTrials whether to show the number of random trials when Monte Carlo method is used
#' @param showStat whether to show the observed statistic
#' @param showAsymptoticX2 whether to include the asymptotic P value corresponding to the Pearson \eqn{X^2} statistic
#' @param showAsymptoticG2 whether to include the asymptotic P value for the \code{LLR} statistic

#' @export
hwdf <- 
function(hwlist, statName=NA, showN=TRUE, showk=TRUE, showMethod=TRUE, showSE=TRUE, showTables=TRUE, showTrials=TRUE, showStat=TRUE, showAsymptoticX2=FALSE, showAsymptoticG2=FALSE) {
	hwlist <- listify(hwlist, statName=statName)
	statNames <- c("LLR", "Prob", "U", "Chisq")
	statID <- which(statNames==hwlist[[1]]$statName)
	P.Value <- unlist(sapply(hwlist, function(x) x$Pvalues[[statID]]))
	method <- sapply(hwlist, function(x) x$method)
	exact <- method=="exact"
	monte <- method=="monte"
	ntrials <- sapply(hwlist, function(x) x$ntrials)
	ntrials <- unlist(ifelse(monte, ntrials, NA))
	observedStat <- unlist(sapply(hwlist, function(x) x$observed[[statID]]))
	tableCount <- sapply(hwlist, function(x) x$tableCount)
	N <- sapply(hwlist, function(x) sum(x$alleles)/2)
	k <- sapply(hwlist, function(x) length(x$alleles))
	denom <- ifelse(monte, ntrials, 1)
	denom <- unlist(denom)
	se.all <- sqrt(abs(P.Value * (1-P.Value)/denom))
	SE <- ifelse(monte, se.all, NA)
	pm <- ifelse(monte, "\u00b1","")
	columns <- list(P.Value)
	names <- paste("P-val(", statNames[[statID]],")", sep="")
	if(showSE && any(monte)){
		names <- c(names," ", "\u00b1 SE")
		columns$pm <- pm
		columns$SE <- SE
	}
	if(showStat) {
		names <- c(names, paste("obs-", statNames[[statID]], sep=""))
		columns$observedStat <- observedStat
	}
	if(showMethod){
		names <- c(names, "Method")
		columns$method <- unlist(method)
	}
	if(showTables && any(method=="exact")){
		names <- c(names, "Tables")
		columns$Tables <- unlist(ifelse(exact, tableCount, NA))
	}
	if(showTrials && any(method=="monte")){
		names <- c(names, "Trials")
		columns$Trials <- ntrials
	}
	if(showN){
		names <- c(names, "N")
		columns$N <- N
	}
	if(showk) {
		names <- c(names, "k")
		columns$k <- k
	}
	if(showAsymptoticX2){
		x2 <- unlist(sapply(hwlist, function(x) x$observed[["Chisq"]]))
		df <- k * (k-1)/2;
		asyX2 <- pchisq(x2, df, lower.tail=F)
		columns$asyX2 <- asyX2
		names <- c(names, "AsymP(X2)")
	}
	if(showAsymptoticG2){
		g2 <- unlist(sapply(hwlist, function(x) (-2) * x$observed[["LLR"]]))
		df <- k * (k-1)/2;
		asyG2 <- pchisq(g2, df, lower.tail=F)
		columns$asyG2 <- asyG2
		names <- c(names, "AsymP(G2)")
	}
	df <- as.data.frame(columns)
	names(df) <- names
	names(columns) <- names
	df
	
	
}
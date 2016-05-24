"npsubtreated" <- function(trait,medication,increase=FALSE) {
	if (length(trait) != length(medication)) stop("trait and medication vector lengths differ")
	df <- data.frame(y.obs=trait,med=medication,y.imp=rep(NA,length(trait)))
	rownames(df) <- paste("c",c(1:length(trait)),sep="")
	meaids <- !(is.na(trait) | is.na(medication))
	if (sum(meaids)<=0) stop("no non-missing values supplied")
	x1 <- df[meaids,]
	levs <- levels(as.factor(x1$med))
	if (length(levs)>2) stop("more than 2 levels in medication variable")
	if (length(levs)<2) stop("less than 2 levels in medication variable")
	if (levs[1] != "0" && levs[2]!="1") stop("medication should be coded with 0 and 1")
	
	if (increase) x1$y.obs <- (-1)*x1$y.obs

	ord <- sort.int(x1$y.obs,decreasing=T,index.return=T)$ix
	x1 <- x1[ord,]
	subtr <- mean(x1$y.obs)
	x1$y.1 <- x1$y.obs - subtr
	x1$y.2 <- rep(NA,dim(x1)[1])
	for (k in c(1:dim(x1)[1])) {
		if (x1$med[k]) {
			if (k>1) 
				x1$y.2[k] <- (sum(x1$y.2[1:(k-1)])+x1$y.1[k])/k
			else 
				x1$y.2[k] <- x1$y.1[k]
		} else {
			x1$y.2[k] <- x1$y.1[k]
		}
	}
	x1$y.imp <- x1$y.2 + subtr
	df[rownames(x1),"y.imp"] <- x1$y.imp

	if (increase) df$y.imp <- (-1)*df$y.imp

	return(df$y.imp)
}

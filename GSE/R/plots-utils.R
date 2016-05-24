.qqplot.pmdadj <- function(object, cutoff=0.99, xlog10=FALSE, ylog10=FALSE){
	p <- object$p
	pmd.adj <- object$pmd.adj
	id <- 1:length(pmd.adj)

	## identify any rows with completely missing data
	id.na <- which( is.na(pmd.adj))
	pmd.na <- pmd.adj[ id.na]
	pmd.adj <- pmd.adj[ setdiff(id, id.na) ] 
	id <- id[ setdiff(id, id.na) ] 

	## plot
	n <- length(pmd.adj)
	id.order <- order(pmd.adj)
	pmd.adj <- pmd.adj[id.order]
	id <- id[id.order]
	pmd.exp <- qchisq((1:n - 0.5)/n, df=p)
	threshold <- qchisq( cutoff, df=p)
	outliers <- factor(ifelse(pmd.adj > threshold, 2, 1), levels=c(1,2), labels=c("N","Y"))

	## output
	x <- pmd.exp
	y <- pmd.adj
	myplot <- ggplot( data=data.frame(x=x, y=y, outliers=outliers), aes(x=x, y=y)) + 
					geom_point(aes(color=outliers)) + 
					geom_hline(aes(yintercept=threshold), col="grey50", lty=2) + 
					xlab("expected") + ylab("observed") + ggtitle("QQ plot")+ theme_bw()
	if( xlog10 ) myplot <- myplot + scale_x_log10()
	if( ylog10 ) myplot <- myplot + scale_y_log10()
	myplot
}

.distplot.pmdadj <- function(object, cutoff=0.99, ylog10=FALSE){
	p <- object$p
	pmd.adj <- object$pmd.adj
	id <- 1:length(pmd.adj)

	## identify any rows with completely missing data
	id.na <- which( is.na(pmd.adj))
	pmd.na <- pmd.adj[ id.na]
	pmd.adj <- pmd.adj[ setdiff(id, id.na) ] 
	id <- id[ setdiff(id, id.na) ] 

	## plot
	threshold <- qchisq( cutoff, df=p)
	outliers <- factor(ifelse(pmd.adj > threshold, 2, 1), levels=c(1,2), labels=c("N","Y"))
	x <- id
	y <- pmd.adj
	myplot <- ggplot( data=data.frame(x=x, y=y, outliers=outliers), aes(x=x, y=y)) + 
					geom_point(aes(color=outliers)) + 
					geom_hline(aes(yintercept=threshold), col="grey50", lty=2) + 
					xlab("case number") + ylab("observed adjusted squared distances") + ggtitle("Index plot")+ theme_bw()
	if( ylog10 ) myplot <- myplot + scale_y_log10()
	myplot
}

.ddplot.pmdadj <- function(object, cutoff=0.99, xlog10=FALSE, ylog10=FALSE){
	p <- object$p
	pmd.adj <- object$pmd.adj
	id <- 1:length(pmd.adj)

	## obtain classical distances
	pmd.adj.EM <- CovEM(object$x)
	pmd.adj.EM <- pmd.adj.EM@pmd.adj

	## identify any rows with completely missing data
	id.na <- which( is.na(pmd.adj))
	pmd.na <- pmd.adj[ id.na]
	pmd.adj <- pmd.adj[ setdiff(id, id.na) ] 
	pmd.adj.EM <- pmd.adj.EM[ setdiff(id, id.na)]
	id <- id[ setdiff(id, id.na) ] 

	## plot
	threshold <- qchisq( cutoff, df=p)
	outliers <- factor(ifelse(pmd.adj > threshold | pmd.adj.EM > threshold, 2, 1), levels=c(1,2), labels=c("N","Y"))
	
	## main plot
	x <- pmd.adj.EM
	y <- pmd.adj
	myplot <- ggplot( data=data.frame(x=x, y=y, outliers=outliers), aes(x=x, y=y)) + 
				geom_point(aes(color=outliers)) + 
				geom_hline(aes(yintercept=threshold), col="grey50", lty=2) + 
				geom_vline(aes(xintercept=threshold), col="grey50", lty=2) + 
				xlab("classical distances") + ylab("robust distances") + ggtitle("Distance-distance plot") + theme_bw()
	if( xlog10 ) myplot <- myplot + scale_x_log10()
	if( ylog10 ) myplot <- myplot + scale_y_log10()
	myplot
}

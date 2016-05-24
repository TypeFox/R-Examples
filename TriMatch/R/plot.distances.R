utils::globalVariables(c('ID','Distance','Model','id'))

#' Barplot for the sum of distances.
#' 
#' @param tmatch the results of \code{\link{trimatch}}.
#' @param caliper a vector indicating where vertical lines should be drawn as a
#'        factor of the standard deviation. Rosenbaum and Rubin (1985) suggested
#'        one quarter of one standard deviation.
#' @param label label the bars that exceed the minimum caliper.
#' @seealso triangle.match
#' @export
distances.plot <- function(tmatch, caliper=.25, label=FALSE) {
	tpsa <- attr(tmatch, 'triangle.psa')
	tmatch$id <- row.names(tmatch)
	tmp <- melt(tmatch[,c('id','D.m1','D.m2','D.m3')], id.vars=1)
	names(tmp) <- c('ID','Model','Distance')
	l <- (sd(tpsa$ps1, na.rm=TRUE) + sd(tpsa$ps2, na.rm=TRUE) + sd(tpsa$ps3, na.rm=TRUE))
	stdev <- sapply(tpsa[,c('ps1','ps2','ps3')], sd, na.rm=TRUE)
	cat(paste("Standard deviations of propensity scores: ",
			  paste(prettyNum(stdev, digits=2), collapse=' + '),
			  " = ", prettyNum(l, digits=2), '\n', sep='' ))
	#sddf <- data.frame(Model=c('1','2','3'), stdev=stdev)
	#sddf$cumsum <- cumsum(sddf$stdev)
	sddf <- data.frame(Model=c('1','2','3'), cumsum=caliper * 1:3)
	levels(tmp$Model) <- levels(sddf$Model)
	p <- ggplot(tmp, aes(x=ID, y=Distance)) +
		geom_bar(aes(fill=Model), stat='identity', width=1) + coord_flip() + 
		#geom_hline(yintercept=3 * caliper, color='black') +
		geom_hline(data=sddf, aes(yintercept=cumsum, color=Model), linetype=2) +
		theme(axis.text.y=element_blank()) +
		xlab(NULL) + xlim(tmatch[order(tmatch$Dtotal),'id'])
	if(label & length(which(tmatch$Dtotal > 3 * (min(caliper)))) > 0) {
		p <- p + geom_text(data=tmatch[tmatch$Dtotal > 3 * (min(caliper)),], 
						   aes(x=id, y=0, label=id),
						   size=3, hjust=1)	
	}
	for(i in seq_along(caliper)) {
		cat(paste("Percentage of matches exceding a distance of ", 
				  prettyNum(caliper[i], digits=2), 
				  " (caliper = ", caliper[i], "): ", 
				  prettyNum(length(which(tmatch$Dtotal > 3 * caliper[i])) / 
				  		  	length(tmatch$Dtotal) * 100, digits=2), '%\n',
				  sep=''))
	}
	return(p)
}

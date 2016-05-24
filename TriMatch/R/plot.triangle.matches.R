utils::globalVariables(c('y','id'))

#' Triangle plot drawing matched triplets.
#' 
#' This plot function adds a layer to \code{\link{plot.triangle.psa}} drawing 
#' matched triplets. If \code{p} is supplied, this function will simply draw on
#' top of the pre-existing plot, otherwise \code{\link{plot.triangle.psa}} will
#' be called first.
#' 
#' If this function calls \code{\link{plot.triangle.psa}}, it will only draw
#' line segements and points for those data rows that were used in the matching
#' procedure. That is, data elements not matched will be excluded from the
#' figure. To plot all segments and points regardless if used in matching, set
#' \code{p = plot(tpsa)}.
#' 
#' @param x matched triplets from \code{link{triangle.match}}.
#' @param sample an number between 0 and 1 representing the percentage of matched
#'        triplets to draw.
#' @param rows an integer vector corresponding to the rows in \code{tmatch} to draw.
#' @param line.color the line color.
#' @param line.alpha the alpha for the lines.
#' @param point.color color of matched triplet points.
#' @param point.size point size for matched triplets.
#' @param p a ggplot to add the match lines. If NULL, then \code{\link{plot.triangle.psa}}.
#' @param ... other parameters passed to \code{\link{plot.triangle.psa}}.
#' @return a \code{ggplot2} graphic.
#' @seealso plot.triangle.psa
#' @seealso triangle.match
#' @method plot triangle.matches
#' @export
plot.triangle.matches <- function(x, 
								  sample=.05, 
								  rows=sample(nrow(tmatch), nrow(tmatch) * sample), 
								  line.color='black',
								  line.alpha=.5,
								  point.color='black',
								  point.size=3,
								  p, ...) {
	tmatch <- x
	tpsa <- attr(tmatch, 'triangle.psa')
	pts.overlay <- data.frame(x=numeric(), y=numeric(), id=integer())
	for(i in rows) {
		tmp <- tpsa[as.integer(tmatch[i,1:3]),]
		tmp$id <- i
		tmp <- melt(tmp[,c('id','ps1','ps2','ps3')], id.vars='id')
		tmp <- tmp[!is.na(tmp$value),]
		tmp <- rbind(
			as.data.frame(matrix(unlist(lapply(tmp[which(tmp$variable == 'ps1'),]$value, segment1)), 
								 ncol=2, byrow=TRUE)),
			cbind(tmp[which(tmp$variable == 'ps2'),]$value, rep(0, 2)),
			as.data.frame(matrix(unlist(lapply(tmp[which(tmp$variable == 'ps3'),]$value, segment2)), 
								 ncol=2, byrow=TRUE))
		)
		tmp <- tmp[c(1:4,6,5),]
		tmp <- rbind(tmp, tmp[1,])
		tmp$id <- i
		pts.overlay <- rbind(pts.overlay, tmp)
	}
	names(pts.overlay) <- c('x','y','id')
	if(missing(p)) {
		if(length(which(names(as.list(match.call(expand.dots=TRUE))) == 'draw.segments')) == 0) {
			p <- plot(tpsa, draw.segments=FALSE, ...)
		} else {
			tpsa2 <- tpsa[which(tpsa$id %in% unlist(c(tmatch[,1:3]))),]
			p <- plot(tpsa2, ...)
		}
	}
	p <- p + geom_path(data=pts.overlay, aes(x=x, y=y, group=id), 
					   alpha=line.alpha, colour=line.color) +
		geom_point(data=pts.overlay, aes(x=x, y=y), colour=point.color, size=point.size)
	return(p)
}	

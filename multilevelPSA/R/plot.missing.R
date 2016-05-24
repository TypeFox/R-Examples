utils::globalVariables(c('y','X2','X1','value'))

#' Returns a heat map graphic representing missinging of variables grouped by
#' the given grouping vector.
#'
#' NOTE: This is an experimental function and the results may vary depending on the
#'       nature of the dataset.
#' 
#' @param x a data frame containing the variables to visualize missingness
#' @param grouping a vector of length nrow(vars) corresponding to how missing will be grouped by
#' @param grid whether to draw a grid between tiles
#' @param widths the ratio of the widths of the heatmap and histogram.
#' @param heights the ratio of the heights of the heatmap and histogram.
#' @param color the color used for indicating missingness.
#' @param ... currently unused.
#' @return a ggplot2 expression
#' @seealso plot.mlpsa
#' @export
missing.plot <- function(x, grouping, grid=FALSE, 
						 widths=c(ggplot2::unit(3, 'null'), ggplot2::unit(1,'inches')), 
						 heights=c(ggplot2::unit(1,'inches'), ggplot2::unit(3, 'null')),
						 color='red',
						 ...) {
	vars = x
	empty <- plyr::empty
		
	colMissing = apply(vars, 2, function(x) sum(is.na(x))) / nrow(vars)
	colMissing = 100 * colMissing
	colMissing = data.frame(x=names(colMissing), y=as.numeric(colMissing))
	
	phist.right = ggplot(colMissing, aes(x=x, y=y, fill=y)) + 
						geom_bar() + coord_flip()
	phist.right = phist.right + xlab(NULL) + ylab(NULL)
	phist.right = phist.right + scale_fill_gradient('Missingness', low='white', high=color, 
						limits=c(0,100), breaks=seq(0, 100, 10), 
						labels=paste(seq(0,100,10), '%', sep=''))
	phist.right = phist.right + geom_text(aes(label=round(y, digits=0)), size=2)
	phist.right = phist.right + ylim(c(0,100))
	phist.right = phist.right + theme(legend.position='none', 
						axis.text.y=element_text(size=6, angle=0, hjust=.5, vjust=.5),
						axis.text.x=element_blank(), 
						axis.ticks=element_blank())
	
	for(g in unique(grouping)) {
		tmp = apply(vars[grouping == g,], 2, function(x) sum(is.na(x))) / 
			nrow(vars[grouping == g,])
		tmp = 100 * tmp
		tmp = as.numeric(tmp)
		colMissing = cbind(colMissing, tmp)
		names(colMissing)[ncol(colMissing)] = g
	}
	rowMissing = apply(colMissing[,3:ncol(colMissing),], 2, mean)
	rowMissing = data.frame(x=names(rowMissing), y=as.numeric(rowMissing))
	phist.top = ggplot(rowMissing, aes(x=x, y=y, fill=y)) + 
					geom_bar() + ylim(c(0,100))
	phist.top = phist.top + scale_fill_gradient('Missingness', low='white', high=color, 
					limits=c(0,100), breaks=seq(0, 100, 10), 
					labels=paste(seq(0,100,10), '%', sep=''))
	phist.top = phist.top + xlab(NULL) + ylab(NULL)
	phist.top = phist.top + theme(axis.text.x=element_text(size=6, angle=-90, hjust=.5, vjust=.5), 
					axis.text.x=element_blank(), axis.ticks=element_blank())
	phist.top = phist.top + geom_text(aes(label=round(y, digits=0)), size=2)
	phist.top = phist.top + theme(legend.position='none', axis.text.y=element_blank())
	
	grouping = as.character(grouping)
	grouping[is.na(grouping)] = 'Unknown'
	grouping = as.factor(grouping)
	testing.NA = matrix(ncol=length(levels(grouping)), nrow=(ncol(vars)))
	for(i in 1:(dim(vars)[2])) {
		testing.NA[i,] = tapply(vars[[i]], grouping, function(x) sum(is.na(x)) / length(x))
	}
	testing.NA = testing.NA * 100
	dimnames(testing.NA) = list(names(vars), unique(grouping))
	testing.NA2 = melt(testing.NA)
	p = ggplot(testing.NA2, aes(x=X2, y=X1, fill=value))
	if(grid) {
		p = p + geom_tile(colour='grey')
	} else {
		p = p + geom_tile()
	}
	p = p + xlab(NULL) + ylab(NULL)
	p = p + theme(axis.ticks=element_blank(), 
				 axis.text.y=element_text(size=6, hjust=1, vjust=.5), 
				 axis.text.x=element_text(size=6, angle=-90, hjust=0, vjust=.5))
	p = p + scale_fill_gradient('Missingness', low='white', high=color, 
				 limits=c(0,100), breaks=seq(0, 100, 10), 
				 labels=paste(seq(0,100,10), '%', sep=''))
	p = p + geom_text(aes(label=round(value, digits=0)), size=2, colour='black')
	p = p + theme(legend.position='none', 
				  axis.text.x=element_blank(), 
				  axis.text.y=element_blank())

	grid_layout = grid.layout(nrow=2, ncol=2, widths=widths, heights=heights)
 	grid.newpage()
 	pushViewport( viewport( layout=grid_layout ) )
 	align.plots(grid_layout, list(p, 2, 1), 
 							list(phist.right, 2, 2), list(phist.top,1,1))
}

utils::globalVariables(c('Desc','level2.Freq','Var.Freq','value'))

#' Heat map representing variables used in a conditional inference tree across level 2 variables.
#' 
#' This figure provides a summary of the covariates used within each level two cluster
#' along with their relative importance. Covariates are listed on the y-axis and 
#' level two clusters along the x-axis. Cells that are shaded indicate that that
#' covariate was present in the conditional. The shade of the color represents
#' the highest level witin the tree that covariate appeared. That is, the darkest
#' color, or depth 1, corresponds to the covariate used at the root of the tree, or
#' the first split.
#' 
#' @param x the results of \code{\link{mlpsa.ctree}}
#' @param colNames the columns to include in the graphic
#' @param level2Col the name of the level 2 column.
#' @param colLabels column labels to use. This is a data frame with two columns, the
#'        first column should match the values in \code{trees} and the second column
#'        the description that will be used for labeling the variables.
#' @param color.low color for variables with greater relative importance as determined
#'        by occurring sooner in the tree (closer to the root split).
#' @param color.high color for variables with less relative importance as determined
#'        by occurring later in the tree (further from the root split).
#' @param color.na color for variables that do not occur in the tree.
#' @param ... currently unused.
#' @return a ggplot2 expression
#' @seealso plot.mlpsa
#' @export
#' @examples
#' \dontrun{
#' require(party)
#' data(pisana)
#' data(pisa.colnames)
#' data(pisa.psa.cols)
#' mlctree = mlpsa.ctree(pisana[,c('CNT','PUBPRIV',pisa.psa.cols)], formula=PUBPRIV ~ ., level2='CNT')
#' student.party = getStrata(mlctree, pisana, level2='CNT')
#' tree.plot(mlctree, level2Col=pisana$CNT)
#' }
tree.plot <- function(x, colNames, level2Col, colLabels=NULL, 
					  color.high="azure", color.low='steelblue', color.na='white',
					  ...) {
	trees = x
	
	if(missing(colNames)) {
		colNames = names(trees[[1]]@data@env$input)
	}
	ncol = length(colNames) + 1
	tree.df <- as.data.frame(matrix(nrow=length(unique(level2Col)), ncol=ncol))
	names(tree.df) = c('level2', colNames)
	tree.df[,2:ncol(tree.df)] = as.numeric(NA)
	tree.df[,1] = unique(level2Col)
	
	convertTree <- function(tree.df.row, node, depth, tree.party) {
		nodeVal = nodes(tree.party, node)[[1]]
		varName = as.character(nodeVal$psplit)[7]
		if(nodeVal$terminal) {
			#Do nothing
		} else {
			tree.df.row[varName] = min(depth, tree.df.row[varName][1,1], na.rm=TRUE)
			tree.df.row = convertTree(tree.df.row, nodes(tree.party, node)[[1]]$left[[1]], 
									  depth+1, tree.party)
			tree.df.row = convertTree(tree.df.row, nodes(tree.party, node)[[1]]$right[[1]], 
									  depth+1, tree.party)
		}
		tree.df.row
	}
	
	for(i in names(trees)) {
		tree.df[which(tree.df[,1] == i),] = convertTree(tree.df[which(tree.df[,1] == i),], 
														1, 1, trees[[i]])
	}
	
	tree.df.m <- melt(tree.df[which(tree.df$level2 %in% unique(level2Col)),], id='level2')
	tree.df.m$level2 = as.factor(as.character(tree.df.m$level2))
	descColName = 'variable'
	if(!is.null(colLabels)) {
		tree.df.m = merge(tree.df.m, colLabels, by.x='variable', 
						  by.y=names(colLabels)[1], all.x=TRUE)
		descColName = names(colLabels)[2]
	}
	
	#Labels
	value.freq = as.data.frame(table(tree.df.m[!is.na(tree.df.m$value),]$variable))
	level2.freq = as.data.frame(table(tree.df.m[!is.na(tree.df.m$value),]$level2))
	tree.df.m = merge(tree.df.m, level2.freq, by.x='level2', by.y='Var1', all.x=TRUE)
	tree.df.m = merge(tree.df.m, value.freq, by.x='variable', by.y='Var1', all.x=TRUE)
	names(tree.df.m)[(ncol(tree.df.m)-1):ncol(tree.df.m)] = c('level2.Freq', 'Var.Freq')
	tree.df.m$Desc = as.character(tree.df.m[,descColName])
	tree.df.m$level2 = as.character(tree.df.m$level2)
	
	p = ggplot(tree.df.m, aes(level2, Desc)) + 
			geom_text(data=tree.df.m[!duplicated(tree.df.m$level2),], 
					  aes(x=level2, y='', label=level2.Freq), size=2, angle=-90, hjust=.5) +
			geom_text(data=tree.df.m[!duplicated(tree.df.m$variable),], 
					  aes(y=Desc, x='', label=Var.Freq), size=2, hjust=.5) +
			geom_tile(aes(fill = value)) +
			scale_fill_gradient("Depth", high=color.high, low=color.low, na.value=color.na) + 
			theme(axis.text.y=element_text(size=6, hjust=0, vjust=.5), 
				 axis.text.x=element_text(size=6, angle=-90, hjust=0, vjust=.5), 
				 axis.ticks=element_blank()) + 
			xlab(NULL) + ylab(NULL)
	return(p)
}

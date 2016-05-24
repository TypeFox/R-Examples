#' @title get FlexTable from a ftable object
#'
#' @description Get a \code{\link{FlexTable}} object from 
#' a \code{\link{ftable}} object.
#' 
#' @param x object to get \code{FlexTable} from
#' @param text.properties default text formatting properties
#' @param padding.left cells paragraphs left padding - 0 or positive integer value.
#' @param padding.right cells paragraphs right padding - 0 or positive integer value.
#' @param ... further arguments, not used. 
#' @return a \code{\link{FlexTable}} object
#' @examples
#' #
#' @example examples/ftable.R
#' @export 
as.FlexTable.ftable = function( x, text.properties = textNormal(), 
	padding.left = 4, padding.right = 4, ... ){
	
	row.vars = attr( x, "row.vars" )
	col.vars = attr( x, "col.vars" )
	rows = rev( expand.grid( rev( row.vars ) ) )
	rows = as.matrix( rows )
	cols = t( as.matrix(rev( expand.grid( rev( col.vars ) ) ) ) )
	
	fftable = format( x, quote = FALSE, method =  "col.compact" )
	data = fftable[ (nrow(cols) + 2):nrow(fftable) , (ncol(rows) + 1):(ncol(fftable))]
	str = cbind( rows, data )
	
	ft = FlexTable( data = str, 
		header.columns = FALSE, 
		add.rownames = FALSE )
	
	for(i in seq_len(nrow(cols))){
		labs = c( names(col.vars[i]), cols[i,] )
		colspan = c(length( row.vars), rep( 1, ncol(data) ) )
		ft = addHeaderRow( ft, labs, colspan = colspan )
	}
	labs = c(dimnames(rows)[[2]], "" )
	ft = addHeaderRow( ft, labs, colspan = c( rep(1, ncol( rows) ), ncol(data) ) )
	
	ft[,1:ncol(rows)] = parRight(padding.right = padding.right, padding.left = padding.left)
	ft[,1:ncol(rows), to = "header"] = parRight(padding.right = padding.right, padding.left = padding.left)
	ft[, (ncol(rows) + 1):(ncol(fftable))] = parCenter(padding.right = padding.right, padding.left = padding.left)
	ft[, (ncol(rows) + 1):(ncol(fftable)), to = "header"] = parCenter(padding.right = padding.right, padding.left = padding.left)
	
	ft[to="header"] = chprop( text.properties, font.weight = "bold" )
	ft[ ] = text.properties
	ft
}

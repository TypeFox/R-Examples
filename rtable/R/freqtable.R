#' @title FlexTable Percentage Table
#'
#' @description Get a \code{\link{FlexTable}} two-way frequency table. 
#' Table contains cell frequencies, cell percentages of the total frequency, 
#' cell percentages of row frequencies and cell percentages of column frequencies.
#' 
#' @param x a \code{table} object
#' @param percent_digits the desired number of digits after 
#' the decimal point for percentages.
#' @param base100 wether to multiply percentages by 100.
#' @param percent whether to include cell percentages in the table.  
#' @param row_percent whether to include cell percentages of 
#' row frequencies in the table.  
#' @param col_percent whether to include cell percentages of column 
#' frequencies in the table.
#' @param label_na label to use for \code{missing} level.
#' @param label_level label to use for level column (only for one way table).
#' @param label_sum label to use for margins.
#' @param label_count label to use for \code{frequencies}.
#' @param label_percent label to use for \code{percentages}.
#' @param label_percent_col label to use for \code{column frequencies}.
#' @param label_percent_row label to use for \code{row frequencies}.
#' @return a \code{\link{FlexTable}} object
#' @examples
#' #
#' @example examples/freqtable.R
#' @export 
freqtable = function( x, percent_digits = 1, 
	label_na = "Missing", label_level = "Level", 
	label_sum = "Sum", 
	label_count = "#", label_percent = "%", 
	base100 = TRUE, percent = T, 
	row_percent = T, col_percent = T, 
	label_percent_col = "col %", 
	label_percent_row = "row %" ) {

	if( !is.table (x) ) stop("x must be a table.")
	
	if( length( dim( x ) ) > 2 ){
		stop("x cannot have more than 2 dimensions.")
	}
	if( length( dim( x ) ) == 1 ){
		ft = freqtable_univariate(x, label_na = label_na, label_level = label_level, 
			label_sum = label_sum, label_count = label_count, 
			label_percent = label_percent, percent_digits = percent_digits, 
			base100 = base100, percent = percent )
	} else if( length( dim( x ) ) == 2 ){
		ft = freqtable_bivariate(x, percent_digits = percent_digits, 
			base100 = base100, percent = percent, 
			row_percent = row_percent, col_percent = col_percent, 
			label_na = label_na, 
			label_sum = label_sum, label_count = label_count, 
			label_percent = label_percent, 
			label_percent_col = label_percent_col, 
			label_percent_row = label_percent_row)
	} 
	
	ft
}


freqtable_univariate = function( c.table, 
		label_na, label_level, 
		label_sum , 
		label_count, label_percent, 
		percent_digits = 1, 
		base100 = TRUE, percent = T ) {
	
	p.table = prop.table( c.table )
	row.vars = names( c.table )
	row.vars[is.na(row.vars)] = label_na
	
	cumsump.table = c( cumsum(p.table), NA )
	
	c.table = c( as.integer( c.table ), as.integer( sum( c.table) ) )
	p.table = c( as.double( p.table ), 1.00 )
	c.table = sprintf( paste0("%.0f"), c.table )
	p.table = sprintf( paste0("%.", percent_digits, "f"), p.table * ifelse(base100, 100, 1) )
	
	row.vars = c( row.vars, label_sum )
	row.vars = factor( x = row.vars, levels = row.vars )
	out = data.frame( levels = row.vars, freq = c.table, prop = p.table)
	names( out ) = c( label_level, label_count, label_percent )
	if( !percent )
		out = out[, c(T,T,F)]
	
	ft = FlexTable( data = out, header.columns = T )
	
	ft = setFlexTableBorders( ft, header = T, body = T, 
			inner.vertical = borderNone(),
			inner.horizontal = borderNone(),
			outer.horizontal = borderNone(), 
			outer.vertical = borderNone())	
	
	parProp = parProperties(padding.bottom = 0, padding.top = 0, 
			padding.left = 3, padding.right = 3, text.align = "center")
	
	ft[, 1] = chprop( parProp, text.align = "right" )
	ft[, 1, to = "header"] = chprop( parProp, text.align = "right" )
	ft[, 2:ncol(out)] = parProp
	ft[, 2:ncol(out), to = "header"] = parProp
	
	ft[ , 1, side = "right"] = borderSolid()
	ft[ , 1, side = "right", to = "header"] = borderSolid()
	ft[ , ncol(out), side = "right"] = borderSolid()
	ft[ , ncol(out), side = "right", to = "header"] = borderSolid()
	ft[ nrow(out), , side = "top"] = borderSolid()
	ft[ 1, , side = "top"] = borderSolid()
	
	ft[,1] = textBold()
	ft[to="header"] = textBold()
	
	ft
}




freqtable_bivariate = function( c.table, percent_digits = 1, 
		base100 = TRUE, percent = T, row_percent = T, col_percent = T, 
		label_na, label_sum, 
		label_count, label_percent, 
		label_percent_col, label_percent_row ) {
	
	p.table = prop.table( c.table )
	pcol.table = prop.table( c.table, margin = 2 )
	prow.table = prop.table( c.table, margin = 1 )
	
	nbcol = ncol( c.table )
	nbrow = nrow( c.table )
	.sum = sum( c.table )
	
	.colsum = colSums( c.table )
	pcol.colsum = .colsum / .sum
	.rowsum = rowSums( c.table )
	pcol.rowsum = .rowsum / .sum
	
	# drop NA from names
	dm = dimnames( c.table )
	dm[[1]][is.na(dm[[1]])] = label_na
	dm[[2]][is.na(dm[[2]])] = label_na
	# add sum level in names
	dm[[1]] = c( dm[[1]], label_sum )
	dm[[2]] = c( dm[[2]], label_sum )
	
	
	# add bottom margin
	c.table = rbind( as.matrix( c.table ), .colsum )
	# add right margin
	c.table = cbind( c.table, c( .rowsum, .sum ) )
	
	# add bottom margin
	p.table = rbind( as.matrix( p.table ), pcol.colsum )
	# add right margin
	p.table = cbind( p.table, c( pcol.rowsum, 1.00 ) )
	
	# add bottom margin
	pcol.table = rbind( as.matrix( pcol.table ), rep(1.00, nbcol ) )
	# add right margin
	pcol.table = cbind( pcol.table, c( pcol.rowsum, 1.00 ) )
	
	# add bottom margin
	prow.table = rbind( as.matrix( prow.table ), pcol.colsum )
	# add right margin
	prow.table = cbind( prow.table, rowSums( prow.table ) )
	
	c.table[,] = sprintf( "%.0f", c.table[,] )
	p.table[,] = sprintf( paste0("%.",percent_digits,"f"), p.table[,] * ifelse(base100, 100, 1) )
	pcol.table[,] = sprintf( paste0("%.",percent_digits,"f"), pcol.table[,] * ifelse(base100, 100, 1) )
	prow.table[,] = sprintf( paste0("%.",percent_digits,"f"), prow.table[,] * ifelse(base100, 100, 1) )
	c.table[c.table %in% c("NaN", "NA")] = ""
	p.table[p.table %in% c("NaN", "NA")] = ""
	pcol.table[pcol.table %in% c("NaN", "NA")] = ""
	prow.table[prow.table %in% c("NaN", "NA")] = ""
	
	dimnames( c.table ) = dm
	dimnames( p.table ) = dm
	dimnames( pcol.table ) = dm
	dimnames( prow.table ) = dm
	
	levels = rep( dm[[1]], sum(c( T, percent, row_percent, col_percent )) )
	levels = factor( x = levels, levels = dm[[1]] )
	labs_stat = c( label_count, label_percent, label_percent_row, label_percent_col )[c( T, percent, row_percent, col_percent )]
	stat.names = rep( labs_stat, each = dim( c.table )[1] )
	
	stat.names = factor( x = stat.names, levels = labs_stat )
	rawdata = c.table
	if( percent ) rawdata = rbind( rawdata, p.table )
	if( row_percent ) rawdata = rbind( rawdata, prow.table )
	if( col_percent ) rawdata = rbind( rawdata, pcol.table )
	rawdata = as.data.frame( rawdata, stringsAsFactors = F )
	row.names( rawdata ) = NULL 
	out = cbind( levels, stat.names, rawdata )
	
	out = out[ order( levels, stat.names ), ]
	
	if( row_percent ){
		out[out[,2]==label_percent_row, ncol(out)] = ""
		out = out[ -which(out[,2]==label_percent_row & out[,1]==label_sum), ]
	}
	if( col_percent ){
		out[out[,2] == label_percent_col, ncol(out)] = ""
		out = out[ -which(out[,2]==label_percent_col & out[,1]==label_sum), ]
	}
	ft = FlexTable( data = out, header.columns = FALSE, 
			add.rownames = FALSE  )
	
	ft = addHeaderRow( ft, c( "", names( out )[-c(1:2)] ), 
			colspan = c(2, rep( 1, ncol(out)-2 ) ) )
	
	ft = spanFlexTableRows(ft, j = 1, runs = as.character( out[,1]) )
	
	ft = setFlexTableBorders( ft, header = T, body = F, 
			inner.vertical = borderSolid(), 
			inner.horizontal = borderSolid(), 
			outer.horizontal = borderSolid(), 
			outer.vertical = borderSolid())	
	ft = setFlexTableBorders( ft, header = F, body = T, 
			inner.vertical = borderSolid(),
			inner.horizontal = borderNone(),
			outer.horizontal = borderSolid(), 
			outer.vertical = borderSolid())	
	
	parProp = parProperties(padding.bottom = 0, padding.top = 0, 
			padding.left = 3, padding.right = 3, text.align = "center")
	
	ft[, 1:2] = chprop( parProp, text.align = "right" )
	ft[, 1:2, to = "header"] = chprop( parProp, text.align = "right" )
	ft[, 3:ncol(out)] = parProp
	ft[, 3:ncol(out), to = "header"] = parProp
	
	line_markers = seq( from = 1, to = nrow(out) , by = sum( c( T, percent, col_percent, row_percent ) ) )
	ft[ c( T, out[ -1, 1] != out[-nrow(out), 1] ), , side = "top"] = borderSolid()
	ft[ line_markers, 1, side = "bottom"] = borderSolid()
	
	ft[,1] = textBold()
	ft[to="header"] = textBold()
	
	ft
}


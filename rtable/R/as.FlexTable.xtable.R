#' @import ReporteRs
#' @import xtable
#' @title get FlexTable from a xtable object
#'
#' @description Get a \code{\link{FlexTable}} object from 
#' a \code{\link{xtable}} object.
#' 
#' @param x \code{xtable} object to get \code{FlexTable} from
#' @param text.properties default text formatting properties
#' @param format.args List of arguments for the formatC function. 
#' See argument \code{format.args} of \code{print.xtable}
#' @param hline.after a vector of numbers between -1 and "nrow(x)", 
#' inclusive, indicating the rows after which a horizontal line should appear.
#' see argument \code{hline.after} of \code{print.xtable}.
#' @param padding.left cells paragraphs left padding - 0 or positive integer value.
#' @param padding.right cells paragraphs right padding - 0 or positive integer value.
#' @param ... further arguments, not used. 
#' @return a \code{\link{FlexTable}} object
#' @examples
#' #
#' @example examples/xtable.R
#' @export 
as.FlexTable.xtable = function( x, text.properties = textNormal(), 
		format.args = list(),
		hline.after = getOption("xtable.hline.after", c(-1,0,nrow(x))), 
		padding.left = 4, padding.right = 4, ... ){
	
	if( ! is.null(hline.after) )
		if (any(hline.after < -1) | any(hline.after > nrow(x))) {
			stop("'hline.after' must be inside [-1, nrow(x)]")
		}
	
	
	cols = matrix("", ncol = ncol( x )+1, nrow = nrow( x ) )
	cols[, 1] = row.names( x )
	pos <- 1
	
	varying.digits <- is.matrix( attr( x, "digits", exact = TRUE ) )

	for(i in 1:ncol(x)) {
		xcol <- x[, i]
		if(is.factor(xcol)) xcol <- as.character(xcol)
		if(is.list(xcol)) xcol <- sapply(xcol, unlist)
		ina <- is.na(xcol)
		is.numeric.column <- is.numeric(xcol)
		
		if(is.character(xcol)) {
			cols[, i+pos] <- xcol
		} else {
			if (is.null(format.args)){
				format.args <- list()
			}
			if (is.null(format.args$decimal.mark)){
				format.args$decimal.mark <- options()$OutDec
			}
			if(!varying.digits){
				curFormatArgs <- c(list( x = xcol, 
					format = ifelse(attr(x, "digits", exact = TRUE )[i+1] < 0, "E",
							attr(x, "display", exact = TRUE )[i+1]),
					digits = abs(attr(x, "digits", exact = TRUE )[i+1])),
					format.args)
				cols[, i+pos] <- do.call("formatC", curFormatArgs)
			}else{
				for( j in 1:nrow( cols ) ) {
					curFormatArgs <- c(list(x = xcol[j],
						format = ifelse(attr(x, "digits", exact = TRUE )[j, i+1] < 0,
									"E", attr(x, "display", exact = TRUE )[i+1]),
						digits = abs(attr(x, "digits", exact = TRUE )[j, i+1])),
						format.args )
					cols[j, i+pos] <- do.call("formatC", curFormatArgs)
				}
			}
		}
		## End Ian Fellows changes
		
		if ( any(ina) ) cols[ina, i+pos] <- ""
		
	}
	
	
	ft = FlexTable( data = cols, header.columns = FALSE, 
			add.rownames = FALSE, 
			header.text.props = text.properties, body.text.props = text.properties )
	ft = addHeaderRow( ft, c( "", names( x ) ) )
	
	ft = setFlexTableBorders( ft, header = T, body = T, 
			inner.vertical = borderNone(), 
			inner.horizontal = borderNone(), 
			outer.horizontal = borderNone(), 
			outer.vertical = borderNone())		

	align = attr(x, "align")
	parProp = parLeft(padding.left = padding.left, padding.right = padding.right)
	
	if( any( align == "|") ){
		new_align = character(0)
		border_right_pos = integer(0)
		do_left_table = FALSE
		do_right_table = FALSE
		for( i in seq_along(align)){
			if( i == 1 && align[i] == "|" ){
				do_left_table = TRUE
			} else if( align[i] == "|" && i < length(align) ){
				border_right_pos = append( border_right_pos, length(new_align) )
			} else if( align[i] == "|" && i == length(align) ){
				do_right_table = TRUE
			} else {
				new_align = append( new_align, align[i] )
			}
		}

		align = new_align
		if( do_left_table ) {
			ft[, 1, side = "left"] = borderSolid( )
			ft[, 1, to = "header", side = "left"] = borderSolid( )
		} 
		if( length( border_right_pos) > 0 ){
			ft[, border_right_pos, side = "right"] = borderSolid( )
			ft[, border_right_pos, to = "header", side = "right"] = borderSolid( )
		}
		if( do_right_table ) {
			ft[, length(align), side = "right"] = borderSolid()
			ft[, length(align), to = "header", side = "right"] = borderSolid()
		} 
		
	}
	
	if (!is.null(hline.after)){
		if( -1 %in% hline.after ){
			hline.after = setdiff( hline.after, -1 )
			ft[, , to = "header", side = "top"] = borderSolid()
		}
		if( 0 %in% hline.after ){
			hline.after = setdiff( hline.after, 0 )
			ft[, , to = "header", side = "bottom"] = borderSolid()
		}
		if( length( hline.after ) > 0 )
			ft[hline.after, , side = "bottom"] = borderSolid()
	}
	
	
	for( i in seq_along(align) ){
		if( align[i]=="r"){
			ft[, i] = chprop(parProp, text.align="right")
			ft[, i, to = "header"] = chprop(parProp, text.align="right")
		} else if( align[i]=="c"){
			ft[, i] = chprop(parProp, text.align="center")
			ft[, i, to = "header"] = chprop(parProp, text.align="center")
		} else {
			ft[, i] = parProp
			ft[, i, to = "header"] = parProp
		} 
	}
	
	ft
	
}

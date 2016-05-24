#' @import tidyr
#' @title creates a FlexTable by transposing variables
#'
#' @description
#' creates a \code{\link{FlexTable}} by transposing variables into observations.
#'
#' @param dataset a tidy data.frame.
#' @param id variable whose values are used to identify transposed variables
#' @param columns variables to transpose
#' @param transpose transpose each \code{columns} by \code{transpose} groups
#' @param columns.transpose transpose also \code{columns}
#' @param color optional, font color, named character vector, names are specifying
#' which element of \code{columns} to colorize and values column names
#' that contains color values to use.
#' @param background.color optional, cell background color, named character vector,
#' names are specifying which element of \code{columns} to colorize and values
#' column names that contains color values to use.
#' @param space.table add a space after each group
#' @param column.label label for column name of stacked \code{columns} when
#' \code{columns.transpose} is \code{FALSE}
#' @return a \code{\link{FlexTable}} object
#' @examples
#' #
#' @example examples/FlexPivot.R
#' @export
FlexPivot = function( dataset, id, transpose, columns,
		columns.transpose = T, color = character(0),
		background.color = character(0),
		space.table = F, column.label = "Stat." ){

	if( missing( dataset ) ) stop("dataset can't be missing")
	if( missing( id ) ) stop("id can't be missing")
	if( missing( transpose ) ) stop("transpose can't be missing")
	if( missing( columns ) ) stop("columns can't be missing")

	if( !is.data.frame( dataset ) ) stop("dataset must be a data.frame")
	if( !is.character( id ) ) stop("id must be a character vector")
	if( !is.character( transpose ) ) stop("transpose must be a character vector")
	if( !is.character( columns ) ) stop("columns must be a character vector")

	if( length( id ) < 1 ) stop("id must contain at least one variable")
	if( length( columns ) < 1 ) stop("columns must contain at least one variable")
	dataset = as.data.frame( dataset )

	color_data = get_color_dataset( color = color, columns = columns,
		dataset = dataset, default_col = "black" )
	dimnames(color_data) = list(NULL, paste0("color_", columns) )

	bgcolor_data = get_color_dataset( color = background.color, columns = columns,
		dataset = dataset, default_col = "white" )
	dimnames(bgcolor_data) = list(NULL, paste0("bgcolor_", columns) )


	color_data = cbind( dataset, color_data, stringsAsFactors = F )
	bgcolor_data = cbind( dataset, bgcolor_data, stringsAsFactors = F )

	check_data = as.data.frame( table( dataset[ , c( transpose, id) ] ) )
	if( any( check_data$Freq > 1 )){
		stop("non unique values found when pivoting dataset.")
	}

	colors = table_transpose( dataset = color_data,
		id = id,
		transpose = transpose,
		columns = paste0("color_", columns),
		columns.transpose = columns.transpose)
	bgcolors = table_transpose( dataset = bgcolor_data,
		id = id,
		transpose = transpose,
		columns = paste0("bgcolor_", columns),
		columns.transpose = columns.transpose)

	data = table_transpose( dataset = dataset,
		id = id,
		transpose = transpose,
		columns = columns,
		columns.transpose = columns.transpose)

	for(i in seq_along(id)){
		colors[,i] = rep("back", nrow(colors))
		bgcolors[,i] = rep("white", nrow(bgcolors))
	}

	add.c = length(columns) > 1

	if( !columns.transpose && add.c ){
		names(data)[is.element(names(data), "column_name")]= "*"
	} else{
		data = data[, !is.element(names(data), "column_name")]
	}

	colheaders = names(data)
	colheaders = strsplit( colheaders, "\\$_\\$" )
	.l = max( sapply( colheaders, length ) )
	colheaders = lapply( colheaders, function(x, n){
		out = character(n)
		out[seq_along(x)] = x
		if( length( x ) < n ) out = rev(out)
		out
	}, n = .l)
	colheaders = do.call( cbind, colheaders)

	if( !columns.transpose && !add.c ){
		baseline = 1:nrow(colheaders)
	} else{
		baseline = 1:(nrow(colheaders)-1)
	}

	if( space.table ){
		attr( data, "headers") = space_header_columns(headers = colheaders, baseline = baseline, content = data )
		attr( data, "color") = space_header_columns(headers = colheaders, baseline = baseline, content = colors )
		attr( data, "bgcolor") = space_header_columns(headers = colheaders, baseline = baseline, content = bgcolors )
	} else {
		data = as.matrix( data )
		mode( data ) = "character"
		data[is.na(data)]=""
		attr( data, "headers") = list( header = colheaders, body = data)
		attr( data, "color") = list( header = colheaders, body = colors)
		attr( data, "bgcolor") = list( header = colheaders, body = bgcolors)
	}
	body = attr( data, "headers")$body
	headers = attr( data, "headers")$header
	if( !columns.transpose && add.c ){
		body[, headers[nrow(headers), ] == "*"]= gsub("^\\.[0-9]+\\_", "", body[, headers[nrow(headers), ] == "*"] )
		headers[nrow(headers), headers[nrow(headers), ] == "*"] = column.label
	}

	ft = FlexTable( body, header.columns = F )

	for(i in seq_len( nrow( headers ) ) ){
		headers[i,] = gsub("^\\.[0-9]+\\_", "", headers[i,] )

		.rle = rle( headers[i,])
		ft = addHeaderRow( ft, value = .rle$values, colspan = .rle$lengths )
	}

	ft = setFlexTableBorders(ft, body = F, header=F, footer=T,
		outer.horizontal = borderNone(), outer.vertical = borderNone(),
		inner.vertical = borderNone(), inner.horizontal = borderNone())

	if( space.table ){
		ft = setFlexTableBorders(ft, body = T, header=F, footer=F,
				outer.horizontal = borderNone(), outer.vertical = borderNone(), inner.vertical = borderNone(), inner.horizontal = borderSolid())
		ft = setFlexTableBorders(ft, body = F, header=T, footer=F,
				outer.horizontal = borderSolid(), outer.vertical = borderNone(),
				inner.vertical = borderNone(), inner.horizontal = borderSolid())
		for(i in seq_len( nrow( headers ) ) ){
			ft[i, headers[i,] == " ", to = "header", side = "bottom"] = borderNone()
		}
		id_cols = 1:max( which( headers[nrow(headers),] %in% id ) )
		for(i in seq_len( nrow( headers )-1 ) ){
			ft[i, id_cols, to = "header", side = "bottom"] = borderNone()
		}
		for(i in seq_len( nrow( headers ) ) ){
			ft[, headers[i,] == " ", side = "bottom"] = borderNone()
		}
	} else {
		ft = setFlexTableBorders(ft, outer.horizontal = borderSolid( ), outer.vertical = borderSolid( ),
				inner.vertical = borderSolid(), inner.horizontal = borderSolid())
	}
	ft[to="header"] = parCenter()
	ft[] = parCenter()

	for( j in id ){
		h = attr( data, "headers")$header

		.w = which( h[nrow(h), ] %in% id  )
		for( ..w in .w )
			ft = spanFlexTableRows(ft, j = ..w, runs = attr( data, "headers")$body[,..w] )
	}

	colors = attr( data, "color")$body
#	colors[is.na(colors) ] = "black"
	colors[colors == " "] = "black"

	bgcolors = attr( data, "bgcolor")$body
	bgcolors[is.na(bgcolors) ] = "white"
	bgcolors[bgcolors %in% c(" ", "")] = "white"

	for(i in seq_len(nrow(colors) ) ){
		for(j in seq_len(ncol(colors) ) ){
			if( !is.na(colors[i,j]) ) try({ft[i,j] = textNormal(color=colors[i,j])}, silent=T)
		}
	}
	for(j in seq_len(ncol(bgcolors) ) ){
		bgc = bgcolors[,j]
		try({ft = setFlexTableBackgroundColors(ft, , j, bgc, to = "body")}, silent=T)
	}

	ft[length(ft), side = "bottom"] = borderSolid()
	ft
}

table_transpose = function( dataset, id, transpose, columns, columns.transpose ){

	if( missing( dataset ) ) stop("dataset can't be missing")
	if( missing( id ) ) stop("id can't be missing")
	if( missing( transpose ) ) stop("transpose can't be missing")
	if( missing( columns ) ) stop("columns can't be missing")

	if( !is.data.frame( dataset ) ) stop("dataset must be a data.frame")
	if( !is.character( id ) ) stop("id must be a character vector")
	if( !is.character( transpose ) ) stop("transpose must be a character vector")
	if( !is.character( columns ) ) stop("columns must be a character vector")

	if( length( id ) < 1 ) stop("id must contain at least one variable")
	if( length( columns ) < 1 ) stop("columns must contain at least one variable")

	subset = dataset[, columns, drop = F ]
	names( subset ) = paste0( ".", seq_len( length(columns) ), "_", columns )

	data = cbind( dataset[, id, drop = F ], subset )

	if( length( transpose ) > 0 ){
		data_transpose = dataset[, transpose, drop = F ]
		data_transpose[,1] = factor( data_transpose[,1], levels = unique(data_transpose[,1]) )
		levels(data_transpose[,1]) = paste0( ".",
				sprintf("%03d", seq_along( levels(data_transpose[,1]) )),
				"_",
				levels(data_transpose[,1]) )
		data_transpose[,1] = as.character(data_transpose[,1])
		data = cbind( data, data_transpose )
	}
	data = gather_(data, key_col = "column_name", value_col = "value", gather_cols = names( subset ), convert = TRUE )
	if( any( is.element( names(data), "variable") ) )
		names(data)[is.element( names(data), "variable")] = "column_name"

	key_col = "keycol"

	if( columns.transpose ){
		newdata = unite_(data, col = key_col,
			from = c(transpose, "column_name"), sep = "$_$" )
		data = spread_(newdata, key_col = key_col, value_col = c("value") )
	} else if(!columns.transpose && length(transpose) < 1 ){
		newdata = data
	} else {
		newdata = unite_(data, col = key_col, from = transpose, sep = "$_$" )
		data = spread_(newdata, key_col = key_col, value_col = c("value") )
	}
	for( i in seq_len( ncol( data ) ) ){
		if( is.factor( data[, i] ) ) data[, i] = as.character(data[, i])
	}

	data
}

space_header_columns = function( headers, baseline = seq_len( nrow(headers) ), content ){

	fakemeta = matrix( "", ncol = 0, nrow = nrow( headers ) )
	fakemeta = cbind(fakemeta, headers[ , 1] )

	fakedata = matrix( "", ncol = 0, nrow = nrow( content ) )
	fakedata = cbind(fakedata, content[ , 1] )

	concern = is.element(seq_len( nrow(headers) ), baseline)

	for(j in setdiff( seq_len( ncol(headers) ), 1 ) ){

		temp = headers[, j] == headers[, j-1]
		temp = ifelse( temp, headers[, j], "$_$" )
		temp = cbind( temp, headers[, j, drop = FALSE] )
		fakemeta = cbind(fakemeta, temp )
	}
	dimnames(fakemeta) = list(NULL, NULL)

	for(j in setdiff( seq_len( ncol(headers) ), 1 ) ){
		temp = rep( "$_$", nrow(content) )
		temp = cbind( temp, content[, j, drop = FALSE] )
		fakedata = cbind(fakedata, temp )
	}
	fakedata[is.na(fakedata)]=""

	mustbetaken = !apply( fakemeta[!concern, , drop = F], 2, function( x ) any( x == "$_$" ) )
	keepspace = apply( fakemeta[concern, , drop = F], 2, function( x ) any( x == "$_$" ) )

	fakemeta = fakemeta[, mustbetaken | keepspace, drop = F ]
	fakemeta[fakemeta=="$_$"]=" "
	fakedata = fakedata[, mustbetaken | keepspace, drop = F]
	fakedata = as.matrix( fakedata )
	fakedata[fakedata=="$_$"]=" "
	list( header = fakemeta, body = fakedata)
}


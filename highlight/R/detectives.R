
#' Simple detective
#' 
#' This detective only uses semantic information to make its 
#' investigation. 
#' @param x output of the parser. The detective is only interested in the 
#'          \samp{token} column of the data.
#' @param \dots ignored
#' @return a vector of styles grouping similar tokens together
#' @examples
#' \dontrun{
#' p <- parse( text = deparse( jitter ), keep.source=TRUE )
#' simple_detective( p )
#' }
#' @export
simple_detective <- function( x, ...){
	
	data <- getParseData( x )
	desc   <- as.character( data[ data[["terminal"]], "token" ] )
	styles <- character( length( desc ) )
	
	styles[ desc == "COMMENT"  ] <- "comment"
	styles[ desc == "ROXYGEN_COMMENT" ] <- "roxygencomment"
	
	styles[ grepl( "^'.*?'$", desc ) ] <- "keyword"
	styles[ desc %in% c( "FUNCTION", "FOR", "IN", "IF", 
		"ELSE", "WHILE", "NEXT", "BREAK", "REPEAT", 
		"AND", "AND2", "OR", "OR2", "GT", 
		"LT", "GE", "LBB", "NE", "SPECIAL", 
		"NS_GET_INT", "NS_GET") ] <- "keyword" 
	
	styles[ desc == "STR_CONST" ] <- "string"
	styles[ desc == "NUM_CONST" ] <- "number"
	
	styles[ desc == "SYMBOL_FUNCTION_CALL" ] <- "functioncall"
	styles[ desc %in% c("SYMBOL_SUB", "EQ_SUB" )  ] <- "argument"
	styles[ desc == "SYMBOL_PACKAGE" ] <- "package"
	
	styles[ desc %in% c("SYMBOL_FORMALS") ] <- "formalargs" 
	styles[ desc %in% "EQ_FORMALS" ] <- "eqformalargs" 
	
	styles[ desc %in% c("EQ_ASSIGN", "LEFT_ASSIGN" )] <- "assignement"
	styles[ desc == "SYMBOL" ] <- "symbol"
	styles[ desc == "SLOT" ] <- "slot"
	styles
	
}


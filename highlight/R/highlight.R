#' Syntax Highlighter for R
#' 
#' Syntax highlighter for R based on output from the R parser
#' 
#' @seealso
#' 	The main function of the package is \code{\link{highlight}}. 
#' 	
#' 	\code{\link{highlight}} delegates rendering the document to 
#' 	\code{\link{renderer}}s, such as the \code{\link{renderer_latex}} 
#' 	or \code{\link{renderer_html}} and is helped by a
#' 	detective to make sense of the results
#' 	from the parser. The package ships a \code{\link{simple_detective}}. 
#' 	
#' 	The package also defines a custom sweave driver 
#' 	(\code{\link{HighlightWeaveLatex}}) for latex based 
#' 	on the standard sweave latex driver (\code{\link[utils]{RweaveLatex}})
#' 	using \code{\link{highlight}} to perform syntax 
#' 	highlighting of R code chunks. 
#' @examples
#' \dontrun{
#' tf <- tempfile()
#' dump( "glm" , file = tf )
#' 
#' # rendering in html
#' highlight( tf, output = stdout(), 
#' 	renderer = renderer_html() )
#' 
#' # rendering in latex
#' highlight( tf, output = stdout(), 
#' 	renderer = renderer_latex() )
#' 
#' # Sweave driver using syntax highlighting
#' if( require( grid ) ){
#' 	v <- vignette( "grid", package = "grid" )$file
#' 	file.copy( v, "grid.Snw" )
#' 	Sweave( "grid.Snw", driver= HighlightWeaveLatex() )
#' 	system( "pdflatex grid.tex" )
#' 	if (.Platform$OS.type == "windows"){ 
#' 		shell.exec( "grid.pdf" )
#' 	} else {
#' 		system(paste(shQuote(getOption("pdfviewer")), "grid.pdf" ), 
#' 			wait = FALSE)
#' 	}
#' }
#' 
#' unlink( tf )
#' }
#' @docType package
#' @name highlight-package
NULL
	
subsetParseData <- function( p, i = 0, styles){
    data <- getParseData(p)
    data$styles <- rep("", nrow(data) )
    data$styles[ data$terminal ] <- styles

    if( is.null(i) || i == 0 ){
        return(data)
    }

    srcref <- attr(p, "srcref")[[i]]
    line1 <- srcref[1L]
    line2 <- srcref[3L]
    col1  <- srcref[5L]
    col2  <- srcref[6L]

    if(line1 == line2){
        data <- data[ data$line1 == line1 & data$col1 >= col1 & data$col2 <= col2, ]
    } else {
        data <- data[
            ( data$line1 > line1  & data$line2 < line2 ) |
            ( data$line1 == line1 & data$col1 >= col1 ) |
            ( data$line2 == line2 & data$col2 <= col2 )
        ,
        ]
    }
    data
}


#' syntax highlighting based on the R parser
#' 
#' The \code{highlight} function performs syntax highlighting based on the 
#' results of the \code{\link[base]{parse}} and the investigation
#' of a detective.
#' 
#' @param file code file to parse. This is only used if the \code{parse.output} is given
#' @param output where to write the rendered text. If this is anything else than the 
#' default (standard output), the \code{\link{sink}} function
#' is used to redirect the standard output to the output.
#' @param detective the detective chooses the style to apply to each token, basing its 
#' investigation on the results of the \code{\link[base]{parse}}
#' @param renderer highlight delegates rendering the information to the renderer. This 
#' package includes html and latex renderers. See \code{\link{renderer_html}}
#' and \code{\link{renderer_latex}}
#' @param encoding encoding to assume for the file. the argument is directly passed 
#' to the \code{\link[base]{parse}}.
#' @param parse.output output from the \code{\link[base]{parse}}. If this is given, the 
#' arguments \code{file} and \code{encoding} are not used
#' @param styles result of the detective investigation. A character vector 
#' with as many elements as there are tokens in the parser output
#' @param expr In case we want to render only one expression and not the full parse
#' tree, this argument can be used to specify which expression
#' to render. The default (NULL) means render all expressions. This 
#' feature is used by the sweave driver shipped with this package. See
#' \code{\link{HighlightWeaveLatex}}
#' @param final.newline logical. Indicates if a newline character is added after all tokens.
#' @param showPrompts if TRUE, the highlighted text will show standard and continue prompt
#' @param prompt standard prompt
#' @param continue continue prompt
#' @param initial.spaces should initial spaces be displayed or skipped.
#' @param size font size. only respected by the latex renderer so far.
#' @param show_line_numbers logical. When TRUE, line numbers are shown in the output.
#' @param \dots additional arguments, currently ignored. 
#'               
#' @return The resulting formatted text is returned invisibly. It is also 
#' written to the output if the output is not \code{NULL}
#' @seealso 
#' 	\code{\link{renderer_html}} and \code{\link{renderer_latex}} are the
#' 	two implementation of renderers currently available in this package. 
#' 	
#' 	\code{\link{simple_detective}} is an example detective which does a very 
#' 	simple investigation.
#' 
#' @examples
#' \dontrun{
#' 	tf <- tempfile()
#' 	dump( "jitter", file = tf )
#' 	highlight( file = tf, detective = simple_detective, 
#' 		renderer = renderer_html( document = TRUE ) )
#' 	highlight( file = tf, detective = simple_detective, 
#' 		renderer = renderer_latex( document = TRUE ) )
#' 	
#' }
#' @export
highlight <- function( file, output = stdout(),
    detective = simple_detective, renderer, encoding = "unknown",
    parse.output = parse( file, encoding = encoding, keep.source = TRUE ),
    styles = detective( parse.output ),
    expr = NULL,
    final.newline = FALSE,
    showPrompts = FALSE,
    prompt = getOption( "prompt" ) ,
    continue = getOption( "continue"),
    initial.spaces = TRUE,
    size = NULL,
    show_line_numbers = FALSE,
    ... ){

    size <- match.arg( size )
    # forcing the arguments in a certain order
    force( parse.output )
    force( styles )
    force( renderer )

    # only terminal symbols matter
    data   <- subsetParseData( parse.output, expr, styles )
    data$top_level <- .Call( "top_level", data$parent, PACKAGE = "highlight" )
    data   <- data[ data[["terminal"]], ]

    # let the renderer do its thing
    data$ftokens <- renderer$formatter(
        tokens = renderer$translator( as.character( data[, "text"] ), size = size ),
        styles = data[, "styles"] )

    # useful to only render a given expression and not all of them.
    # this is mainly used in the sweave driver
    startline <- if( !is.null( expr ) ) as.integer( data[1, "line1" ] ) else 1L

    line_numbers <- seq( startline, max(data$line2))
    width <- max( nchar( line_numbers ) )
    line_numbers <- renderer$formatter(
        sprintf( sprintf( "%%0%dd  ", width ), line_numbers ),
        rep( "line", length(line_numbers) )
    )

    # paste everything together in C++ using Rcpp
    highlighted_text <- c( if( !is.null(renderer$header) ) renderer$header(),
        .Call( "get_highlighted_text",
            data,
            startline,
            max(data$line2) ,
            renderer$space(),
            renderer$newline(),
            if( showPrompts) renderer$formatter( renderer$translator( prompt, size = size ) , "prompt" ) else "",
            if( showPrompts) renderer$formatter( renderer$translator( continue, size = size ) , "prompt" ) else "",
            initial.spaces = initial.spaces,
            line_numbers,
            isTRUE( show_line_numbers),
            PACKAGE = "highlight"
        ),
        if( !is.null(renderer$footer) ) renderer$footer() )

    # maybe write the result to the output connection
    if( !is.null(output) ){
        try( writeLines( highlighted_text, output, sep = "" ), silent = TRUE )
    }
    invisible( highlighted_text )
}
fm <- formals(highlight)
fm[[ which( names(fm) == "size") ]] <- LATEX_SIZES
formals( highlight ) <- fm

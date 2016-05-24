
highlight_supported_languages <- function(){
    files <- list.files( 
        system.file( "highlight", "langDefs", package = "highlight" ),   
        pattern = "lang$" )
    gsub( "[.]lang$", "", files )
}

#' List of themes supported by external_highlight
#' 
#' List of themes supported by \code{\link{external_highlight}}
#' 
#' @return A character vector with the names of the themes
#' @export
highlight_themes <- function(){
    files <- list.files( 
        system.file( "highlight", "themes", package = "highlight" ),   
        pattern = "style$" )
    gsub( "[.]style$", "", files )
}

#' List of available output types supported by external_highlight
#' 
#' List of available output types supported by \code{\link{external_highlight}}
#' 
#' @return A character vector with the list of supported types
#' @export
highlight_output_types <- function(){
    c("HTML","XHTML","TEX","LATEX","RTF","XML","ANSI","XTERM256",
        "HTML32", "SVG","BBCODE" )    
}

highlight_theme <- function( theme = "emacs" ){
    if( missing(theme) ){
        theme <- highlight_themes()[1L]
    } else {
        theme <- match.arg( theme, highlight_themes() )
    }
    system.file( "highlight", "themes", sprintf( "%s.style", theme ), package = "highlight" )
}

highlight_lang <- function( lang = highlight_supported_languages() ){
    if( missing(lang)){
        stop( "no language" )
    } else {
        lang <- match.arg(lang, highlight_supported_languages() )
    }
    system.file( "highlight", "langDefs", sprintf("%s.lang", lang), package = "highlight" ) 
}

highlight_type <- function(type = highlight_output_types() ){
    if( missing( type ) ){ type <- "HTML" }
    type <- match.arg( type, highlight_output_types() )
    match( type, highlight_output_types() ) - 1L
}

#' Multi-language source code highlighter
#' 
#' Multi-language source code highlighter
#' 
#' @param file Source file to highlight
#' @param outfile Destination of the highlighted code. 
#'       When \code{NULL}, the code is simply returned as a character vector
#' @param theme One of the themes. See \code{\link{highlight_themes}} for the list
#'              of available themes.
#' @param lang The language in which the code is to be interpreted. If this argument
#'             is not given, it will be deduced from the file extension.
#' @param type Output format. See \code{\link{highlight_output_types}} for the list 
#'             of supported output types.
#' @param line_numbers if \code{TRUE}, the result will include line numbers
#' @param doc if \code{TRUE}, the result is a stand alone document, otherwise, just a 
#'            portion to include in a document
#' @param code If given, then the source code is not read from the file
#' 
#' @return Nothing if \code{outfile} is given, with the side effect of writing into the file. 
#' The result as a character vector if outfile is NULL
#' @seealso \code{\link{highlight}} to highlight R code using the information from the parser
#' @export
external_highlight <- function( file, 
    outfile = stdout(), 
    theme = "kwrite",
    lang  = NULL , 
    type  = "HTML", 
    line_numbers = FALSE, 
    doc = TRUE, 
    code
){
        
    if( !missing(code) ){
        file <- sprintf( "%s.%s", tempfile(), lang )
        writeLines( code, file )    
    }
    type  <- highlight_type(type)
    theme <- highlight_theme(theme) 
    
    lang <- highlight_guess_language(file, lang = lang)
    lang <- highlight_lang(lang)
    
    using_tempfile <- is.null(outfile) || !is.character(outfile)
    output_file <- if( using_tempfile ) tempfile() else outfile
    .Call( "HighlightMain", file, output_file, type, theme, lang, 
        isTRUE(line_numbers), 
        isTRUE(doc), 
        PACKAGE = "highlight"
        )
    code <- readLines(output_file)
    
    w <- which( code == "\\mbox{}") ; code <- code[ - tail(w,1) ]
    w <- tail( grep( "\\\\\\\\$", code ), 1 )
    code[w] <- gsub( "\\\\\\\\$", "", code[w] ) 
    
    if( !is.null(outfile) ) writeLines( code, outfile )
    invisible(code)
}


highlight_extensions <- function(){
    txt <- readLines( system.file( "highlight", "filetypes.conf", package = "highlight" ) )
    
    df <- do.call( rbind, lapply( grep( "^[$]ext" , txt, value = TRUE ), function(x) {
        
        extensions <- strsplit( sub( "^.*=", "",  x), " ")[[1]]
        language   <- sub("^.*[(](.*)[)].*$", "\\1", x  )
        
        data.frame( 
            lang = rep( language, length(extensions)+1L ), 
            ext = c( language, extensions ), 
            stringsAsFactors= FALSE 
       )
    } ) )
    
    
    files <- list.files( system.file( "highlight", "langDefs", package = "highlight" ), pattern = "[.]lang$" )
    languages <- sub( "[.]lang$", "", files )
    
    missings <- setdiff( languages, unique( df$lang ) )
    df <- rbind( df, data.frame( lang = missings, ext = missings, stringsAsFactors = FALSE ) )
    
    df <- df[ order(df$lang), ]
    
    
}

highlight_guess_language <- function(file, lang = NULL){
    if( is.null(lang)) lang <- sub( "^.*[.]([^.]*)$", "\\1", file )
    if( lang == "" ) stop( "no extension" ) 
    
    df <- highlight_extensions()
    id <- match( lang, df$ext )
    if( is.na(id ) ) stop( "unknown extension" )
   
    df[ id, "lang" ]
}


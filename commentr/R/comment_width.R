



################################################################################
#                                                                              #
#                           Exported helper function                           #
#                                                                              #
################################################################################



#' Function to specify the width of the comment
#'
#' Helper function to the commet family of functions. Can take both numerical values as the width
#' in characters or a prespecified standard through a character string.
#'
#' The width is specified in order to fit the comment on a printed page or in a script file showed in RStudio
#' (hence  the -5 in \code{getOption("width") - 5} to adjust for the line numbering in the script file compared to the console.
#' It is also possible to use a global option \code{comment_width}. This is not a standard option and has to be manually
#' specified (for example in a .Rprofile-file).
#'
#' @param width is the prefered width. Could be numeric (number of characters), a paper size
#' (currently \code{"a4portrait"}, \code{"a4landscape"}, \code{"a3portrait"} or \code{"a3landscape"}),
#'  \code{"script_width"} (= \code{getOption("width") - 5}) or \code{"option"} (to get data from global option \code{"comment_width"}).
#'  Default is \code{"option"} but if a global option does not exist, "a4portrait".
#' @param a4portrait_width specifies the number of characters that can be printed on a single line on a a4 paper
#' in portrait orientation. This value is usually 80 (or in the range from 60 to 75). Here a lower value is set by default
#' due to experimentation on the authors own computer. Please contact the author if this value seems strange.
#' Note however that it is usully more sufficient to use a global option for the \code{width} parameter than to
#' manually change this value (whih is however possible for increased flexibility). This value depends on margins and font size
#' when printing.
#'
#' @return An integer specifiing the text width (in number of characters) to be used in a comment.
#' @export
#' @examples
#' comment_width()
#' comment_width(42)
#' comment_width("a4portrait")
#'
#'\dontrun{
#' # We can set a global option for the comment_width
#' options(comment_width = 80)
#' header_comment("Test", "A small test")
#' }
#'
comment_width <- function(width = "option", a4portrait_width = 80){
  
  # Check that width is a valid argument
  valid_width_arguments <- c("option", "script_width", "a4portrait", "a4landscape", "a3portrait", "a3landscape")
  if ( !(is.numeric(width) | width %in% valid_width_arguments)){
    stop("Invalid width argument! See ?comment_width for help!")
  }
  
  # Use global option if possible, otherwise use a4portrait as default
  if (width == "option" & !is.null(getOption("comment_width"))){
    width <- getOption("comment_width")
  } else if (width == "option" & is.null(getOption("comment_width"))){
    width <- "a4portrait"
  }
  
  # a4portrait is usually recommended to 80 but 53 was the value found when experimenting
  a4portrait   <- a4portrait_width
  golden_ratio <- 8 / 5
  
  # If width_text is numeric, it should be translated into a numeric value.
  # If it is not a character, it is preserved and returned unchanged
  text2num <- function(width_text){
    if (is.character(width_text)){
      switch(width_text,
             "script_width" = getOption("width") - 5,
             "a4portrait"   = a4portrait,
             "a4landscape"  = a4portrait * golden_ratio,
             "a3portrait"   = a4portrait * golden_ratio,
             "a3landscape"  = 2 * a4portrait)
    } else{
      width_text
    }
  }
  
  if (width %in% valid_width_arguments){
    width <- text2num(width)
  }
  
  
  min(width, text2num("script_width"))
  
}



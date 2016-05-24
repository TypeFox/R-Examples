


#' @rdname comment
#' @export
block_comment <- function(
  description         = "",
  empty_lines_first   = 1,
  empty_lines_last    = empty_lines_first,
  allign              = "center",
  token               = "#",
  html                = FALSE,
  clipboard           = TRUE,
  verbose             = TRUE,
  ...
) {
  
  ##  Specify width
  width <- comment_width(...)
  
  # Translate the allign argument to the side argument for str_pad
  # The names would be confusing without this recoding
  side <- switch(allign, "right" = "left", "left" = "right", "center" = "both")
  
  ## Create the empty lines to add first and last
  elf <- paste(rep(empty_line(width, token = token), empty_lines_first), collapse = "\n")
  ell <- paste(rep(empty_line(width, token = token), empty_lines_last), collapse = "\n")
  
  msg <- paste(
    full_line(width, token = token),
    elf,    
    ## Conditions used to specify if the text should be on just one alligned line or split
    ## onto several lines.
    if (str_length(description) <= width - 4){
      paste(token, str_pad(description, width - 4, side = side), token)
    } else{
      x <- strwrap(description, width = width - 2, prefix = token, indent = 1, exdent = 1)
      paste(str_pad(x, width - 2, side = "right"), token, collapse = "\n")
    },
    ell,
    full_line(width, token = token),
    sep = "\n"
  )
  
  out(msg, html, clipboard, verbose)
  
}

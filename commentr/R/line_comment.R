
#' @rdname comment
#' @export
line_comment <- function(title = "", ..., clipboard = TRUE, verbose = TRUE, token = "#", html = FALSE){
  
  # Space before and after if any text
  if (title != ""){
    title <- paste0(" ", title, " ")
  }
  
  # Use block_comment and give warning if to long text
  if (stringr::str_length(title) > comment_width(...)){
    warning("Comment too long to fit on one line. 'block_comment' used instead!")
    block_comment(title, token = token, html = html, 
                  clipboard = clipboard, verbose = verbose)
  
  # Otherwise give comment as one liner
  } else{
   msg <- paste(
    stringr::str_pad(title, comment_width(...), side = "both", pad = token),
    sep = "\n")
   out(msg, html, clipboard, verbose)
  }
}

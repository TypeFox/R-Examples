

################################################################################
#                                                                              #
#                               Helper functions                               #
#                                                                              #
################################################################################


#' Helper function to create string with repeated hashtags
#' @param ... arguments passed to \code{comment_width}
#' @param token a character specifying the comment symbol. "#" by default.
#' @keywords internal
full_line <- function(..., token = "#"){
  paste0(rep_len(token, comment_width(...)), collapse = "")
}


#' Helper function to create string with hashtags and spaces
#' @param ... arguments passed to \code{comment_width}
#' @param token a character specifying the comment symbol. "#" by default.
#' @keywords internal
empty_line <- function(..., token = "#"){
  paste0(token, str_pad(" ", comment_width(...) - 2, side = "both"), pad = token, collapse = "")
}


#' Comment start for R Markdown/HTML
#' @param html should the comment be used in HTML (FASLE by default)
#' @keywords internal
comment_start <- function(html = FALSE){
  if (html) "<!--" else ""
}

#' Comment end for R Markdown/HTML
#' @param html should the comment be used in HTML (FASLE by default)
#' @keywords internal
comment_end <- function(html = FALSE){
  if (html) "-->" else ""
}


#' Handle output of the comment
#' @param msg the comment to return
#' @param clipboard,verbose,html See \code{\link{commentr}}
#' @keywords internal
out <- function(msg, html, clipboard = TRUE, verbose = TRUE){
 
  mac <- Sys.info()[['sysname']] == "Darwin"
    
  msg <- paste(
    comment_start(html),
    msg,
    comment_end(html),
    sep = "\n"
  )
    
  ## if on Mac OSX and flagged to true, then copy to clipboard
  if(clipboard &&  mac) {
    con <- pipe("pbcopy", "w")
    writeLines(msg, con)
    close(con)
    message("The comment has been copied to clipboard and can be pasted into a script file!") 
  }
  
  if (verbose){
    writeLines(msg)
  }
}


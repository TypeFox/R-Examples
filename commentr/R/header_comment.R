
#' @rdname comment
#' @export

header_comment <- function(
  title,
  description = "",
  author      = getOption("name"),
  contact     = getOption("contact"),
  client      = author,
  date_created= format(Sys.time(), "%Y-%m-%d"),
  date_updated= date_created,
  source      = getwd(),
  tab         = 17,
  token       = "#",
  html        = FALSE,
  clipboard   = TRUE,
  verbose     = TRUE,
  ...
){
  
  ##  Error if missing author/contact
  if (is.null(author) | is.null(contact)){ 
    stop("Author and/or contact not specified!")
  }
  
  ##  Width of the comment
  width <- comment_width(...)
  
  ##  Line with text and label
  text_line <- function(which_text, text) {
    pre  <- str_pad(paste0(token, " ", which_text, ":"), tab, side = "right")
    text <- paste0(pre, text)
    paste0(stringr::str_pad(text, width - 1, side = "right"), token)
  }
  
  ##  Line with text only
  text_paragraph <- function(text){
    x <- strwrap(text, width = width - 2, indent = tab - 1, exdent = tab - 1, prefix = token)
    paste(stringr::str_pad(x, width - 2, side = "right"), token, collapse = "\n")
  }
  
  ##  Linebreak too long text line
  text_cond_paragraph <- function(which_text, text){
    if (stringr::str_length(text) < width - tab){
      text_line(which_text, text)
    } else{
      part1 <- substr(text, 1, width - tab - 3)
      part2 <- substring(text, width - tab - 2)
      paste(
        text_line(which_text, part1),
        text_paragraph(part2),
        sep = "\n"
      )
    }
  }
  
  
  msg <- paste(
    full_line(width, token = token),
    empty_line(width, token = token),
    text_cond_paragraph("Purpose", title),
    empty_line(width, token = token),
    text_cond_paragraph("Author", author),
    text_cond_paragraph("Contact", contact),
    text_cond_paragraph("Client", client),
    empty_line(width, token = token),
    text_cond_paragraph("Code created", date_created),
    text_cond_paragraph("Last updated", date_updated),
    text_cond_paragraph("Source", source),
    empty_line(width, token = token),
    text_cond_paragraph("Comment", description), # writeLines called through text_cond_paragraph
    empty_line(width, token = token),
    full_line(width, token = token),
    sep = "\n"
  )

  ## Handle the output
  out(msg, html, clipboard, verbose)
}


#' Generate numeric XPath expression
#'
#' @title num_xpath: Generate numeric XPath expression
#' @param data the header XPath
#' @export
num_xpath <- function(data) UseMethod("num_xpath")

#' @export
num_xpath.default <- function(data){
  stop("Unknown input", call. = F)
}

#' @export
num_xpath.character <- function(data){
  return(data)
}

#' @export
num_xpath.numeric <- function(data){

  index.count <- data - 1
  index.xpath <- sapply(1:length(index.count), function(pos) {
    sprintf("count(preceding::tr) = %s", index.count[pos])
  })
  index.xpath <- sprintf("//tr[%s]", paste(index.xpath, collapse = " or "))

  return(index.xpath)
}

#' @export
num_xpath.list <- function(data){

  index.xpath <- lapply(data, num_xpath)
  return(index.xpath)
}


#' Return trindex given an XPath
#' @param table.Node the table node
#' @param xpath XPath
#' @noRd
get_trindex <- function(xpath, table.Node) UseMethod("get_trindex")

get_trindex.default <- function(xpath, table.Node){
  stop("Unknow XPath", call. = FALSE)
}

get_trindex.NULL <- function(xpath, table.Node){
  return(logical())
}

get_trindex.character <- function(xpath, table.Node){

  tr.index <- unlist(XML::xpathSApply(table.Node, xpath, XML::xmlGetAttr, "HTMLTABtrindex"))
  tr.index <- as.numeric(tr.index)
  return(tr.index)
}

get_trindex.list <- function(xpath, table.Node){

  tr.index <- lapply(seq_along(xpath), function(id) {
    XML::xpathSApply(table.Node, xpath[[id]], function(x) {
      val <- XML::xmlGetAttr(XML::xmlParent(x), "HTMLTABtrindex")
      val <- as.numeric(val)
      return(val)
    })
  })

  return(tr.index)
}

#' Return header xpath
#'
#' @param table.Node the table node
#' @param header an information for the header rows
#' @return a character vector of XPath statements
#' @noRd
get_head_xpath <- function(header, table.Node) UseMethod("get_head_xpath")

get_head_xpath.default <- function(header, table.Node) {
  stop("Unknown header information", .call = F)
}

get_head_xpath.character <- function(header, table.Node) {
  return(header)
}

get_head_xpath.numeric <- function(header, table.Node) {

  if(header[1] == 0){
    header <- NULL
  } else{
    header <- num_xpath(header)
  }

  return(header)
}

get_head_xpath.NULL <- function(header, table.Node){

  thead <- has_tag(table.Node, "//thead")
  thead.th <- has_tag(table.Node, "//thead/tr[th]")
  thead.td <- has_tag(table.Node, "//thead/tr[td]")

  tr <- has_tag(table.Node, "//tr")
  th <- has_tag(table.Node, "//tr[th and not(./td)]")
  td <- has_tag(table.Node, "//tr[td and not(./th)]")


  if(thead) {
    header.xpath <- '//thead/tr'
    return(c(header.xpath))
  }

  if (!thead && th){
    header.xpath <- "//tr[th and not(./td)]"
    return(c(header.xpath))
  }

  if (!thead && !th){
    header.xpath <- "//tr[position() = 1]"
    message("Neither <thead> nor <th> information found. Taking first table row for the header. If incorrect, specifiy header argument.")
    return(c(header.xpath))
  }
}


#' Return body xpath
#'
#' @param table.Node the table node
#' @param body an information for the body rows
#' @return a character vector of XPath statements
#' @noRd
get_body_xpath <- function(body, table.Node) UseMethod("get_body_xpath")

get_body_xpath.default <- function(body, table.Node){
  stop("Unknow body information", .call = F)
}

get_body_xpath.numeric <- function(body, table.Node){
  body <- num_xpath(body)
  return(body)
}

get_body_xpath.character <- function(body, table.Node){
  return(body)
}

get_body_xpath.NULL <- function(body, table.Node){
  tbody <- has_tag(table.Node, "//tbody")
  tbody.th <- has_tag(table.Node, "//tbody/tr[th]")
  tbody.td <- has_tag(table.Node, "//tbody/tr[td]")

  tr <- has_tag(table.Node, "//tr")
  th <- has_tag(table.Node, "//tr[th and not(./td)]")
  td <- has_tag(table.Node, "//tr[td and not(./th)]")

  if(tbody){
    body.xpath <- "//tbody/tr"
    return(c(body.xpath))
  } else {
    body.xpath <- "//tr[./td]"
    return(c(body.xpath))
  }

}

#' Assemble XPath expressions for header and body
#'
#' @param table.Node the table node
#' @param header a vector that contains information for the identification of the header row(s). A numeric vector can be specified where each element corresponds to the table rows. A character vector may be specified that describes an XPath for the header rows. If left unspecified, htmltable tries to use semantic information from the HTML code
#' @param body a vector that specifies which table rows should be used as body information. A numeric vector can be specified where each element corresponds to a table row. A character vector may be specified that describes an XPath for the body rows. If left unspecified, htmltable tries to use semantic information from the HTML code
#' @param complementary logical, should htmltab ensure complementarity of header, inbody header and
#'    body elements (default TRUE)?
#' @return a character vector of XPath statements
#' @noRd
identify_elements <- function(table.Node, header, body, complementary = T){

  header_MAIN <- header[[1]]
  header_INBODY <- header[-1]

  # switch compementary
  if(!is.null(header_MAIN[1])){
    if(header_MAIN[1] == 0){complementary <- FALSE}
  }

  # Produce XPaths
  header_MAIN.xpath <- get_head_xpath(table.Node = table.Node, header = header_MAIN)
  body.xpath <- get_body_xpath(table.Node = table.Node, body = body)
  header_INBODY.xpath <- num_xpath(header_INBODY)

  # Receive HTMLTABtrindex
  header_MAIN.trindex <- get_trindex(header_MAIN.xpath, table.Node)
  header_INBODY.trindex <- get_trindex(header_INBODY.xpath, table.Node)
  body.trindex <- get_trindex(body.xpath, table.Node)

  # Check INBODY correct
  ifstop(length(header_INBODY.xpath) > 0 && length(header_INBODY.trindex[[1]]) == 0,
         "Your specified inbody header cells could not be identified.", call. = F)

  # Complementarity
  if(isTRUE(complementary)){

    if(is.null(body)){
      body.trindex <- setdiff(body.trindex, c(header_MAIN.trindex, unlist(header_INBODY.trindex)))
    }

    if(is.null(header_MAIN) && !is.null(body)){
      header_MAIN.trindex <- setdiff(header_MAIN.trindex, c(body.trindex, unlist(header_INBODY.trindex)))
    }

  }

  trindex <- list(header = header_MAIN.trindex,
             inbody = header_INBODY.trindex,
             body = body.trindex)

  xpath <- list(header = header_MAIN.xpath,
                inbody = header_INBODY.xpath,
                body = body.xpath)

  return(list(trindex = trindex, xpath = xpath))
}

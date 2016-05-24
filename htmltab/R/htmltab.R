#' Assemble a data frame from HTML table data
#'
#' Robust and flexible methods for extracting structured information out of HTML tables
#'
#' @export
#' @param doc the HTML document which can be a file name or a URL or an already parsed document
#'    (by XML's parsing functions)
#' @param which a vector of length one for identification of the table in the document. Either a
#'    numeric vector for the tables' rank or a character vector that describes an XPath for the table
#' @param header the header formula, see details for specifics
#' @param headerFun a function that is executed over the header cell nodes
#' @param headerSep a character vector that is used as a seperator in the construction of the table's
#'    variable names (default ' >> ')
#' @param body a vector that specifies which table rows should be used as body information. A numeric
#'    vector can be specified where each element corresponds to a table row. A character vector may be
#'    specified that describes an XPath for the body rows. If left unspecified, htmltab tries to use
#'    semantic information from the HTML code
#' @param bodyFun a function that is executed over the body cell nodes
#' @param complementary logical, should htmltab ensure complementarity of header, inbody header and
#'    body elements (default TRUE)?
#' @param fillNA character vector of symbols that are replaced by NA (default c(''))
#' @param rm_superscript logical, should superscript information be removed from header and body cells
#'    (default TRUE)?
#' @param rm_footnotes logical, should semantic footer information be removed (default TRUE)?
#' @param rm_nodata_cols logical, should columns that have no alphanumeric data be removed (default TRUE)?
#' @param rm_escape a character vector that, if specified, is used to replace escape sequences in header
#'    and body cells (default ' ')
#' @param rm_invisible logical, should nodes that are not visible be removed (default TRUE)?
#' @param rm_whitespace logical, should leading/trailing whitespace be removed from cell values (default TRUE)?
#' @param colNames a character vector of column names, or a function that can be used to replace specific
#'    column names (default NULL)
#' @param ... additional arguments passed to HTML parsers
#'
#' @return An R data frame
#' @details The header formula has the following format: level1 + level2 + level3 + ... .
#' level1 specifies the main header dimension (column names). This information must
#' be for rows. level2 and deeper signify header dimensions that appear throughout the body.
#' Those information muste be for cell elements, not rows. Header information may be
#' one of the following types:
#'
#'\itemize{
#' \item the NULL value (default). No information passed, htmltab will try to identify
#' header elements through heuristics (heuristics only work for the main header)
#' \item A numeric vector that retrieves rows in the respective position
#' \item A character string of an XPath expression
#' \item A function that when evaluated produces a numeric or character vector
#' \item 0, when the process of finding the main header should be skipped (only works
#' for main header)
#' }
#'
#'
#' @author Christian Rubba <\url{http://www.christianrubba.com}>
#' @references \url{https://github.com/crubba/htmltab}
#' @examples
#' \dontrun{
#'# When no spans are present, htmltab produces output identical to XML's readHTMLTable()
#'
#'  url <- "http://en.wikipedia.org/wiki/World_population"
#'  xp <- "//caption[starts-with(text(),'World historical')]/ancestor::table"
#'  htmltab(doc = url, which = xp)
#'
#'  popFun <- function(node) {
#'    x <- XML::xmlValue(node)
#'    gsub(',', '', x)
#'  }
#'
#'  htmltab(doc = url, which = xp, bodyFun = popFun)
#'
#' #This table lacks header information. We provide them through colNames.
#' #We also need to set header = 0 to indicate that no header is present.
#' doc <- "http://en.wikipedia.org/wiki/FC_Bayern_Munich"
#' xp2 <- "//td[text() = 'Head coach']/ancestor::table"
#' htmltab(doc = doc, which = xp2, header = 0, encoding = "UTF-8", colNames = c("name", "role"))
#'
#' #htmltab recognizes column spans and produces a one-dimension vector of variable information,
#' #also removes automatically superscript information since these are usually not of use.
#'
#'  doc <- "http://en.wikipedia.org/wiki/Usage_share_of_web_browsers"
#'  xp3 <-  "//table[7]"
#'  bFun <- function(node) {
#'    x <- XML::xmlValue(node)
#'    gsub('%$', '', x)
#'  }
#'
#'  htmltab(doc = doc, which = xp3, bodyFun = bFun)
#'
#'
#' #When header information appear throughout the body, you can specify their
#' #position in the header formula
#'
#' htmltab("https://en.wikipedia.org/wiki/Arjen_Robben", which = 3,
#' header = 1:2 + "//tr/th[@@colspan='3' and not(contains(text(), 'Club'))]")
#' }

htmltab <- function(doc,
                    which = NULL,
                    header = NULL,
                    headerFun = function(node)XML::xmlValue(node),
                    headerSep = " >> ",
                    body = NULL,
                    bodyFun = function(node)XML::xmlValue(node),
                    complementary = T,
                    fillNA = NA,
                    rm_superscript = T,
                    rm_escape = " ",
                    rm_footnotes = T,
                    rm_nodata_cols = T,
                    rm_invisible = T,
                    rm_whitespace = T,
                    colNames = NULL,
                    ...){

  # Deparse
  header <- deparse(substitute(header), width.cutoff = 500L)
  body <- deparse(substitute(body), width.cutoff = 500L)
  ev_header <- eval_header(arg = header)
  ev_body <- eval_body(arg = body)

  # Check Inputs & Clean Up & Normalize tr --------
  table.Node <- check_type(doc = doc, which = which, ...)
  table.Node <- rm_nuisance(table.Node = table.Node,
                            rm_superscript = rm_superscript,
                            rm_footnotes = rm_footnotes,
                            rm_invisible = rm_invisible)
  table.Node <- normalize_tr(table.Node = table.Node)

  #Produce XPath for header and body and add class information
  LL <- identify_elements(table.Node = table.Node,
                          header = ev_header,
                          body = ev_body,
                          complementary = complementary)

  # Create Table ---------------------------

  #Retrieve Head Elements
  header.names <- make_header(trindex = LL$trindex$header,
                              table.Node = table.Node,
                              headerSep = headerSep,
                              headerFun = headerFun,
                              rm_escape = rm_escape,
                              rm_whitespace = rm_whitespace)


  # Create Body ---------------------------

  #Get Body Cell Nodes
  cells <- get_cells(table.Node = table.Node)
  #row_type <- get_row_type(table.Node = table.Node)


  #Extract and transform body cell elements
  vals <- get_cell_element(cells, elFun = bodyFun,
                           rm_escape = rm_escape,
                           rm_whitespace = rm_whitespace)

  #Produce rowspans and colspans lists from body cell
  body.rowspans <- get_span(cells, span = "rowspan")
  body.colspans <- get_span(cells, span = "colspan")

  #Produce table body
  tab <- expand_body(vals, colspans = body.colspans, rowspans = body.rowspans)


  # Finish ---------------------------

  #Produce DF
  tab <- make_colnames(df = tab,
                       header.names = header.names,
                       colNames = colNames,
                       header.xpath = LL$xpath$header)

  #Replace empty vals by NA
  tab[is.na(tab)] <- fillNA

  tab <- as.data.frame(tab, stringsAsFactors = F)

  # Inbody header
  tab <- create_inbody(tab = tab, table.Node = table.Node,
                       trindex = LL$trindex$inbody,
                       xpath = LL$xpath$inbody)

  # Subset
  tab <- tab[LL$trindex$body, ]

  #Check if there are no data columns
  if(isTRUE(rm_nodata_cols)){
    tab <- rm_empty_cols(df = tab)
  }

  return(tab)
}

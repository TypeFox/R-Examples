#' Retrieve table head rows
#'
#' @param table.Node the table node
#' @return the header element
#' @noRd
get_head <- function(trindex, table.Node) UseMethod()

get_head <- function(table.Node, trindex) {

  bracket <- paste(sprintf("@HTMLTABtrindex = '%s'", trindex), collapse = " or ")
  xpath <- sprintf("//tr[%s]", bracket)

  head_MAIN <- XML::xpathSApply(table.Node, path = xpath)

  return(head_MAIN)
}


#' Produce table header
#' @noRd
make_header <- function(trindex, table.Node, headerSep, headerFun,
                           rm_escape, rm_whitespace) UseMethod("make_header")

make_header.logical <- function(trindex, table.Node, headerSep, headerFun,
                                rm_escape, rm_whitespace){
  return(logical(0))
}

make_header.numeric <- function(trindex, table.Node, headerSep, headerFun,
                             rm_escape,
                             rm_whitespace){

  head <- get_head(table.Node = table.Node, trindex = trindex)

  header.colspans <- get_span(head, span = "colspan", tag = "td | th")
  header.rowspans <- get_span(head, span = "rowspan", tag = "td | th")

  header.names <- get_cell_element(head, tag = "td | th",
                                   elFun = headerFun,
                                   rm_escape = rm_escape,
                                   rm_whitespace = rm_whitespace)

  #Span head
  header.names <- span_header(header.names,
                              header.colspans,
                              header.rowspans,
                              headerSep = headerSep)

  return(header.names)
}


#' Creates header using span information
#'
#' @param header.names a list of header names
#' @param header.colspans a list of header colspans
#' @param header.rowspans a list of header rowspans
#' @param headerSep a character vector that is used as a seperator in the construction of the table's variable names (default value ' >> ')
#' @return a vector of header column names
#' @noRd
span_header <- function(header.names, header.colspans, header.rowspans, headerSep) {

  #has no header information
  if(length(header.names) == 0){
    header.name.table <- vector()
    return(header.name.table)
  }

  #Remove rows which have all empty cells
  empty.rows <- which(sapply(header.names, function(x) all(x == "")))
  if(!is_empty(empty.rows)){
    header.names <- header.names[-empty.rows]
    header.colspans <- header.colspans[-empty.rows]
    header.rowspans <- header.rowspans[-empty.rows]
  }

  #return empty header
  if(length(header.names) == 0){
    header.name.table <- vector()
    return(header.name.table)
  }

  header.name.table <- expand_header(header.names, header.colspans, header.rowspans)

  header.name.table <- lapply(header.name.table, function(col) col[!is.na(col)])

  header.name.table <- unlist(lapply(header.name.table, function(col) paste(col, collapse = headerSep)))

  return(header.name.table)
}

#' Expand the header
#' @noRd
expand_header <- function(vals, colspans, rowspans){

  body.table <- list()
  col <- 1
  n.row <- length(vals)

  repeat{

    # Break when there are no header information or when last column is missspecified
    if(is_empty(vals[[1]])){break} # || length(vals) > 1 && length(vals[[2]]) == 0

    body.row <- vector()

    for(row in 1:n.row) {

      length.col <- colspans[[row]][1]
      length.row <- rowspans[[row]][1]
      name <- vals[[row]][1]
      name <- ifelse(grepl("[[:alnum:][:punct:]]+", name), name, NA)

      if(is.na(length.col)) break

      #Remove cell info
      colspans[[row]] <- colspans[[row]][-1]
      rowspans[[row]] <- rowspans[[row]][-1]
      vals[[row]] <- vals[[row]][-1]

      #Expand along columnes
      colspans[[row]] <- append(colspans[[row]], rep(1, length.col - 1), 0)
      rowspans[[row]] <- append(rowspans[[row]], rep(length.row, length.col - 1), 0)
      vals[[row]] <- append(vals[[row]], rep(name, length.col - 1), 0)

      #Expand along rows
      length(vals[[1]])
      these.rows <- row:(row + (length.row - 1))
      rowspans[these.rows] <- lapply(rowspans[these.rows], append, 1, after = 0) #an erster Stelle
      colspans[these.rows] <- lapply(colspans[these.rows], append, 1, after = 0)
      vals[these.rows] <- lapply(vals[these.rows], append, NA, after = 0) #vals[these.rows] <- lapply(vals[these.rows], append, NA, after = 0)

      #remove the first colum info
      colspans[[row]] <- colspans[[row]][-1] #check for colspans different
      rowspans[[row]] <- rowspans[[row]][-1] #check for colspans different
      vals[[row]] <- vals[[row]][-1]

      #Add cell info to column name vector
      body.row <- c(body.row, name)
    }

    #Break if
    if(col > 1 && length(body.row) < length(body.table[[col - 1]])) break

    body.table[[col]] <- vector()
    body.table[[col]] <- body.row

    col <- col + 1
  }

  return(body.table)
}


#' Extracts cells elements
#'
#' @param cells a list of cell nodes
#' @param tag a character vector that provides information used in the XPath expression to extract the correct elements
#' @param elFun a function that is executed over the header/body cell nodes
#' @param rm_escape a character vector that, if specified, is used to replace escape sequences in header and body cells (default value ' ')
#' @param rm_whitespace logical, should leading/trailing whitespace be removed from cell values (
#'    default value TRUE)?
#' @return the body element
#' @noRd
get_cell_element <- function(cells, tag = "td | th", elFun, rm_escape, rm_whitespace) {

  cell.element <- lapply(cells, function(tr) {
    XML::xpathSApply(tr, tag, elFun)
  })

  if(!is.null(rm_escape)) {
    cell.element <- lapply(cell.element, function(el) gsub("([[:alpha:]])-[\b\n\t\r]([[:alpha:]])", "\\1\\2", el))
    cell.element <- lapply(cell.element, function(el) gsub("[\b \n \t \r]", rm_escape, el))
  }

  if(isTRUE(rm_whitespace)){
    cell.element <- lapply(cell.element, function(el) rm_str_white(el))
  }
  return(cell.element)
}


#' Extracts rowspan information
#'
#' @param cells a list of cell nodes
#' @param span a character for the span element name
#' @param tag a character vector that provides information used in the XPath expression to extract the correct elements
#' @return A list of row information from the cells
#' @noRd
get_span <- function(cells, span, tag = "td | th"){

  span.val <- lapply(cells, function(tr) {
    XML::xpathSApply(tr, tag, function(node) {
      val <- XML::xmlGetAttr(node, span)
      val <- ifelse(is.null(val) || val == "0", 1, val) #Check Firefox for colspan == 0
      val <- as.numeric(val)
      return(val)
    })
  })

  return(span.val)
}

#' Extracts header elements
#'
#' @param cells a list of cell nodes
#' @param tag a character vector that provides information used in the XPath expression to extract the correct elements
#' @return A list of header information from the cells
#' @noRd
get_header_elements <- function(cells, tag = "td | th"){

  header_elements <- lapply(cells, function(tr) {
    XML::xpathSApply(tr, tag, function(node) {
      if(XML::xmlName(node) != "sup") {
        value <- XML::xmlValue(node)
      }
      return(value)
    })
  })
  return(header_elements)
}

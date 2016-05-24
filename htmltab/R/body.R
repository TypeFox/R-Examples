
#' Creates body using span information
#'
#' @param vals a list of body cell values
#' @param colspans a list of body colspans
#' @param rowspans a list of body rowspans
#' @param header.names the header name vector
#' @return a matrix of the assembled body
#' @noRd
span_body <- function(vals, colspans, rowspans, header.names) {

  #Create empty container matrix
  if(length(header.names) == 0){
    n.cols <- sum(as.numeric(unlist(colspans[[1L]]))) #get idea from first row column length, should be right usually
  } else {
    n.cols <- length(header.names) #use header dim info
  }
  n.rows <- length(rowspans) #guess, better make flexible #lapply(rowspans, function(x) x[1]) %>% unlist %>% as.numeric %>% sum

  mat <- matrix(NA, ncol = n.cols, nrow = n.rows)
  #row.cont <- vector()

  col = 1L

  while(col <= n.cols){ # col start

    row <- 1L

    while(row <= n.rows){ # row start

      cel.val <- vals[[row]][col]
      col.span.length <- colspans[[row]][col]
      row.span.length <- rowspans[[row]][col]

      #This block controls for undefined col- and rowspans (hacky)
      if(is.na(row.span.length) && row < n.rows){
        col.span.length <- 1L
        row.span.length <- 1L
        cel.val <- mat[row, col -1L]
      }

      if(row.span.length < 2L) {
        mat[row, col] <- cel.val
      } else {

        if(row != n.rows){ #Control for situation: specified rows (main), last row demands col/rowspans
          index <- (row + 1L) : (row + row.span.length - 1L)

          for(counter in index){
            vals[[counter]] <- append(vals[[counter]], cel.val, (col - 1L )) #append col val is after
            rowspans[[counter]] <- append(rowspans[[counter]], 1L, (col - 1L))
            colspans[[counter]] <- append(colspans[[counter]], 1L, (col - 1L))
          }
          rowspans[[row]] <- rowspans[[row]][-col]
          rowspans[[row]] <- append(rowspans[[row]], 1L, (col-1L))

          mat[row:(row + row.span.length-1L), col] <- cel.val
        }
      }

      if(col.span.length > 1L){
        vals[[row]] <- append(vals[[row]], rep(cel.val, (col.span.length - 1L)), col) #append
        colspans[[row]] <- colspans[[row]][-col]
        colspans[[row]] <- append(colspans[[row]], rep(1, col.span.length), (col-1L))
        rowspans[[row]] <- rowspans[[row]][-col]
        rowspans[[row]] <- append(rowspans[[row]], rep(row.span.length, col.span.length), (col-1L))
      }

      row <- row + 1L
    } #row end

    col <- col + 1L
  } #col end

  return(mat)

} #function end



#' Expand the body
#' @noRd
expand_body <- function(vals, colspans, rowspans){

  body.table <- list()
  col <- 1
  n.row <- length(vals)

  repeat{

    #Break when there are no header information or when last column is missspecified
    if(is_empty(vals[[1]])){break} # || length(vals) > 1 && length(vals[[2]]) == 0

    body.row <- vector()

    for(row in 1:n.row) {

      length.col <- colspans[[row]][1]
      length.row <- rowspans[[row]][1]
      name <- vals[[row]][1]
      #name <- ifelse(name == "", NA, name)
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
      these.rows <- row:(row + (length.row - 1))
      rowspans[these.rows] <- lapply(rowspans[these.rows], append, 1, after = 0) #an erster Stelle
      colspans[these.rows] <- lapply(colspans[these.rows], append, 1, after = 0)
      vals[these.rows] <- lapply(vals[these.rows], append, name, after = 0) #vals[these.rows] <- lapply(vals[these.rows], append, NA, after = 0)

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

  # Cbind all
  tab <- do.call("cbind", body.table)

  return(tab)
}


#' Extract body cell values
#'
#' @param table.Node the table node
#' @return list of body information
#' @family get_head
#' @noRd
get_cells <- function(table.Node) {

  cells <- XML::xpathSApply(table.Node, path = "//tr")

  ifstop(is_empty(cells), sprintf("No body generated. Body is empty.
                 Try passing information to the body argument. Body XPath was '%s'.", body), call. = FALSE)

  return(cells)
}


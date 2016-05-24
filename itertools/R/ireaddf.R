ireaddf <- function(filenames, n, start=1, col.names, chunkSize=1000) {
  opencol <- function(fname) {
    # Extract the type of data from the file name
    m <- regexpr('(factor|character|integer|double)', fname)
    if (m < 1)
      stop('illegal file name: ', fname)
    type <- substr(fname, m, attr(m, 'match.length') + m - 1L)

    if (type == 'factor') {
      conn <- file(fname, 'rb')
      if (start > 1)
        seek(conn, where=4 * (start - 1))
      lfile <- sub('\\..+$', '.lev', fname)
      if (lfile == fname)
        lfile <- paste(fname, '.lev', sep='')
      lvls <- readLines(lfile)
    } else if (type == 'character') {
      conn <- file(fname, 'rt')
      if (start > 1)
        readLines(conn, n=start-1)
      lvls <- NULL
    } else if (type == 'integer') {
      conn <- file(fname, 'rb')
      if (start > 1)
        seek(conn, where=4 * (start - 1))
      lvls <- NULL
    } else if (type == 'double') {
      conn <- file(fname, 'rb')
      if (start > 1)
        seek(conn, where=8 * (start - 1))
      lvls <- NULL
    } else {
      stop('error: type = ', type)
    }

    list(type=type, conn=conn, lvls=lvls)
  }

  # Read the column data from the specified file
  readcol <- function(colinfo, n) {
    # Read the column data according to the data type
    if (colinfo$type == 'factor') {
      x <- readBin(colinfo$conn, what='integer', n=n)
      factor(colinfo$lvls[x], levels=colinfo$lvls)
    } else if (colinfo$type == 'character') {
      readLines(colinfo$conn, n=n)
    } else if (colinfo$type == 'integer') {
      readBin(colinfo$conn, what=colinfo$type, n=n)
    } else if (colinfo$type == 'double') {
      readBin(colinfo$conn, what=colinfo$type, n=n)
    } else {
      stop('internal error: type = ', colinfo$type)
    }
  }

  # Handle the col.names argument
  if (missing(col.names)) {
    col.names <- cnames(filenames)
  } else {
    if (any(duplicated(col.names)))
      stop('col.names must contain all unique names')
    if (length(col.names) != length(filenames))
      stop('col.names must be the same length as filenames')
  }

  # Open the column data files and return column information
  # in a list of (type, conn, lvls) lists
  names(filenames) <- col.names
  columndata <- lapply(filenames, opencol)
  stopped <- FALSE

  # Define the "next element" function needed for every iterator
  nextEl <- function() {
    # First check if we've already stopped
    if (stopped)
      stop('StopIteration')

    # Next check if there are any more rows to read
    if (n == 0) {
      for (colinfo in columndata)
        close(colinfo$conn)
      stopped <<- TRUE
      stop('StopIteration')
    }

    # Read the columns from the files into a list of equal length vectors
    df <- lapply(columndata, readcol, n=min(chunkSize, n))

    # Get the length of all columns
    xn <- sapply(df, length)
    mn <- min(xn)

    # Check if we ran out of date in at least one column data file
    if (mn == 0) {
      # Close all connections and throw a StopIteration
      for (colinfo in columndata)
        close(colinfo$conn)
      stopped <<- TRUE
      stop('StopIteration')
    }

    # Decrement n by the number of rows that we just read
    n <<- n - mn

    # Convert the list of vectors into a data frame if possible
    if (all(xn == mn)) {
      attr(df, 'row.names') <- .set_row_names(mn)
      class(df) <- 'data.frame'
    }

    df
  }

  # Construct and return the iterator object
  obj <- list(nextElem=nextEl)
  class(obj) <- c('abstractiter', 'iter')
  obj
}

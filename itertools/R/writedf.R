# Construct default column names from the specified file names
cnames <- function(filenames) {
  m <- regexpr('^(.+)_', filenames)
  if (any(m < 1)) {
    paste('X', seq_along(filenames), sep='.')
  } else {
    nms <- substr(filenames, m, attr(m, 'match.length') + m - 2L)
    if (any(duplicated(nms)))
      paste('X', seq_along(filenames), sep='.')
    else
      nms
  }
}

# This is a writedf combiner factory function
writedf.combiner <- function(filenames) {
  opencol <- function(i) {
    # Extract the type of data from the file name
    m <- regexpr('(factor|character|integer|double)', filenames[i])
    if (m < 1)
      stop('illegal file name: ', filenames[i])
    type <- substr(filenames[i], m, attr(m, 'match.length') + m - 1L)

    if (type == 'factor') {
      conn <- file(filenames[i], 'wb')
      lfile <- sub('\\..+$', '.lev', filenames[i])
      if (lfile == filenames[i])
        lfile <- paste(filenames[i], '.lev', sep='')
    } else if (type == 'character') {
      conn <- file(filenames[i], 'wt')
      lfile <- NULL
    } else if (type == 'integer') {
      conn <- file(filenames[i], 'wb')
      lfile <- NULL
    } else if (type == 'double') {
      conn <- file(filenames[i], 'wb')
      lfile <- NULL
    } else {
      stop('error: type = ', type)
    }

    list(type=type, conn=conn, lfile=lfile)
  }

  combine <- function(...) {
    # Write each column of "df" to a file
    dfs <- list(...)
    for (icol in seq_along(columndata)) {
      colinfo <- columndata[[icol]]

      if (colinfo$type == 'factor') {
        for (df in dfs) {
          if (! is.null(df)) {
            levsvar <- sprintf('levels.%d', icol)
            levs <- get(levsvar)
            if (is.null(levs)) {
              assign(levsvar, levels(df[[icol]]), inherits=TRUE)
            } else {
              if (any(levs != levels(df[[icol]]))) {
                # XXX should this be an error?
                # XXX should a try to fix the problem?
                warning('inconsistent levels found for column ', icol)
              }
            }
            writeBin(as.integer(df[[icol]]), colinfo$conn)
          }
        }
      } else if (colinfo$type == 'character') {
        for (df in dfs) {
          if (! is.null(df))
            writeLines(df[[icol]], colinfo$conn)
        }
      } else if (colinfo$type == 'integer') {
        for (df in dfs) {
          if (! is.null(df))
            writeBin(df[[icol]], colinfo$conn)
        }
      } else if (colinfo$type == 'double') {
        for (df in dfs) {
          if (! is.null(df))
            writeBin(df[[icol]], colinfo$conn)
        }
      } else {
        stop('unsupported type: ', colinfo$type)
      }
    }
    NULL
  }

  closeallcolumns <- function() {
    # Close all column data files
    for (i in seq_along(columndata)) {
      colinfo <- columndata[[i]]
      close(colinfo$conn)
      if (! is.null(colinfo$lfile)) {
        levsvar <- sprintf('levels.%d', i)
        levs <- get(levsvar)
        if (is.null(levs)) {
          warning(sprintf('not writing levels file for column %d\n', i))
        } else {
          writeLines(levs, colinfo$lfile)
        }
      }
    }
  }

  columndata <- lapply(seq_along(filenames), opencol)
  for (i in seq_along(filenames))
    assign(sprintf('levels.%d', i), NULL)

  # Construct and return the combiner object
  obj <- list(combine=combine, close=closeallcolumns)
  class(obj) <- c('combiner')
  obj
}

# Set plot area to export high resolution pictures
AQSysHR <- function (HR) {
  par.defaults <- par(no.readonly = TRUE)
  par(par.defaults)
  if (HR == TRUE) {
    clwd = 10
    par(
      mar = c(5,6,5,3) + 0.5,
      cex = 2.5,
      cex.lab = 2.5,
      cex.axis = 2.5,
      cex.main = 2.5,
      cex.sub = 2.5,
      lwd = clwd,
      mgp = c(4,2,0)
    )
  }else{
    clwd = 1
    par(
      mar = c(5,4,4,2) + 0.1,
      cex = 1,
      cex.lab = 1,
      cex.axis = 1,
      cex.main = 1,
      cex.sub = 1,
      lwd = clwd
    )
  }
  return(clwd)
}
# Display message when package is loaded
.onAttach <- function(libname, pkgname) {
  packageStartupMessage(
    '
    Be aware that LLSR is a collaborative package that still in
    development and your help is essential.\n
    If you found any bugs or have a suggestion, do not hesitate and
    contact us on https://github.com/eqipehub/LLSR/issues.\n
    You also can fork this project directly from github and commit
    improvements to us (https://github.com/eqipehub/LLSR).\n
    The information used in the database was obtained free of charge
    but it might be copyrighted by third parties and references must
    be included appropriately.\n
    Please use LLSR.info() to read more details about the current
    package version.', domain = NULL, appendLF = TRUE
  )
}
# merges data contained in all sheets from the same workbook which are named using
# the pattern datasource_sheetname_YX in its name.
# the order in the end is optional by now but will soon be part of the main function
AQSys.merge <- function(wrbk, sheets) {
  # initialize variables
  sys.nrow <- NULL
  sys.mrow <- NULL
  sys.ncol <- NULL
  sys.data <- NULL
  # make a loop through all sheets that satisfy the condition
  # the integrity of the worksheets will be evaluated before the merge
  # and determine which system in the workbook have bigger dataset (mrow)
  for (nSh in grep("datasource_", sheets)) {
    # determine the number of row and columns in the worksheet
    sys.nrow <- nrow(readWorksheet(wrbk, nSh, header = FALSE))
    sys.ncol <- ncol(readWorksheet(wrbk, nSh, header = FALSE))
    # initialize variables for first run only
    if (is.null(sys.mrow))
      sys.mrow <- sys.nrow
    if (sys.nrow > sys.mrow)
      sys.mrow <- sys.nrow
    # Each system must have two columns. If the total number in a sheet is odd
    #  it triggers an error (check AQSys.err.R for details)
    if (is.odd(sys.ncol))
      AQSys.err("2")
  }
  #
  for (nSh in grep("datasource_", sheets)) {
    #
    sys.temp <- readWorksheet(wrbk, nSh, header = FALSE)
    sys.nrow <- nrow(sys.temp)
    # populate all rows are initialized with NA
    if (sys.nrow < sys.mrow) {
      sys.temp[sys.mrow,] <- NA
    }
    # if no system has been added, just add the first sheet
    if (is.null(sys.data)) {
      sys.data <- sys.temp
      # but if sys.data have data, convert it to list and concatenate it with data
      # from the current sheet. Then convert it to dataframe and store it.
    } else {
      sys.data <-
        as.data.frame(c(sys.data, sys.temp), stringsAsFactors = FALSE)
    }
    
  }
  # return all data merged into a single dataframe
  invisible(sys.data)
}

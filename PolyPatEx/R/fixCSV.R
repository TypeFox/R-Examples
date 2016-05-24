
##  Function fixCSV
##
##' Tidies up a Comma Separated Value (CSV) file, ensuring that each
##' row of the table in the file contains the same number of commas,
##' and no empty rows are left below the table.
##'
##' \code{fixCSV} tidies up a Comma Separated Value (CSV) file
##' to ensure that the CSV file contains a strictly rectangular block
##' of data for input into R (ignoring any preliminary comment rows
##' via the \code{skip=} argument).
##'
##' CSV formatted files are a plain text file format for tabular data,
##' in which cell entries in the same row of a table are separated by
##' commas.  When such files are exported from other applications such
##' as spreadsheet software, the software has to decide whether any
##' empty cells to the right-hand side of, or below, the table or
##' spreadsheet should be represented by trailing commas in the CSV
##' file.  Such decisions can result in a \sQuote{ragged} table in the
##' CSV file, in which some rows contain fewer commas (\sQuote{short
##' rows}) or more commas (\sQuote{long rows}) than others, or where
##' empty rows below the table are included as comma-only rows in the
##' CSV file.
##'
##' While R's \code{\link{read.table}} and related functions can
##' sensibly extend short rows as needed, ragged tables in a CSV file
##' can still result in errors, unwanted empty rows (below the table)
##' or unwanted columns (to the right of the table) when the data is
##' loaded into R.
##'
##' \code{fixCSV} reads in a specified CSV file and removes or adds
##' commas to rows, to ensure that each row in the body of the table
##' contains the same number of cells as the header row of the table.
##' Any empty rows below the table are also removed.  The resulting
##' table is then written back to file, either to a new file with
##' \sQuote{FIXED} added to the filename (argument
##' \code{overwrite=FALSE}, the default) or overwriting the original
##' file (\code{overwrite=TRUE} - the original file is copied to a
##' \code{.BAK} file before being overwritten).
##'
##' Note that:
##' \itemize{
##'
##' \item The table of data in the CSV file \emph{must} contain a
##' header row of the correct length, since this row is used to
##' determine the correct number of columns for the table.  Note: if
##' this header row is too short, then subsequent rows will be
##' truncated to match the length of the header, so beware.
##' Misspecification of the \code{skip=} argument (see below) can
##' similarly lead to such corruption of the \sQuote{fixed} file.
##'
##' \item In the header row, any trailing commas representing empty
##' cells to the right of the (non-empty) header entries are first
##' removed before determining the correct number of columns for the
##' table.  Thus the length of the header row (and hence the assumed
##' width of the entire table) is determined by the \emph{right-most
##' non-empty cell} in the header row.
##'
##' \item \code{fixCSV} does not remove empty cells, rows or columns
##' within the interior (or on the left side) of the table - it is
##' concerned only with the right and bottom boundaries of the table.
##'
##' \item A \code{skip=} argument is included to tell \code{fixCSV} to
##' ignore the specified number of comment rows preceding the header
##' row.  Such rows are simply copied over into the output file
##' unchanged.  The default for this parameter is \code{skip=0}, so
##' that the first row in the data file is assumed to be the header
##' row. As noted above, misspecification of this argument can
##' seriously corrupt the output.
##'
##' \item \code{fixCSV} can overwrite your data file(s) (via
##' \code{overwrite=TRUE}), and althought it makes a backup of your
##' original file, you should still make sure that you have a separate
##' backup of your data file in a safe place before using this
##' function!  The author of this code takes no responsibility for any
##' data loss or corruption as a result of the use of this routine...
##'
##' }
##'
##' @title Tidy a comma separated value (CSV) file
##' @param file character: the name of the CSV file to be
##' \sQuote{fixed}.
##' @param skip integer: the number of lines in the CSV file to skip
##' before the header row of the table.  The skipped lines are copied
##' directly to the output file unchanged.  The default is
##' \code{skip=0}, implying that the header row is the first row of
##' the CSV file.
##' @param overwrite logical: Write output to a separate,
##' \sQuote{FIXED} file (\code{overwrite=FALSE}, the default), or
##' overwrite the original file (\code{overwrite=TRUE})?  If
##' \code{overwrite=TRUE}, the original file is copied to a
##' \code{.BAK} file before being overwritten.
##' @author Alexander Zwart (alec.zwart at csiro.au)
##' @export
##' @examples
##' \dontrun{
##'
##' ## Assuming CSV file 'alleleDataFile.csv' exists in the current
##' ## directory.  The following overwrites the CSV file - make sure
##' ## you have a backup!
##'
##' fixCSV("alleleDataFile.csv",overwrite=TRUE)
##'
##' }
##'
fixCSV <- function(file, skip=0, overwrite=FALSE) {
  ##
  foundShortLines <- FALSE
  foundLongLines <- FALSE
  foundEmptyTrailingLines <- FALSE
  ##
  vv <- readLines(file)
  if (skip > 0) { ## skip the preamble, if any
    preamble <- vv[1:skip]
    vv <- vv[-(1:skip)]
  }
  ##
  ##Strip any trailing delimiters or whitespace in header
  vv[1] <- sub("[ ,]+$","",vv[1])
  ##
  ##strsplit() neglects any trailing empty strings (see ?strsplit
  ## regarding matches and the beginning and end of a string).  I find
  ## that this complicates matters, unless I add an extra delimiter
  ## to the end of each string in vv.  Inelegant, perhaps, but the
  ## simplest solution...
  vv[-1] <- paste(vv[-1],",",sep="")
  ##
  ## Detect trailing empty rows (spaces and commas only) and remove.
  emptyLines <- grep("^[ ,]+$",vv)
  if ( (length(emptyLines) > 0) && (max(emptyLines)==length(vv)) ) {
    ##Found at least one empty line at END of the vector
    foundEmptyTrailingLines <- TRUE
    ##A trick courtesy of D Lovell :  Generate a sequence of length
    ## length(emptyLines), working BACK from length(vv), then
    ## reverse it:
    ss <- rev(seq(from=length(vv), by=-1,
                  along.with=emptyLines))
    ##Now, emptyLines entries that match the sequence ss are the
    ## trailing empty rows, so remove 'em...
    vv <- vv[-emptyLines[emptyLines == ss]]
  }
  ##
  ##Now check (and fix) the lengths of the remaining rows
  ss <- strsplit(vv,split=",")
  ll <- sapply(ss,length)
  headerLength <- ll[1]
  for (line in 2:length(ss)) {
    if (ll[line] < headerLength) {
      foundShortLines <- TRUE
      ss[[line]] <- c(ss[[line]],rep("",headerLength-ll[line]))
    } else if (ll[line] > headerLength) {
      foundLongLines <- TRUE
      ##Strip out the extra entries
      ss[[line]] <- ss[[line]][1:headerLength]
    }
  }
  ##Convert back to vector
  vv <- sapply(ss,
               function(tvec) {
                 paste(tvec,collapse=",")
               })
  ##
  if (foundShortLines || foundLongLines || foundEmptyTrailingLines) {
    ## Add back the preamble, if any:
    if (skip>0) {
      vv <- c(preamble,vv)
    }
    if (foundShortLines) {
      cat("\n\n Note : FixCSV found short rows in",
          file,"\n")
      cat("        These have been extended to match the header\n")
    }
    if (foundLongLines) {
      cat("\n\n Note : FixCSV found long rows in",
          file,"\n")
      cat("        The extra entries have been stripped\n")
    }
    if (foundEmptyTrailingLines) {
      cat("\n\n Note : FixCSV found empty trailing rows in",
          file,"\n")
      cat("        The extra rows have been removed\n")
    }
    rr <- regexpr(".([^.]+)$",file)
    ##Alternative, from Regexp cookbook - may be more robust:
    ## regexpr('.[^.\\/:*?"<>|\r\n]+$',file)
    ##Had to drop a leading '\' to get the cookbook version to
    ## work...
    fname <- substr(file,1,rr-1)
    extn <- substr(file,rr+1,nchar(file))
    if (overwrite) {  ##Back up the original
      file.copy(from=file,
                to=paste(fname,".BAK",sep=""),
                overwrite=TRUE)
    } else {
      if(rr == 1) {
        file <- paste(file,"FIXED",sep="")
      } else {
        file <- paste(fname,"FIXED.",extn,sep="")
      }
    }
    writeLines(vv,file)
    cat("\n        Output written to",file,"\n\n")
  } else {
    cat("\n No problems found in",file,"\n\n")
  }
  return(invisible(NULL))
}

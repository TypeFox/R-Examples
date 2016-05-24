#' Parse input files or folders with immune receptor repertoire data.
#'
#' @aliases repLoad 
#'
#' @description
#' Load the immune receptor repertoire data from the given input: either a file name, a list of file names, a name of the folder with repertoire files,
#' or a list of folders with repertoire files. The folder / folders must contain only files with the specified format.
#' Input files could be either text files or archived with gzip ("filename.txt.gz") or bzip2 ("filename.txt.bz2").
#' For a general parser of table files with cloneset data see \code{\link{parse.cloneset}}.
#' 
#' Parsers are available for:
#' MiTCR ("mitcr"), MiTCR w/ UMIs ("mitcrbc"), MiGEC ("migec"), VDJtools ("vdjtools"), 
#' ImmunoSEQ ("immunoseq" or 'immunoseq2' for old and new formats respectively),
#' MiXCR ("mixcr"), IMSEQ ("imseq") and tcR ("tcr", data frames saved with the `repSave()` function).
#' 
#' Output of MiXCR should contain either all hits or best hits for each gene segment.
#' 
#' Output of IMSEQ should be generated with parameter "-on". In this case there will be no positions of aligned gene segments in the output data frame
#' due to restrictions of IMSEQ output.
#' 
#' tcR's data frames should be saved with the `repSave()` function.
#' 
#' For details on the tcR data frame format see \link{parse.file}.
#' 
#' @param .path Character vector with path to files and / or folders.
#' @param .format String that specifies the input format.
#' 
#' @seealso \link{parse.file}
#' 
#' @examples 
#' \dontrun{
#' datalist <- repLoad(c("file1.txt", "folder_with_files1", "another_folder"), "mixcr")
#' }
repLoad <- function (.path, .format = c("mitcr", "migec")) {
  res <- list()
  
  for (i in 1:length(.path)) {
    if (dir.exists(.path[i])) {
      res <- c(res, parse.folder(.path[i], .format[1]))
    } else if (file.exists(.path[i])) {
      res <- c(res, list(parse.file(.path[1], .format[1])))
    } else {
      cat('Can\'t find folder or file:\t"', .path[i], '"', sep = '', end = '\n')
    }
  }
  
  res
}


#' Save tcR data frames to disk as text files or gzipped text files.
#' 
#' @description
#' Save repertoire files to either text files or gzipped text files.
#' You can read them later by \code{repLoad} function with \code{.format = "tcr"}.
#' 
#' @param .data Either tcR data frame or a list of tcR data frames.
#' @param .format "txt" for simple tab-delimited text tables, "gz" for compressed (gzipped) tables.
#' @param .names Names of output files. By default it's an empty string so names will be taken from names of the input list.
#' @param .folder Path to the folder with output files.
#' 
#' @seealso \link{repLoad}
repSave <- function (.data, .format = c("txt", "gz"), .names = "", .folder = "./") {
  if (has.class(.data, 'data.frame')) { .data <- list(Sample = .data) }
  
  .folder <- paste0(.folder, "/")
  
  postfix <- ".txt"
  filefun <- function (...) file(...)
  if (.format[1] == "gz") { 
    postfix <- ".txt.gz"
    filefun <- function (...) gzfile(...)
  }
  
  if (.names[1] == "") {
    .names = paste0(.folder, names(.data), postfix)
  } else {
    if (length(.data) != length(.names)) {
      cat("Number of input data frames isn't equal to number of names\n")
      return(NULL)
    } else {
      .names = paste0(.folder, .names, postfix)
    }
  }
  
  for (i in 1:length(.data)) {
    cat("Writing", .names[i], "file...\t")
    fc <- filefun(description = .names[i], open = "w")
    write.table(.data[[i]], fc, quote = F, row.names = F, sep = '\t')
    close(fc)
    cat("Done.\n")
  }
}
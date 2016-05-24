#' Read and parse .ini file to list
#'
#' @param filepath file to parse
#' @param encoding Encoding of filepath parameter, will default to system
#' encoding if not specifield
#'
#' @details Lines starting with '#' or ';' are comments and will not be parsed
#'
#' @seealso \code{\link{write.ini}}
#'
#' @return List with length equivalent to number of [sections], each section is
#' a new list
#'
#' @examples
#' ## Create a new temp ini for reading
#' iniFile <- tempfile(fileext = '.ini')
#'
#' sink(iniFile)
#' cat("; This line is a comment\n")
#' cat("# This one too!\n")
#' cat("[    Hello World]\n")
#' cat("Foo = Bar          \n")
#' sink()
#'
#' ## Read ini
#' checkini <- read.ini(iniFile)
#'
#' ## Check structure
#' checkini
#' checkini$`Hello World`$Foo
#'
#' @export
#'
read.ini <- function(filepath, encoding = getOption("encoding")) {

  sectionREGEXP <- '^\\s*\\[\\s*(.+?)\\s*]'
  # match section and capture section name

  keyValueREGEXP <- '^\\s*[^=]+=.+'
  # match "key = value" pattern

  ignoreREGEXP <- '^\\s*[;#]'
  # match lines with ; or # at start

  # amazing lack of trim at old versions of R
  trim <- function(x) sub('^\\s*(.*?)\\s*$', '\\1', x)

  ini <- list()
  con <- file(filepath, open = 'r', encoding = encoding)
  on.exit(close(con))

  while ( TRUE ) {

    line <- readLines(con, n = 1, encoding = encoding, warn = F)
    if ( length(line) == 0 ) {
      break
    }

    if ( grepl(ignoreREGEXP, line) ) {
      next
    }

    if ( grepl(sectionREGEXP, line) ) {
      matches <- regexec(sectionREGEXP, line)
      lastSection <- regmatches(line, matches)[[1]][2]
    }

    if ( grepl(keyValueREGEXP, line) ) {
      tempKeyValue <- trim(unlist(strsplit(line, "=")))
      key <- tempKeyValue[1]
      value <- tempKeyValue[2]

      ini[[ lastSection ]] <- c(ini[[ lastSection ]], list(key = value))
      names(ini[[ lastSection ]])[ match('key', names(ini[[ lastSection ]])) ] <- key
    }

  }

  ini
}

#' Write list to .ini file
#'
#' @param x List with structure to be write at .ini file.
#'
#' @param filepath file to write
#' @param encoding Encoding of filepath parameter, will default to system
#' encoding if not specifield
#'
#' @seealso \code{\link{read.ini}}
#'
#' @examples
#' ## Create a new temp ini for writing
#' iniFile <- tempfile(fileext = '.ini')
#'
#' ## Create a new list holding our INI
#' newini <- list()
#' newini[[ "Hello World" ]] <- list(Foo = 'Bar')
#'
#' ## Write structure to file
#' write.ini(newini, iniFile)
#'
#' ## Check file content
#' \dontrun{
#' file.show(iniFile)
#' }
#'
#' @export
#'
write.ini <- function(x, filepath, encoding = getOption("encoding")) {

  con <- file(filepath, open = 'w', encoding = encoding)
  on.exit(close(con))

  for(section in names(x) ) {
    writeLines( paste0('[', section, ']'), con)
    for (key in x[ section ]) {
      writeLines( paste0(names(key), ' = ', key), con)
    }
    writeLines("", con)
  }

}

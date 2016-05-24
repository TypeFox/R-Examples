#' Convert an ISI or SCOPUS Export file into a data frame
#'
#' It converts a SCOPUS and Thomson Reuters' ISI Web of Knowledge export file and create a data frame from it, with cases corresponding to articles and variables to Field Tag in the original file.
#'
#' Actually the function allows to convert both SCOPUS/ISI files in bibtext format and just ISI files in plain text format.
#'
#' @param file is a character array containing data read from an ISI WoK Export file (in plain text or bibtex format) or SCOPUS Export file (exclusively in bibtex format).
#' @param dbsource is a character indicating the bibliographic database. \code{dbsource} can be \code{"isi"} or \code{"scopus"}. Default is \code{dbsource = "isi"}.
#' @param format is a character indicating the format of the SCOPUS and Thomson Reuters' ISI Web of Knowledge export file. \code{format} can be \code{"bibtex"} or \code{"plaintext"}. Default is \code{format = "bibtex"}.
#' @return a data frame with cases corresponding to articles and variables to Field Tag in the original export file.
#'
#' data frame columns are named using the standard ISI WoS Field Tag codify. The main field tags are:
#'
#' \tabular{lll}{
#' \code{AU}\tab   \tab Authors\cr
#' \code{TI}\tab   \tab Document Title\cr
#' \code{SO}\tab   \tab Publication Name (or Source)\cr
#' \code{JI}\tab   \tab ISO Source Abbreviation\cr
#' \code{DT}\tab   \tab Document Type\cr
#' \code{DE}\tab   \tab Authors' Keywords\cr
#' \code{ID}\tab   \tab Keywords associated by SCOPUS or ISI database \cr
#' \code{AB}\tab   \tab Abstract\cr
#' \code{C1}\tab   \tab Author Address\cr
#' \code{RP}\tab   \tab Reprint Address\cr
#' \code{CR}\tab   \tab Cited References\cr
#' \code{TC}\tab   \tab Times Cited\cr
#' \code{PY}\tab   \tab Year\cr
#' \code{SC}\tab   \tab Subject Category\cr
#' \code{UT}\tab   \tab Unique Article Identifier\cr
#' \code{DB}\tab   \tab Database\cr}
#'
#' for a complete list of filed tags see: \href{https://images.webofknowledge.com/WOK46/help/WOS/h_fieldtags.html}{ISI WoS Field Tags}
#' @examples
#' # An ISI or SCOPUS Export file can be read using \code{\link{readLines}} function:
#'
#' # largechar <- readLines('filename.txt')
#'
#' # filename.txt is an ISI or SCOPUS Export file in plain text or bibtex format.
#' # The file have to be saved without Byte order mark (U+FEFF) at the beginning
#' # and EoF code at the end of file.
#' # The original file (exported by ISI or SCOPUS search web site) can be modified
#' # using an advanced text editor like Notepad++ or Emacs.
#'
#' #  biblio <- readLines('~/extdata/bibliometrics_articles.txt')
#'
#' data(biblio)
#'
#' biblio_df_df <- convert2df(file = biblio, dbsource = "isi", format = "bibtex")
#'
#' @seealso \code{\link{scopus2df}} for converting SCOPUS Export file (in bibtex format)
#' @seealso \code{\link{isibib2df}} for converting ISI Export file (in bibtex format)
#' @seealso \code{\link{isi2df}} for converting ISI Export file (in plain text format)
#' @family converting functions

convert2df<-function(file,dbsource="isi",format="bibtex"){

  if (length(setdiff(dbsource,c("isi","scopus")))>0){
    cat("\n 'dbsource' argument is not properly specified")
    cat("\n 'dbsource' argument has to be a character string matching 'isi or 'scopus'.\n")}
  if (length(setdiff(format,c("plaintext","bibtex")))>0){
    cat("\n 'format' argument is not properly specified")
    cat("\n 'format' argument has to be a character string matching 'plaintext or 'bibtex'.\n")}

  file=iconv(file, "latin1", "ASCII", sub="")
  switch(dbsource,
    isi={
      switch(format,
             bibtex={M=isibib2df(file)},
             plaintext={M=isi2df(file)}
      )},
    scopus={M=scopus2df(file)
    }
)

  return(M)

}

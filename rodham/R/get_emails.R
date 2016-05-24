#' Get emails and its contents
#'
#' @description Get the content of Hillary Rodham Clinton's emails by release.
#'
#' @param release Name of the batch of release of emails; see details.
#' @param save.dir Directory where to save the extracted text defaults to
#' \code{getwd()}
#' @param extractor Full path to pdf extractor (text to pdf), see details.
#'
#' @details Below are the valid values for \code{release}; follows the
#' \href{http://graphics.wsj.com/hillary-clinton-email-documents/}{WSJ} naming
#' convention.
#' \itemize{
#' \item Benghazi
#' \item June
#' \item July
#' \item August
#' \item September
#' \item October
#' \item November
#' \item January 7
#' \item January 29
#' \item February 19
#' \item february 29
#' \item December
#' \item Non-disclosure
#' }
#'
#' The \code{extractor} argument is the full path to your pdftotext.exe
#' extractor; visit \href{http://www.foolabs.com/xpdf/download.html}{xpdf} to
#' download or try \code{\link{get_xpdf}} which attempts to download and
#' unzip the text to pdf extractor. See examples.
#'
#' @return Fetches email zip file from the WSJ and extract text files in
#' \code{save.dir}, returns full path to directory that contains parsed txt
#' files.
#'
#' @examples
#' \dontrun{
#' # get xpdf extractor
#' ext <- get_xpdf()
#'
#' # create
#' dir.create("./emails")
#'
#' # get emails released in august
#' emails_aug <- get_emails(release = "August", save.dir = "./emails",
#'                      extractor = ext)
#'
#' # use manually downloaded extractor
#' ext <- "C:/xpdfbin-win-3.04/bin64/pdftotext.exe"
#'
#' # get emails related to Benghazi released in December
#' emails_bengh <- get_emails(release = "Benghazi", extractor = ext,
#'                            save.dir = "./emails")
#'
#' files <- list.files(emails_bengh)
#' content <- lapply(1:length(files), function(x){
#'    readLines(paste0(emails_bengh, "/", files[[x]]))
#' })
#' }
#'
#' @seealso \code{\link{get_xpdf}}
#'
#' @author John Coene \email{jcoenep@gmail.com}
#'
#' @export
get_emails <- function(release, save.dir = getwd(), extractor){
  if (missing(extractor)) {
    stop("Missing extractor, see get_extractor")
  }
  if (!dir.exists(save.dir)) {
    stop("save.dir does not exist")
  }
  if(!file.exists(extractor)){
    stop("wrong path to extractor")
  }
  v <- suppressWarnings(system(paste(extractor, "-v"))) # test extractor
  if (!v %in% c(99, 0, 1, 2, 3)) {
    stop("incorrect path to extractor, see get_xpdf")
  }
  if (!dir.exists(save.dir)) {
    stop("save.dir does not exist")
  }
  uri <- checkRelease(release) # check release input and return URL
  temp_zip <- tempfile(fileext = ".zip") # create temp
  download.file(uri, destfile = temp_zip) # download
  temp_dir <- tempdir() # create dir
  unzip(zipfile = temp_zip, exdir = temp_dir) # unzip
  unlink("temp_zip", recursive = TRUE) # delete temp once unzipped
  files <- list.files(temp_dir)  # list files
  files <- files[grep("pdf", files)] # only take pdf files
  dest <- gsub(".pdf", ".txt", files) # name destinations for extraction
  cat("Extracting content from", length(files), "pdf files...\n")
  pb <- txtProgressBar(style = 3)
  for (i in 1:length(files)) {
    pdf <- paste0(temp_dir,"\\", files[i])
    txt <- paste0(save.dir,"/", dest[i])
    system(paste(extractor, "-nopgbrk" , pdf, txt, sep = " "),
           wait = TRUE)
    setTxtProgressBar(pb, i/length(files))
  }
  unlink("temp_dir", recursive = TRUE) # delete temp once extracted
  close(pb)
  message("emails .txt files saved at: ", save.dir)
  return(save.dir)
}

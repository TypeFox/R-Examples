#' Extract text from PDF files and return a word-occurrence data.frame.
#'
#' \code{getPDF} returns a word-occurrence data.frame from PDF files.
#' It needs \code{XPDF} in order to run (http://www.foolabs.com/xpdf/download.html),
#' and uses \code{parallel} to perform parallel computation.
#'
#' @param myPDFs A character vector containing PDF file names.
#' @param minword An integer specifying the minimum number of letters per word
#'   into the returned data.frame.
#' @param maxword An integer to specifying the maximum number of letters per
#'   word into the returned data.frame.
#' @param minFreqWord An integer specifying the minimum word frequency into the
#'   returned data.frame.
#' @param pathToPdftotext A character containing an alternative path to XPDF
#'   \code{pdftotext} function, see Details section.
#' @details \code{getPDF} uses \code{XPDF pdftotext} function to extract the
#'   content of PDF files into a TXT file. If  \code{pdftotext} is not in the
#'   \code{PATH}, an alternative is to provide the full path of the program into
#'   the \code{pathToPdftotext} parameter.
#' @return A list of list with word-occurrence data.frame and file name.
#' @examples
#' \dontrun{
#' getPDF(myPDFs = "mypdf.pdf")
#' }
#' @export
getPDF <- function(myPDFs, minword = 1, maxword = 20, minFreqWord = 1, pathToPdftotext = ""){
  ncores <- parallel::detectCores()
  if(length(myPDFs)<ncores){ncores <- length(myPDFs)}
  cl <- parallel::makeCluster(ncores)
  # parallel::clusterExport(cl = cl, varlist = c("minword","maxword","minFreqWord"))  ### for testing purposes
  parallel::clusterExport(cl = cl, varlist = c("pathToPdftotext", "preProcTxt", "postProcTxt"), envir = environment())
  parallel::clusterEvalQ(cl, {library(tm); library(SnowballC)})
  on.exit(parallel::stopCluster(cl))
  d <- parallel::parLapply(cl, myPDFs, function(i){ # lapply to all PDF files
    if (Sys.info()[1]=="Windows"){ # extract txt from PDF with Windows
      tryCatch(system(paste0("\"pdftotext\" \"", i, "\""), wait = TRUE),
        error = system(paste("\"", pathToPdftotext, "\" \"", i, "\"", sep = ""), wait = TRUE))
      # system(paste("\"", pathToPdftotext, "\" \"", i, "\"", sep = ""), wait = TRUE)
    }else{
      if (Sys.info()[1]=="Linux"){
        system(paste("pdftotext" , i, sep = " "), wait = TRUE) # extract txt from PDF with Linux
      }else{ ### MacOS = "/"
        system(paste("pdftotext" , i, sep = " "), wait = TRUE) # extract txt from PDF with Mac OSX
      }
    }
    filetxt <- paste(strsplit(i, split = "\\.")[[1]][1], ".txt", sep = "")
    txt <- preProcTxt(filetxt)
    if (length(txt)<5){
      print(paste0("File ", filetxt, " does not contain text and will be ignored."))
      txt <- NULL
      file.remove(filetxt)
    }
    d1 <- postProcTxt(txt = txt, minword = minword, maxword = maxword, minFreqWord = minFreqWord)
    d2 <- strsplit(filetxt, split = "\\.")[[1]][1]
    d <- list(wc = d1, name = d2)
    if(nrow(d[[1]])==0){d <- NULL}
    file.remove(filetxt) # remove file
    return(d)
  })
  d <- d[vapply(d, Negate(is.null), NA)]
  return(d)
}

#' Extract text from TXT files and return a word-occurrence data.frame.
#'
#' @param myTXTs A character vector containing TXT file names (or complete path
#'   to these files).
#' @return A list of list with word-occurrence data.frame and file name.
#' @examples
#' data("loremIpsum")
#' loremIpsum01 <- loremIpsum[1:100]
#' loremIpsum02 <- loremIpsum[101:200]
#' loremIpsum03 <- loremIpsum[201:300]
#' loremIpsum04 <- loremIpsum[301:400]
#' loremIpsum05 <- loremIpsum[401:500]
#' subDir <- "RESULTS"
#' dir.create(file.path(getwd(), subDir), showWarnings = FALSE)
#' write(x = loremIpsum01, file = "RESULTS/loremIpsum01.txt")
#' write(x = loremIpsum02, file = "RESULTS/loremIpsum02.txt")
#' write(x = loremIpsum03, file = "RESULTS/loremIpsum03.txt")
#' write(x = loremIpsum04, file = "RESULTS/loremIpsum04.txt")
#' write(x = loremIpsum05, file = "RESULTS/loremIpsum05.txt")
#' wordOccuFreq <- getTXT(myTXTs = list.files(path = paste0(getwd(), 
#'   "/RESULTS/"), pattern = "loremIpsum", full.names = TRUE))
#' file.remove(list.files(pattern = "loremIpsum"))
#' @export
getTXT <- function(myTXTs){
	lapply(myTXTs, function(i){
		txt <- preProcTxt(i)
		d <- postProcTxt(txt = txt)
		df <- list(wc = d, name = i)
		return(df)
	})
}

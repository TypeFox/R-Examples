#' Download attachment files of a dataset in BEFdata portal
#'
#' This function will download attachment files associated with a dataset into specified directory.
#'
#' @param dataset The ID of dataset you like to fetch the attachments from. You find the ID on the dataset page
#'	  on the BEFdata portal.
#' @param directory The directory to store attachment files to. By default it creates a folder called "downloads"
#' 	  under the current working directory. The default can be changed by bef.options.
#' @param curl If using in a loop, call getCurlHandle() first and pass the returned value
#'         in here (avoids unnecessary footprint).
#' @param \dots  arguments passed to \code{\link[RCurl]{getURLContent}}.
#' @return A data frame of file information is returned invisibly. NULL is returned when
#'         the dataset has no attachement files.
#' @export bef.portal.get.attachments bef.get.attachments bef.get.attachments_for bef.portal.get.attachments_for
#' @aliases bef.get.attachments bef.get.attachments_for bef.portal.get.attachments_for

bef.portal.get.attachments <- bef.get.attachments <- bef.get.attachments_for <- bef.portal.get.attachments_for <- function(dataset, directory = bef.options('download_dir'), curl = getCurlHandle(), ...) {
  dataset_url = dataset_url(dataset, "freeformat", user_credentials = bef.options('user_credentials'))
  freeformats_csv = getURLContent(dataset_url, curl = curl, ...)
  if (getCurlInfo(curl)$response.code != 200) {
    stop("Dataset not found or not accessible. Please check your credentials and make sure you have access to it")
  }
  files = read.csv(text = freeformats_csv, stringsAsFactors=F)
  if (nrow(files)) {
    if (!file.exists(directory)) dir.create(directory)
    files$path = file.path(directory, sapply(files$Filename, suggest_filename, dir = directory, USE.NAMES=T))
    for (i in seq_len(nrow(files))) {
      f = CFILE(files$path[i], mode="wb")
      cat(sprintf("Saving %s => %s\n", sQuote(files$Filename[i]), sQuote(files$path[i])))
      flush.console()
      curlPerform(url = files$URL[i], writedata = f@ref)
      close(f)
    }
    invisible(files)
  } else {
    warning(paste("No attachement files for dataset:", dataset))
    invisible(NULL)
  }
}

## For package updates information

#' @importFrom utils packageDescription
hyfoUpdates <- function(){
  page <- readLines('http://yuanchao-xu.github.io/hyfo/')
  updatesLine <- grep('id=\\"updates"', page)
  versionLine <- updatesLine + 2
  
  version <- unlist(strsplit(page[versionLine], split = ' '))[2]
  version_local <- packageDescription("hyfo")$Version
  
  
  # the first tow digit is the most important part of the version
  version12 <- unlist(strsplit(version, split = "[.]"))[1:2]
  version_local12 <- unlist(strsplit(version_local, split = "[.]"))[1:2]
  
  sameVersion <- version12 == version_local12
  
  if (any(sameVersion == FALSE)) {
    # generate message
    version_msg <- strsplit(strsplit(page[versionLine], split = '<p>')[[1]][2], split = '</p>')[[1]]
    infoLine_start <- versionLine + 2
    infoLine_end <- grep('<p>For historical releases and the introduction of updates about each version', page) - 1
    info_msg <- character()
    for (infoLine in infoLine_start:infoLine_end) {
      info_line <- strsplit(strsplit(page[infoLine], split = '>')[[1]][2], split = '<')[[1]][1]
      if (!is.na(info_line)) info_msg <- c(info_msg, info_line)
    }
    
    install_msg <- 'More details on http://yuanchao-xu.github.io/hyfo/'
    
    message_out <- paste(version_msg, paste(info_msg, collapse = '\n'), install_msg, sep = '\n')
  } else message_out <- NULL
  return(message_out)
}

.onAttach <- function(libname, pkgname) {
  message_out <- suppressWarnings(try(hyfoUpdates(), silent = TRUE))
  if (!is.null(message_out)) {
    if (grepl('Version', message_out)) {
      packageStartupMessage(message_out)
    }
  }
}

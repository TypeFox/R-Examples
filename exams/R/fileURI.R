fileURI <- function(file) {
  ## see mime types at e.g.
  ## http://www.freeformatter.com/mime-types-list.html
  f_ext <- tolower(file_ext(file))
  if(f_ext %in% c("bmp", "png", "jpg", "jpeg", "gif",
    "csv", "raw", "txt", "xls", "xlsx", "zip", "pdf", "doc", "docx",
    "rda", "dta")) {
    mime <- switch(file_ext(file),
      "bmp" = "image/bmp",
      "png" = "image/png",
      "jpg" = "image/jpeg",
      "jpeg" = "image/jpeg",
      "gif" = "image/gif",
      "csv" = "text/csv",
      "raw" = "text/plain",
      "txt" = "text/plain",
      "xls" = "application/vnd.ms-excel",
      "xlsx" = "application/vnd.openxmlformats-officedocument.spreadsheetml.sheet",
      "zip" = "application/zip",
      "pdf" = "application/pdf",
      "doc" = "application/msword",
      "docx" = "application/vnd.openxmlformats-officedocument.wordprocessingml.document",
      "rda" = "application/octet-stream",
      "dta" = "application/octet-stream",
    )
    rval <- base64enc::dataURI(file = file, mime = mime)
  } else {
    owd <- getwd()
    setwd(dirname(file))
    zip(zipfile = zipname <- paste(file_path_sans_ext(basename(file)), "zip", sep = "."),
      files = basename(file))
    rval <- base64enc::dataURI(file = zipname, mime = "application/zip")
    setwd(owd)
  }
  rval
}

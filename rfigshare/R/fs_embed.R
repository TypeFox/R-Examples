
#' get image url
#' 
#' get image url
#' @param id a (public) figshare figure id number 
#' @param debug logical, enable debugging mode?
#' @return a url to the image file
#' @details this is currently an unstable hack of html parsing
#' @import XML httr
#' @export 
fs_image_url <- function(id, debug = FALSE) {
  a <- fs_details(id)
  b <- GET(a$doi)
  if(debug | b$status_code != 200)
    b
  else {
    doc <- htmlParse(content(b, as = "text"))
    path <- xpathSApply(doc, "//div[@class='filesdownload' and @id='download_all']/a/@href")[[1]]
#  resp <- GET(path)
#  if(resp$headers$`content-type` != "image/png")
    path
  }
}


#' Upload a figure to figshare and return the url
#' 
#' @param file path to an image file
#' @return a url to the image file
#' @details use with opts_knit$set(upload.fn = fs_embed)
#' @export
fs_embed <- function(file) {
  ## Read in title, and details from an environment?  from chunk? 
  title <- file 
  description <- "embedded file automatically uploaded to figshare from R"
  id <- fs_new_article(title, description, type = "figure", visibility = "public")
  fs_image_url(id)
}
#
#

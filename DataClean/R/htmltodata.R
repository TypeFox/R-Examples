#' "htmltodata" function is used to transfer information from html files to R or xlsx files
#'
#' @title htmltodata
#'
#' @param path is the path of the file that you want to import into R and then export.
#'
#' @importFrom XML htmlTreeParse getNodeSet xmlGetAttr xmlValue
#' @return The return data are a list include all text results of submitters' answers.
#' @author Xiaorui(Jeremy) Zhu
#' @export
htmltodata <- function(path){
  ## try to change the path into characters
  cha.path <- as.character(path)
  ### Use XML package to transfer html file into text
  root <- htmlTreeParse(cha.path)
  ### use getNodeSet to find right position of answers.
  getTitle <- getNodeSet(root, '//head//title')

  ### here for the file title
  title <- sapply(getTitle, xmlValue)

  values <- getNodeSet(root, '//body//div')

  text_values <- sapply(values, xmlValue) ### Just get those text answers, image lost
  ### image saving!!!!!!
  image_path <- getNodeSet(root, '//body//div[@class="field-value"]//img[@src]')
  image_saved_path <- sapply(image_path, function(els) xmlGetAttr(els, "src"))
  if (length(image_saved_path)>=1) {
    image <- paste(image_saved_path, sep = '') }
  else {
    image <- "No Image data"}
  list(title = title, text_values = text_values, Image_path = image)
}


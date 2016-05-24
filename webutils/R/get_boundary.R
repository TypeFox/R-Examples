get_boundary <- function(content_type){
  # Check for multipart
  if(!grepl("multipart/form-data; boundary=", content_type, fixed=TRUE))
    stop("Content type is not multipart/form-data: ", content_type)

  # Extract bounary
  sub("multipart/form-data; boundary=", "", content_type, fixed=TRUE)
}

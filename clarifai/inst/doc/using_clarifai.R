## ---- eval=F-------------------------------------------------------------
#  # install.packages("devtools")
#  devtools::install_github("soodoku/clarifai")

## ---- eval=F-------------------------------------------------------------
#  library(clarifai)

## ---- eval=F-------------------------------------------------------------
#  secret_id(c("client_id", "secret"))

## ---- eval=F-------------------------------------------------------------
#  get_token()

## ---- eval=F-------------------------------------------------------------
#  get_info()

## ---- eval=F-------------------------------------------------------------
#  res <- tag_image_urls("http://www.clarifai.com/img/metro-north.jpg")
#  head(res)

## ---- eval=F-------------------------------------------------------------
#  path <- system.file("inst/extdata/", package = "clarifai")
#  filep <- paste0(path, "/metro-north.jpg")
#  tag_images(filep)

## ---- eval=F-------------------------------------------------------------
#  feedback(file_path="path_to_image", feedback_type="add_tags", feedback_value="suggested_tag")


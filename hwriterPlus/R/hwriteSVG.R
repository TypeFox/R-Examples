hwriteSVG <- function(image.url, page = NULL,
                      width = 500, height = 500, id = "",
                      attributes = "border = '0'", ...) {

    ## Base string
    strng <- "<!--[if !IE]>-->
              <object data = @imageData type = 'image/svg+xml'
	              width = @imageWidth height = @imageHeight
                      id = @imageID @imageAttributes> <!--<![endif]-->
              <!--[if lt IE 9]>
              <object src = @imageData classid = 'image/svg+xml'
	              width = @imageWidth height = @imageHeight
                      id = @imageID @imageAttributes> <![endif]-->
              <!--[if gte IE 9]>
              <object data = @imageData type = 'image/svg+xml'
	              width = @imageWidth height = @imageHeight
                      id = @imageID @imageAttributes> <![endif]-->
              </object>"

    ## Substitute values in the string
    imageData <- paste("'", image.url, "'", sep = "")
    strng <- gsub("@imageData", imageData, strng)
    imageWidth <- paste("'", width, "'", sep = "")
    strng <- gsub("@imageWidth", imageWidth, strng)
    imageHeight <- paste("'", height, "'", sep = "")
    strng <- gsub("@imageHeight", imageHeight, strng)
    strng <- gsub("@imageID", id, strng)
    strng <- gsub("@imageAttributes", attributes, strng)

  ## final
  hwrite(strng, page, ...)
}

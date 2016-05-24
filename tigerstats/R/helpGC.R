#' @title Quick Vignettes in the Viewer

#' @description A convenience function to show package vignettes.  Vignette will show in a 
#' viewer selected by the front-end (e.g, the Viewer pane in R Studio), or in a 
#' browser if the "viewer" option is not set.
#' 
#' @rdname helpGC
#' @usage helpGC(topic,package="tigerstats")
#' @param topic filename of the vignette, exclusive of the .html extension
#' @param package Name of the installed package containing the desried vignette
#' @return side effects
#' @export
#' @author Homer White (hwhite0@@georgetowncollege.edu)
#' @examples
#' \dontrun{
#' helpGC(lmGC)
#' }
helpGC <- function(topic,package="tigerstats") {

# The following is outdated, as the rstudio package is no longer
# available:
#   if (!("rstudio" %in% rownames(installed.packages()))) {
#     return(cat("You need to be using R Studio for this helpGC() to help you!\n"))
#   }
  
  topic <- as.character(substitute(topic))
  vigList <- vignette(package=package)$results
  vigItems <- vigList[,"Item"]
  libPaths <- vigList[,"LibPath"]
  
  if (!(topic %in% vigItems)) {
    stop(paste0("Sorry, there is no vignette for ",topic," in package ",package,"."))
  } else {
    frow <- which(vigItems==topic)
    vigPath <- paste0(libPaths[frow],"/",package,"/doc/",topic,".html")
  }
  
  inCon  <- file(vigPath,"r")
  inStuff <- readLines(inCon)
  close(inCon)
  tempDir <- tempfile()
  dir.create(tempDir)
  htmlFile <- file.path(tempDir, "index.html")
  outCon <- htmlFile
  writeLines(inStuff,con=outCon)
  viewer <- getOption("viewer")
  if (!is.null(viewer))
    viewer(htmlFile)
  else
    utils::browseURL(htmlFile)
}


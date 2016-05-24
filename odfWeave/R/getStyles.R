getStyles <- function() get("odfStyles", envir = .odfEnv)

setStyles <- function(style)
{
  styleNames <- names(style)
  required <- c(
                "paragraph", "input", "output", "table", "cell", 
                "header", "cellText", "headerText", "page", "figureFrame")
  hasTypes <- required %in% styleNames
  
  if(any(!hasTypes))
    {
      missingList <- vector(mode = "list", length = sum(!hasTypes))
      names(missingList) <- required[!hasTypes]
      missingList <- lapply(missingList, function(x) "")
      style <- c(style, missingList)
    }
  
#   if(any(!hasTypes)) 
#      stop(
#         cat("There must be styles for:",
#            paste(required, collapse = ", "), "\n"))
#   typeCounts <- table(styleNames)
#   if(any(typeCounts > 1)) stop("only one style name can be declared here")
  assign( "odfStyles",  style, pos = .odfEnv)
}


getStyleDefs <- function() get("styleDefs", envir = .odfEnv)

setStyleDefs <- function(def)
{
  styleTypes <- unique(unlist(lapply(def, function(x) x$type)))
  required <- c("Paragraph", "Table Cell", "Table", "Bullet List", "Page", "Figure Frame")
  hasTypes <- required %in% styleTypes
  if(any(!hasTypes)) 
    stop(
         cat("There must be at least one style definintion for:",
             paste(required, collapse = ", "), "\n"))
  assign( "styleDefs",  def, pos = .odfEnv)
}

getImageDefs <- function() get("imageDefs", envir = .odfEnv)




setImageDefs <- 
  function (def, verbose = TRUE) 
{
  if (!(all(names(def) %in% c("type", "device", "plotHeight", 
                              "plotWidth", "dispHeight", "dispWidth", "args")))) 
    stop("invalid arguments were included. see ?setImageDefs")
  if (!is.null(def$args)) {
    if (!is.list(def$args)) 
      stop("args must be a list")
    if (any(names(def$args) %in% c("width", "height"))) 
      stop("these options should be specified without using the args list")
  }
  if (def$device %in% c("pdf", "svg")) 
    stop("pdf and svg formats not supported by OpenOffice")

  if (def$device == "postscript" & verbose) {
    psArgs <- names(def$args) %in% c("horizontal", "onefile", 
                                     "paper")
    psNote <- paste("you will probabiliy need to set", "\nhorizontal = FALSE, onefile = FALSE,", 
                    "and paper = \"special\" to", "generate ps graphics for OpenOffice\n")
    if (length(psArgs) == 0 | any(!psArgs)) 
      cat(psNote)
  }
  
  if (def$device %in% c("png", "bmp", "jpeg", "CairoPNG", "CairoJPEG", "CairoTIFF") & verbose) {
    if (def$plotHeight <= 30 | def$plotWidth <= 30) 
      cat(paste("an image size of", def$plotHeight, 
                "pixels by", def$plotWidth, "pixels has been requested.\n"))
  }
  else {
    if (def$plotHeight >= 30 | def$plotWidth >= 30) 
      cat(paste("an image size of", def$plotHeight, 
                "inches by", def$plotWidth, "inches has been requested.\n"))
  }
  flush.console()
  assign("imageDefs", def, pos = .odfEnv)
}

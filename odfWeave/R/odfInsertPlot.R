noCaptionXML <- function(name, width, height, units, file, anchor)
{
   # this is the code for odfWeave prior to version 1.5
   plotString <- c(
         "    <draw:frame ",
         paste("      draw:name=\"", name, "\"", sep = ""),
         "      text:anchor-type=\"paragraph\"",
         paste("      svg:width=\"", width, units, "\"", sep = ""),
         paste("      svg:height=\"", height, units, "\"", sep = ""),
         "      draw:z-index=\"0\">",
         "      <draw:image",
         paste("         xlink:href=\"Pictures/", basename(file), "\"", sep = ""),
         "         xlink:type=\"simple\"",
         "         xlink:show=\"embed\"",
         "         xlink:actuate=\"onLoad\"/>",
         "    </draw:frame>")
   
   out <- c(anchor[1], plotString, anchor[2], "\n")
   paste(out, collapse = "\n")
}

withCaptionXML <- function(name, width, height, units, file, anchor, caption)
{
   styles <- getStyles()
   frameDefs <-  getStyleDefs()[styles$figureFrame][[1]]

   # we need to add a little space between the frame and the image
   # use 5% of the width and height
   hFudge <- height * 0.05
   wFudge <- width * 0.05

   # try to check for padding and add to total size
   if(length(grep(units, frameDefs$padding) > 1))
   {
      paddingData <- strsplit(frameDefs$padding, units)[[1]]
      paddingAmount <- as.numeric(paddingData[1])     
      wFudge <- (width + paddingAmount) * 0.05
      hFudge <- (height + paddingAmount) * 0.05
            
   } else paddingAmount <- 0
   
   plotList <- list()
   plotList <- c(plotList, anchor[1])
   plotList <- c(plotList, '<draw:frame')
   plotList <- c(plotList, sprintf(' draw:style-name="%s"', styles$figureFrame))
   plotList <- c(plotList, sprintf(' draw:name="frame%d"', floor(runif(1) * 1000)))
   plotList <- c(plotList, sprintf(' svg:height="%s%s"', height + paddingAmount + hFudge, units))
   # XXX this seems to help, but needs some more thought
   plotList <- c(plotList, sprintf(' svg:width="%s%s"', width + paddingAmount+ wFudge, units))
   plotList <- c(plotList, ' svg:x="0in"')
   plotList <- c(plotList, ' svg:y="0in"')
   plotList <- c(plotList, sprintf(' text:anchor-type="%s"', frameDefs$frameAnchor))      

#   plotList <- c(plotList, paste(' text:anchor-type="', anchorType, '"', sep = ""))
   plotList <- c(plotList, '>')
   if (!is.null(caption))
   {
      plotList <- c(plotList, '<draw:text-box')
      plotList <- c(plotList, sprintf(' fo:min-height="%s%s"', height, units))
      plotList <- c(plotList, '>')
      plotList <- c(plotList, '<text:p text:style-name="Standard">')
      plotList <- c(plotList, '<draw:frame')
      plotList <- c(plotList, sprintf(' draw:name="%s"', name))
      plotList <- c(plotList, sprintf(' svg:width="%s%s"', width, units))
      plotList <- c(plotList, sprintf(' svg:height="%s%s"', height, units))
      plotList <- c(plotList, sprintf(' text:anchor-type="%s"', frameDefs$imageAnchor))      
      
      plotList <- c(plotList, '>')
   }
   plotList <- c(plotList, '<draw:image')
   plotList <- c(plotList, sprintf(' xlink:href="Pictures/%s"', basename(file)))
   plotList <- c(plotList, ' xlink:type="simple"')
   plotList <- c(plotList, ' xlink:show="embed"')
   plotList <- c(plotList, ' xlink:actuate="onLoad"/>')
   if (!is.null(caption))
   {
      plotList <- c(plotList, '</draw:frame>')
      plotList <- c(plotList, paste(caption$label, ' '))
      plotList <- c(plotList, '<text:sequence ')
      plotList <- c(plotList, sprintf(' style:num-format="%s"', caption$numformat))
      plotList <- c(plotList, sprintf(' text:formula="ooow:%s"', caption$formula))
      plotList <- c(plotList, ' text:name="Illustration"')
      plotList <- c(plotList, sprintf(' text:ref-name="%s"', caption$refname))
      if (!is.null(caption$numlettersync) && caption$numlettersync)
         plotList <- c(plotList, ' text:num-letter-sync="true"')
      plotList <- c(plotList, '>')
      plotList <- c(plotList, '999')
      plotList <- c(plotList, '</text:sequence>')
      plotList <- c(plotList, sprintf(': %s', caption$caption))
      plotList <- c(plotList, '</text:p>')
      plotList <- c(plotList, '</draw:text-box>')
   }
   plotList <- c(plotList, '</draw:frame>')
   plotList <- c(plotList, anchor[2])
   paste(plotList, collapse='\n')   
}

"odfInsertPlot" <-
function(file, height, width,
   units = "in",
   anchor = c('<text:p>', '</text:p>'),
   name = paste("graphics", floor(runif(1) * 1000), sep = ""),
   externalFile = FALSE,
   dest = paste(getwd(), "/Pictures", sep = ""),
   caption = NULL)
{
  
   if(getExt(file) %in% c("pdf")) stop("graphics format not supported")
   adjust <- if (!is.null(caption)) caption$adjust else 0.0

   if(is.null(caption))
   {
      out <- noCaptionXML(name, width, height, units, file, anchor)
   } else {
      out <- withCaptionXML(name, width, height, units, file, anchor, caption)
   }

   if(externalFile)
   {
      newPath <- paste(dest, "/", basename(file), sep = "")
      if(!file.exists(dest)) stop(paste(dest, "does not exist"))
      if(!file.copy(file, newPath,  overwrite = TRUE)) stop("Error copying file")
   }

   # append the basename of the image file to .odfEnv$picVector
   .odfEnv$picVector <- c(.odfEnv$picVector, basename(file))

   out
}


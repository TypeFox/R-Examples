figGen <- function(plotName = NULL)
{

   deviceArgs <- getImageDefs()
   # remove the non-device arguments
   deviceArgs <-  deviceArgs[!(names(deviceArgs) %in% "anchorType")]

   # make up a random name with prefix "rPlot"
   if(is.null(plotName)) plotName <- paste("rPlot", floor(runif(1) * 10000), 
      ".", deviceArgs$type, sep = "")

   if(!is.null(deviceArgs$args))
   {
      deviceArgs <- c(deviceArgs, deviceArgs$args)
      deviceArgs$args <- NULL
   }
   
   thisDevice <- deviceArgs$device
   deviceArgs$device <- deviceArgs$type <- NULL
   deviceArgs$dispHeight <- deviceArgs$dispWidth <- NULL
   
   deviceArgs$height <- deviceArgs$plotHeight
   deviceArgs$width <- deviceArgs$plotWidth
   deviceArgs$plotHeight <- deviceArgs$plotWidth <- NULL   
   
   if(thisDevice %in% c("postscript", "bitmap")) 
   {
      deviceArgs$file <- plotName   
   } else {
      deviceArgs$filename <- plotName
   }
   
   
   do.call(thisDevice, deviceArgs) 
}


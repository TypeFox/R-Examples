adjustImageSize <- function(x, y, scale = 1, res = 96)
{
   # x, y are width, height in inches
   # res is image resolution
   # scale is a scaling factor to easily adjust
   #   pch and text sizes

   temp <- getImageDefs()   
      
   if(temp$device %in% c("bmp", "jpeg", "png", "CairoPNG", "CairoJPEG", "CairoTIFF"))
   {
      temp$plotHeight <- res * y * scale
      temp$plotWidth <- res * x * scale
   } else {
      temp$plotHeight <- y * scale
      temp$plotWidth <- x * scale   
   }

   temp$dispHeight <- y
   temp$dispWidth <- x
   setImageDefs(temp)
   rm(temp)
   invisible()
}



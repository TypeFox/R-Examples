"pkgVersions" <-
function(type = "string", ncol = 4)
{
   if(!(type %in% c("string", "matrix", "data frame")))
      stop("type must be either: string, matrix or data frame")
   pkgs <- grep("package", search(), value = TRUE)
   pkgs <- lapply(strsplit(pkgs, ":"), function(data) data[2])
   pkgs <- lapply(pkgs, packageDescription, fields = c("Package", "Version"))
   pkgs <- lapply(pkgs, unlist)
   pkgDF <- data.frame(do.call("rbind", pkgs))
   pkgString <- paste(pkgDF[,1], " (", pkgDF[,2], ")", sep = "")
   pkgString <- pkgString[order(tolower(pkgString))]
   
   out <- switch(type,
      "data frame" = pkgDF,
      string = listString(pkgString, verbose = FALSE),
      matrix =
      {
         sizeDelta <- length(pkgString)%%ncol
         if(sizeDelta > 0)   out <- c(pkgString, rep("", ncol - length(pkgString) %% ncol))
            else out <- pkgString
         out <- matrix(out, ncol = ncol)
         out
      })

   out
}

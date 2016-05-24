##########################
#    Open user manual    #
##########################

"oppdf"<-function (file, bg = TRUE) 
{   
  file=file.path(paste(.libPaths(), "/OligoSpecificitySystem/doc/OSSUM.pdf",sep=""))
  OST <- .Platform$OS.type
  if (OST == "windows") 
    shell.exec(file)
    else if (OST == "unix")
    {
      bioCOpt <- getOption("BioC")
      pdf <- getOption("pdfviewer")
      if (is.null(pdf)) 
      {
         warning(paste("pdfViewer is set to:", pdf, "which does not seem to exist.  Please", 
                "run the command setOptionPdfViewer()"))
         return(FALSE)
      }
      cmd <- paste(pdf, file)
      if (bg) 
        cmd <- paste(cmd, "&")
        system(cmd)
    }
  return(TRUE)
}
"optxt"<-function(n){
if (n==1) file=file.path(paste(.libPaths(), "/OligoSpecificitySystem/doc/db1.txt",sep=""))
if (n==2) file=file.path(paste(.libPaths(), "/OligoSpecificitySystem/doc/db2.txt",sep=""))
if (n==3) file=file.path(paste(.libPaths(), "/OligoSpecificitySystem/doc/db3.txt",sep=""))
  fil<<-file
  OST <- .Platform$OS.type
  if (OST == "windows")
    shell.exec(file)
    else if (OST == "unix")
    {
      bioCOpt <- getOption("BioC")
      pdf <- getOption("pdfpager")
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
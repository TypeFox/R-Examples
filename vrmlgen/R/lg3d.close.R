
# terminates a Livegraphics3D-environment and write the
# corresponding 3D-scene to a VRML file
`lg3d.close` <-
function () 
{
    
    # read global Livegraphics3D parameters
    vrmlgenEnv <- get(".vrmlgenEnv",envir=.GlobalEnv)
    filename <- vrmlgenEnv$filename
    htmlout <- vrmlgenEnv$html
    VRMLDir <- vrmlgenEnv$VRMLDir
    hheight <- vrmlgenEnv$hheight
    hwidth <- vrmlgenEnv$hwidth
    scale <- vrmlgenEnv$scale
    col.bg <- vrmlgenEnv$col
    ambientlight <- vrmlgenEnv$ambientlight
        
    # catch calling of vrml.close without previous
    # vrml.open-call as an error and show warning message    
    curwarn.setting <- as.numeric(options("warn"))    
    options(warn = 2)
    
    if (data.class(result <- try(rm(".vrmlgenEnv", envir = .GlobalEnv), 
        TRUE)) == "try-error") {
        options(warn = curwarn.setting)
        warning("You are trying to apply vrml.close() without having used vrml.open() before")
    }
    
    options(warn = curwarn.setting)
    
    
    # create output files in the current directory
    curdir <- getwd()
    setwd(VRMLDir)
    
    write(paste("\n}, Boxed -> False, Axes -> False, AmbientLight->GrayLevel[", 
        ambientlight, "], Lighting -> True, BoxRatios -> Automatic, PlotRange -> All ]\n", 
        sep = ""), file = filename, append = TRUE)
        
    file.copy(filename, paste(curdir, filename, sep = "/"), overwrite = TRUE)
    file.remove(filename)
    
    datadir <- system.file("extdata", package = "vrmlgen")
    if (data.class(result<-try(file.copy(file.path(datadir, "live.jar"), file.path(curdir, "live.jar")), TRUE))=="try-error")
    {
      warning("\nCannot copy file live.jar from vrlmgen-folder to current directory. You might need to copy the file manually.")
    }
    
    setwd(curdir)
    if (!is.null(htmlout)) {
    
    		# create HTML output
        cat("<HTML>", file = htmlout, append = FALSE)
        cat("<HEAD><TITLE>VRMLGen-visualization</TITLE></HEAD><BODY>", 
            file = htmlout, append = TRUE)
        cat(paste("<APPLET ARCHIVE=\"live.jar\" CODE=\"Live.class\" WIDTH=", 
            hwidth, " HEIGHT=", hheight, " ALIGN=LEFT>", sep = ""), 
            file = htmlout, append = TRUE)
        coln <- col2rgb(col.bg)
        cat(paste("<PARAM NAME=\"BGCOLOR\" VALUE=\"", rgb(red = coln[1, 
            ]/255, green = coln[2, ]/255, blue = coln[3, ]/255), 
            "\">", sep = ""), file = htmlout, append = TRUE)
        cat(paste("<PARAM NAME=\"MAGNIFICATION\" VALUE=", scale, 
            ">", sep = ""), file = htmlout, append = TRUE)
        cat(paste("<PARAM NAME=\"INPUT_FILE\" VALUE=\"", filename, 
            "\">", sep = ""), file = htmlout, append = TRUE)
        cat("</APPLET></BODY>", file = htmlout, append = TRUE)
        cat("</HTML>", file = htmlout, append = TRUE)
    }
    
    # show success message
    message(paste("\nOutput file \"", filename, "\" was generated in folder ", 
        getwd(), ".\n", sep = ""))
}


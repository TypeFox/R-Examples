
# terminates a VRML-environment and write the
# corresponding 3D-scene to a VRML file
`vrml.close` <-
function () 
{
            
    # read global VRML parameters
    vrmlgenEnv <- get(".vrmlgenEnv",envir=.GlobalEnv)
    filename <- vrmlgenEnv$filename
    htmlout <- vrmlgenEnv$html
    VRMLDir <- vrmlgenEnv$VRMLDir
    hheight <- vrmlgenEnv$hheight
    hwidth <- vrmlgenEnv$hwidth
    
    
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
    
    file.copy(filename, paste(curdir, filename, sep = "/"), overwrite = TRUE)
    file.remove(filename)
    
    setwd(curdir)
    if (!is.null(htmlout)) {
    
    		# create HTML output
        cat("<HTML>", file = htmlout, append = FALSE)
        cat("<HEAD><TITLE>VRMLGen-visualization</TITLE></HEAD><BODY><br>", 
            file = htmlout, append = FALSE)
        cat(paste("<object type=\"x-world/x-vrml\" data=\"", 
            filename, "\" ", sep = ""), file = htmlout, append = TRUE)
        cat(paste("width=\"", hwidth, "\" height=\"", hheight, 
            "\"><br>", sep = ""), file = htmlout, append = TRUE)
        cat(paste("<param name=\"src\" value=\"", filename, "\"><br>", 
            sep = ""), file = htmlout, append = TRUE)
        cat("Your browser cannot display VRML files.<br>Please INSTALL a VRML-plugin or open the file in an external VRML-viewer.</object><br>", 
            file = htmlout, append = TRUE)
        cat("</BODY></HTML>", file = htmlout, append = TRUE)
    }
    
    # show success message
    message(paste("\nOutput file \"", filename, "\" was generated in folder ", 
        getwd(), ".\n", sep = ""))
}


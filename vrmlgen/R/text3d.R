
# plot text in a 3D-scene in the VRML- or Livegraphics3D-format
`text3d` <-
function (x, y = NULL, z = NULL, text, col = "black", fontweight = "normal", 
    fontfamily = "sans", hyperlink = NULL, rot = c(0, 1, 0, 0), 
    scale = 1) 
{
    
    # verify if text3d is called within a vrmlgen-environment
    # initialized with lg3d.open() or vrml.open()           
    if (exists(".vrmlgenEnv")) {
                        
        curdir <- getwd()
        vrmlgenEnv <- get(".vrmlgenEnv",envir=.GlobalEnv)
        setwd(vrmlgenEnv$VRMLDir)
        
        # extract the coordinate vectors from the input
        
        xyz_parse <- xyz.coords(x, y, z)
        x <- xyz_parse$x
        y <- xyz_parse$y
        z <- xyz_parse$z
        
        # check: does length of text vector match
        # to number of coordinates        
        if (length(text) != length(x)) 
            stop("\nNumber of coordinates does not match to length of string vector in text3d-function!\n")
        
        # create VRML output with lower-level drawing function
        if (vrmlgenEnv$type == "vrml") {
            
            # apply font settings
            if (fontweight == "normal") 
                fontweight <- "NONE"                
            else fontweight <- toupper(fontweight)
            
            fontfamily < toupper(fontfamily)
            
            # draw text
            .vrmltext(vrmlgenEnv$scale * x, vrmlgenEnv$scale * 
                y, vrmlgenEnv$scale * z, text, vrmlgenEnv$filename, 
                col, scale, rot, fontfamily, fontweight, hyperlink)
        }
        
        # create Livegraphics3D output with lower-level drawing function
        else {
        
        		# apply font settings
            if (fontweight == "bold") 
                fontweight <- "Bold"
            else fontweight <- "Plain"
            
            if (fontfamily == "sans") 
                fontfamily <- "Arial"
            if (fontfamily == "serif") 
                fontfamily <- "Times"
            
            # draw text
            .lg3dtext(x, y, z, text, vrmlgenEnv$filename, col, 
                scale, fontfamily, fontweight, hyperlink)
        }
        setwd(curdir)
    }
    
    # show warning, if user did not call vrml.open() or lg3d.open()
    else {
        warning("\nYou first need to call the wrapper-function vrml.open() or lg3d.open() in order to use the text3d-function.")
    }
}


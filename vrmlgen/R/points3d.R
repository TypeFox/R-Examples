
# plot data points in a 3D-scene in the VRML- or Livegraphics3D-format
`points3d` <-
function (x, y = NULL, z = NULL, col = "black", pointstyle = "s", 
    transparency = 0, hyperlinks = NULL, scale = 1) 
{
    
    # verify if points3d is called within a vrmlgen-environment
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
        
        # compute RGB color
        rcol <- (col2rgb(col)/255)[1]
        gcol <- (col2rgb(col)/255)[2]
        bcol <- (col2rgb(col)/255)[3]
        
        # call low-level plotting functions
        # for VRML or Livegraphics3D output
        
        if (vrmlgenEnv$type == "vrml") {
            
            .vrmlpoints(x, y, z, vrmlgenEnv$filename, rcol, 
                gcol, bcol, pointstyle, hyperlinks, vrmlgenEnv$scale, 
                scale, transparency)
        }
        else {
        
            .lg3dpoints(x, y, z, vrmlgenEnv$filename, rcol, 
                gcol, bcol, pointstyle, hyperlinks, 1, scale)
        }
        
        setwd(curdir)
    }
    
    # show warning, if user did not call vrml.open() or lg3d.open()
    else {
        message("\nYou first need to call the wrapper-function vrml.open() or lg3d.open() in order to use the points3d-function.")
    }
}


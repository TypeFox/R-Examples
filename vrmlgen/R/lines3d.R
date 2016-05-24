
# draw lines in a 3D-scene in the VRML- or Livegraphics3D-format
`lines3d` <-
function (x, y = NULL, z = NULL, col = "black", lwd = 1) 
{
    
    # verify if lines3d is called within a vrmlgen-environment
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
            
            # for VRML the start- and end-points
            # of lines must be converted into
            # start-points + rotation vectors
            
            for (i in seq(1, (length(x)-1), by = 2)) {
                              
                  # compute line length  
                  len <- sqrt((x[i + 1] - x[i])^2 + (y[i + 1] - 
                    y[i])^2 + (z[i + 1] - z[i])^2)
                    
                  
                  # get center of line
                  center <- 0.5 * (c(x[i], y[i], z[i]) + c(x[i + 
                    1], y[i + 1], z[i + 1]))
                  
                  
                  rot <- NULL  
                  if (len > 0) {
                  
                  	# compute rotation of vector to point
                  	# to the target position
                    target_pos <- c(x[i + 1], y[i + 1], z[i + 
                      1]) - c(x[i], y[i], z[i])
                      
                    target_pos <- target_pos/len
                    ey <- c(0, 1, 0)
                    
                    # compute rotation axis
                    rot_axis <- c(ey[2] * target_pos[3] - ey[3] * 
                      target_pos[2], ey[3] * target_pos[1] - 
                      ey[1] * target_pos[3], ey[1] * target_pos[2] - 
                      ey[2] * target_pos[1])
                      
                    # compute rotation angle
                    sine <- sqrt(rot_axis[1]^2 + rot_axis[2]^2 + 
                      rot_axis[3]^2)
                      
                    cosine <- crossprod(ey, target_pos)
                    angle <- atan2(sine, cosine)
                    
                    if ((abs(angle) < 1e-04) || (abs(angle - 
                      2 * pi) < 1e-04)) {
                      
                      rot <- NULL
                    }
                    else {
                      if (abs(sine) < 1e-04) 
                        rot_axis <- c(1, 0, 0)
                      rot <- c(rot_axis, angle)
                    }
                  }
                  
                  x[i] <- center[1]
                  y[i] <- center[2]
                  z[i] <- center[3]
                  
                  # draw line in VRML-format
                                    
                  if (is.null(rot)) 
                    write(paste("Transform {\n\ttranslation ", 
                      x[i] * vrmlgenEnv$scale, y[i] * vrmlgenEnv$scale, 
                      z[i] * vrmlgenEnv$scale, "\n\tchildren Shape {\n\t\tappearance Appearance { material Material {\n\t\tdiffuseColor ", 
                      rcol, gcol, bcol, " } }\n\tgeometry Cylinder { height ", 
                      len * vrmlgenEnv$scale, "radius", 0.04 * 
                        lwd, " }\n\t}\n}", sep = " "), file = vrmlgenEnv$filename, 
                      append = TRUE)
                  else write(paste("Transform {\n\ttranslation ", 
                    x[i] * vrmlgenEnv$scale, y[i] * vrmlgenEnv$scale, 
                    z[i] * vrmlgenEnv$scale, "\n\trotation", 
                    rot[1], rot[2], rot[3], rot[4], "\n\tchildren Shape {\n\t\tappearance Appearance { material Material {\n\t\tdiffuseColor ", 
                    rcol, gcol, bcol, " } }\n\tgeometry Cylinder { height ", 
                    len * vrmlgenEnv$scale, "radius", 0.04 * 
                      lwd, " }\n\t}\n}", sep = " "), file = vrmlgenEnv$filename, 
                    append = TRUE)
            }
        }
        else {
            for (i in seq(1, (length(x)-1), by = 2)) {
                                  
                  # draw line in Livegraphics3D-format
                  
                  write(paste("Thickness[", 0.005 * lwd/vrmlgenEnv$scale, 
                    "], RGBColor[", rcol, ",", gcol, ",", bcol, 
                    "], Line[{{", x[i] * vrmlgenEnv$scale, ",", 
                    y[i] * vrmlgenEnv$scale, ",", z[i] * vrmlgenEnv$scale, 
                    "},{", x[i + 1] * vrmlgenEnv$scale, ",", 
                    y[i + 1] * vrmlgenEnv$scale, ",", z[i + 
                      1] * vrmlgenEnv$scale, "}}],", sep = " "), 
                    file = vrmlgenEnv$filename, append = TRUE)

            }
        }
        setwd(curdir)
    }
    
    # show warning, if user did not call vrml.open() or lg3d.open()
    else {
        message("\nYou first need to call the wrapper-function vrml.open() or lg3d.open() in order to use the lines3d-function.")
    }
}


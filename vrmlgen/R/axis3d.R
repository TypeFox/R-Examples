
# add coordinate axes to a 3D-scene in the VRML- or Livegraphics3D-format
`axis3d` <-
function (lab.axis = c("X-axis", "Y-axis", "Z-axis"), filename = NULL, 
    type = "vrml", col.lab = "black", col.axis = "black", cex.lab = 1, 
    local_scale = 1, global_scale = 1) 
{
    
    # verify if axis3d is called within a vrmlgen-environment
    # initialized with lg3d.open() or vrml.open()           
    curdir <- getwd()
    if (exists(".vrmlgenEnv")) {
    
    		# read global 3D scene settings
    		vrmlgenEnv <- get(".vrmlgenEnv",envir=.GlobalEnv)
        type <- vrmlgenEnv$type
        filename <- vrmlgenEnv$filename
        global_scale <- vrmlgenEnv$scale
        setwd(vrmlgenEnv$VRMLDir)
    }
    
    # set the axis and label colors
    lab_rcol <- (col2rgb(col.lab)/255)[1]
    lab_gcol <- (col2rgb(col.lab)/255)[2]
    lab_bcol <- (col2rgb(col.lab)/255)[3]
    
    ax_rcol <- (col2rgb(col.axis)/255)[1]
    ax_gcol <- (col2rgb(col.axis)/255)[2]
    ax_bcol <- (col2rgb(col.axis)/255)[3]
    
    # write the axes in the specified format
    if (type == "vrml") {
    
    		# write axis labels
        write(paste("Transform {\n\ttranslation ", global_scale * 
            local_scale + 0.8, " 0 0\n\tscale ", cex.lab * 0.28, 
            cex.lab * 0.28, cex.lab * 0.28, "\n\trotation 0.00000 0.70711 0.70711 3.14159\n\tchildren Shape {\n\t\tappearance Appearance { material Material {diffuseColor ", 
            lab_rcol, lab_gcol, lab_bcol, "  } }\n\tgeometry Text { string \"", 
            lab.axis[1], "\"\nfontStyle FontStyle {\nfamily [\"SANS\"]\njustify \"MIDDLE\"\n }\n }\n\t}\n}",
            sep = " "), file = filename, 
            append = TRUE)
            
        write(paste("Transform {\n\ttranslation 0 ", global_scale * 
            local_scale + 0.8, " 0\n\tscale ", cex.lab * 0.28, 
            cex.lab * 0.28, cex.lab * 0.28, "\n\trotation 0.57735 0.57735 0.57735 2.09440\n\tchildren Shape {\n\t\tappearance Appearance { material Material {diffuseColor ", 
            lab_rcol, lab_gcol, lab_bcol, "  } }\n\tgeometry Text { string \"", 
            lab.axis[2], "\"\nfontStyle FontStyle {\nfamily [\"SANS\"]\njustify \"MIDDLE\"\n }\n }\n\t}\n}",
            sep = " "), file = filename, 
            append = TRUE)
            
        write(paste("Transform {\n\ttranslation 0 0 ", global_scale * 
            local_scale + 0.8, "\n\tscale ", cex.lab * 0.28, 
            cex.lab * 0.28, cex.lab * 0.28, "\n\trotation 0.28108 0.67860 0.67860 2.59356\n\tchildren Shape {\n\t\tappearance Appearance { material Material {diffuseColor ", 
            lab_rcol, lab_gcol, lab_bcol, "  } }\n\tgeometry Text { string \"", 
            lab.axis[3], "\"\nfontStyle FontStyle {\nfamily [\"SANS\"]\njustify \"MIDDLE\"\n }\n }\n\t}\n}",
            sep = " "), file = filename, 
            append = TRUE)
        
        # write axes   
        write(paste("Transform {\n\ttranslation ", 0.5 * global_scale * 
            local_scale, " 0 0\n\trotation 0 0 1 1.5708\n\tchildren Shape {\n\t\tappearance Appearance { material Material {\n\t\tdiffuseColor ", 
            ax_rcol, ax_gcol, ax_bcol, " } }\n\tgeometry Cylinder { height ", 
            global_scale * local_scale, " radius 0.04 }\n\t}\n}", 
            sep = " "), file = filename, append = TRUE)
            
        write(paste("Transform {\n\ttranslation 0 ", 0.5 * global_scale * 
            local_scale, " 0\n\trotation 0 0 1 0\n\tchildren Shape {\n\t\tappearance Appearance { material Material { diffuseColor ", 
            ax_rcol, ax_gcol, ax_bcol, " } }\n\tgeometry Cylinder { height ", 
            global_scale * local_scale, " radius 0.04 }\n\t}\n}", 
            sep = " "), file = filename, append = TRUE)
            
        write(paste("Transform {\n\ttranslation 0 0 ", 0.5 * 
            global_scale * local_scale, "\n\trotation 1 0 0 1.5708\n\tchildren Shape {\n\t\tappearance Appearance { material Material { diffuseColor ", 
            ax_rcol, ax_gcol, ax_bcol, " } }\n\tgeometry Cylinder { height ", 
            global_scale * local_scale, " radius 0.04 }\n\t}\n}", 
            sep = " "), file = filename, append = TRUE)
    }
    else {
        
        # write axis labels
        write(paste("RGBColor[", lab_rcol, ",", lab_gcol, ",", 
            lab_bcol, "], Text [ \"", lab.axis[1], "\", {", local_scale + 
                0.5, ",0,0 }],\n", sep = ""), file = filename, 
            append = TRUE)
            
        write(paste("RGBColor[", lab_rcol, ",", lab_gcol, ",", 
            lab_bcol, "], Text [ \"", lab.axis[2], "\", {0,", 
            local_scale + 0.5, ",0 }],\n", sep = ""), file = filename, 
            append = TRUE)
            
        write(paste("RGBColor[", lab_rcol, ",", lab_gcol, ",", 
            lab_bcol, "], Text [ \"", lab.axis[3], "\", {0, 0,", 
            local_scale + 0.5, " }],\n", sep = ""), file = filename, 
            append = TRUE)
         
        # write axes  
        write(paste("Thickness[0.005], RGBColor[", ax_rcol, ",", 
            ax_gcol, ",", ax_bcol, "], Line[{{0,0,0},{", local_scale, 
            ",0,0}}],", sep = " "), file = filename, append = TRUE)
            
        write(paste("Thickness[0.005], RGBColor[", ax_rcol, ",", 
            ax_gcol, ",", ax_bcol, "], Line[{{0,0,0},{0,", local_scale, 
            ",0}}],", sep = " "), file = filename, append = TRUE)
            
        write(paste("Thickness[0.005], RGBColor[", ax_rcol, ",", 
            ax_gcol, ",", ax_bcol, "], Line[{{0,0,0},{0,0,", 
            local_scale, "}}],", sep = " "), file = filename, 
            append = TRUE)
    }
    
    if (exists(".vrmlgenEnv")) 
        setwd(curdir)
}


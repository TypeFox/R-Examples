
# visualize 3D meshes and parametric functions in the VRML- or Livegraphics3D-format.
`mesh3d` <-
function (xfun = "sin(v)*cos(u)", yfun = "sin(v)*sin(u)", zfun = "cos(v)", 
    param1 = "u", param2 = "v", range1 = c(0, 2 * pi), range2 = c(0, 
        pi), size1 = 30, size2 = 30, type = "vrml", x = NULL, 
    y = NULL, z = NULL, edges = NULL, obj_infile = NULL, filename = "out.wrl", 
    write_obj = FALSE, cols = "red", scalefac = 4, autoscale = ifelse(is.null(obj_infile), 
        "independent", "equicenter"), lab.axis = c("X-axis", 
        "Y-axis", "Z-axis"), col.axis = "black", showaxis = TRUE, 
    col.lab = "black", col.bg = "white", cex.lab = 1, htmlout = NULL, 
    hwidth = 1200, hheight = 800, vrml_navigation = "EXAMINE", 
    vrml_transparency = 0, vrml_fov = 0.785, vrml_pos = rep(scalefac + 
        4, 3), vrml_dir = c(0.19, 0.45, 0.87, 2.45), lg3d_ambientlight = 0.5) 
{
    
    # if called within a higher-level VRML or 
    # Livegraphics3D-environment, get the plot type
    if (exists(".vrmlgenEnv")) {
    		vrmlgenEnv <- get(".vrmlgenEnv",envir=.GlobalEnv)
        type <- vrmlgenEnv$type
    }
    
    # call lower-level plotting functions
    
    if ((type == "vrml") || write_obj) {
        
        .vmesh(xfun = xfun, yfun = yfun, zfun = zfun, param1 = param1, 
            param2 = param2, range1 = range1, range2 = range2, 
            size1 = size1, size2 = size2, x = x, y = y, z = z, 
            edges = edges, infile = obj_infile, filename = filename, 
            write_obj = write_obj, cols = cols, scalefac = scalefac, 
            autoscale = autoscale, lab.axis = lab.axis, col.axis = col.axis, 
            showaxis = showaxis, col.lab = col.lab, col.bg = col.bg, 
            cex.lab = cex.lab, navigation = vrml_navigation, 
            transparency = vrml_transparency, fov = vrml_fov, 
            pos = vrml_pos, dir = vrml_dir, htmlout = htmlout, 
            hwidth = hwidth, hheight = hheight)
    }
    else {
        if (!is.null(obj_infile)) {
            stop("\nReading input files in obj-format is only supported in combination with the VRML-output format.")
        }
        
        .lmesh(xfun = xfun, yfun = yfun, zfun = zfun, param1 = param1, 
            param2 = param2, range1 = range2, range2 = range2, 
            size1 = size1, size2 = size2, filename = filename, 
            cols = cols, scalefac = scalefac, autoscale = autoscale, 
            lab.axis = lab.axis, col.axis = col.axis, showaxis = showaxis, 
            col.lab = col.lab, col.bg = col.bg, cex.lab = cex.lab, 
            ambientlight = lg3d_ambientlight, htmlout = htmlout, 
            hwidth = hwidth, hheight = hheight)
    }
}


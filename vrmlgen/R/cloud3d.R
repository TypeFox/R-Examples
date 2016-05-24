
# create 3D scatter plot visualizations in the VRML- or Livegraphics3D-format
`cloud3d` <-
function (x, y = NULL, z = NULL, labels = rownames(data), filename = "out.wrl", 
    type = "vrml", pointstyle = c("s", "b", "c"), metalabels = NULL, 
    hyperlinks = NULL, cols = rainbow(length(unique(labels))), 
    scalefac = 4, autoscale = "independent", lab.axis = c("X-axis", 
        "Y-axis", "Z-axis"), col.axis = "black", showaxis = TRUE, 
    col.lab = "black", col.bg = "white", cex.lab = 1, htmlout = NULL, 
    hwidth = 1200, hheight = 800, showlegend = TRUE, vrml_navigation = "EXAMINE", 
    vrml_showdensity = FALSE, vrml_fov = 0.785, vrml_pos = rep(scalefac + 
        4, 3), vrml_dir = c(0.19, 0.45, 0.87, 2.45), vrml_transparency = 0, 
    lg3d_ambientlight = 0.5) 
{
    
    # if called within a higher-level VRML or 
    # Livegraphics3D-environment, get the plot type
    if (exists(".vrmlgenEnv")) {
    		vrmlgenEnv <- get(".vrmlgenEnv",envir=.GlobalEnv)
        type <- vrmlgenEnv$type
    }
    
    
    # call low-level plotting functions
    # for VRML or Livegraphics3D output
    
    if (type == "vrml") {
        .vcloud(x, y, z, labels = labels, filename = filename, 
            pointstyle = pointstyle, metalabels = metalabels, 
            hyperlinks = hyperlinks, cols = cols, scalefac = scalefac, 
            autoscale = autoscale, lab.axis = lab.axis, col.axis = col.axis, 
            col.lab = col.lab, col.bg = col.bg, cex.lab = cex.lab, 
            navigation = vrml_navigation, transparency = vrml_transparency, 
            fov = vrml_fov, pos = vrml_pos, dir = vrml_dir, htmlout = htmlout, 
            hwidth = hwidth, hheight = hheight, showlegend = showlegend, 
            showdensity = vrml_showdensity)
    }
    else {
        .lcloud(x, y, z, labels = labels, filename = filename, 
            pointstyle = pointstyle, metalabels = metalabels, 
            hyperlinks = hyperlinks, cols = cols, scalefac = scalefac, 
            autoscale = autoscale, lab.axis = lab.axis, col.axis = col.axis, 
            col.lab = col.lab, col.bg = col.bg, cex.lab = cex.lab, 
            htmlout = htmlout, hwidth = hwidth, hheight = hheight, 
            showlegend = showlegend, ambientlight = lg3d_ambientlight)
    }
}


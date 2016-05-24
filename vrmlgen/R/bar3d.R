
# create 3D barplots and height map visualizations
# in the VRML- or Livegraphics3D-format.
`bar3d` <-
function (data, row.labels = rownames(data), col.labels = colnames(data), 
    metalabels = NULL, filename = "out.wrl", type = "vrml", space = 0.5, 
    cols = rainbow(length(as.matrix(data))), rcols = NULL, ccols = NULL, 
    origin = c(0, 0, 0), scalefac = 4, lab.axis = c("X-axis", 
        "Y-axis", "Z-axis"), lab.vertical = FALSE, col.axis = "black", 
    showaxis = TRUE, autoscale = TRUE, ignore_zeros = TRUE,
    col.lab = "black", col.bg = "white", cex.lab = 1, 
    cex.rowlab = 1, cex.collab = 1, htmlout = NULL, hwidth = 1200, 
    hheight = 800, showlegend = TRUE, vrml_navigation = "EXAMINE", 
    vrml_transparency = 0, vrml_fov = 0.785, vrml_pos = rep(scalefac + 
        4, 3), vrml_dir = c(0.19, 0.45, 0.87, 2.45), lg3d_ambientlight = 0.5) 
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
        .vbar(data, row.labels, col.labels, metalabels = metalabels, 
            filename = filename, space = space, cols = cols, 
            rcols = rcols, ccols = ccols, origin = origin, scalefac = scalefac, 
            lab.axis = lab.axis, lab.vertical = lab.vertical, 
            col.axis = col.axis, showaxis = showaxis, autoscale = autoscale, 
    				ignore_zeros = ignore_zeros, col.lab = col.lab, col.bg = col.bg,
    				cex.lab = cex.lab, cex.rowlab = cex.rowlab, 
            cex.collab = cex.collab, fov = vrml_fov, pos = vrml_pos, 
            dir = vrml_dir, htmlout = htmlout, hwidth = hwidth, 
            hheight = hheight, showlegend = showlegend)
    }
    else {
        .lbar(data, row.labels, col.labels, filename = filename, 
            space = space, cols = cols, rcols = rcols, ccols = ccols, 
            origin = origin, scalefac = scalefac, lab.axis = lab.axis, 
            col.axis = col.axis, showaxis = showaxis, autoscale = autoscale, 
    				ignore_zeros = ignore_zeros, col.lab = col.lab, 
            col.bg = col.bg, cex.lab = cex.lab, cex.rowlab = cex.rowlab, 
            cex.collab = cex.collab, ambientlight = lg3d_ambientlight, 
            htmlout = htmlout, hwidth = hwidth, hheight = hheight, 
            showlegend = showlegend)
    }
}


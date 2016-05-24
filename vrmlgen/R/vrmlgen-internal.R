# helper function to apply the automatic coloring of
# labeled datapoints in different plotting functions
`.colorselection` <-
function (type = "cloud", index, cols, numpoints, labels, numlabels) 
{
    
    # automatic coloring for bar plot functions
    if (type == "bar") {
        if (length(cols) == length(labels)) {
            rcol <- (col2rgb(cols[index])/255)[1]
            gcol <- (col2rgb(cols[index])/255)[2]
            bcol <- (col2rgb(cols[index])/255)[3]
        }
        else if (length(cols) == 1) {
            rcol <- (col2rgb(cols[1])/255)[1]
            gcol <- (col2rgb(cols[1])/255)[2]
            bcol <- (col2rgb(cols[1])/255)[3]
        }
        else {
            rcol <- (col2rgb("black")/255)[1]
            gcol <- (col2rgb("black")/255)[2]
            bcol <- (col2rgb("black")/255)[3]
        }
    }
    
    # automatic coloring for 3D cloud and mesh plots
    else {
        if (!length(cols)) {
            rcol <- (col2rgb(rainbow(numpoints)[index])/255)[1]
            gcol <- (col2rgb(rainbow(numpoints)[index])/255)[2]
            bcol <- (col2rgb(rainbow(numpoints)[index])/255)[3]
        }
        else if ((length(cols) == 1) || ((length(cols) < numpoints) && 
            !length(labels))) {
            rcol <- (col2rgb(cols[1])/255)[1]
            gcol <- (col2rgb(cols[1])/255)[2]
            bcol <- (col2rgb(cols[1])/255)[3]
        }
        else if (length(labels) && (length(cols) == length(unique(labels))) && 
            (type != "mesh")) {
            rcol <- (col2rgb(cols[numlabels[index]])/255)[1]
            gcol <- (col2rgb(cols[numlabels[index]])/255)[2]
            bcol <- (col2rgb(cols[numlabels[index]])/255)[3]
        }
        else {
            rcol <- (col2rgb(cols[index])/255)[1]
            gcol <- (col2rgb(cols[index])/255)[2]
            bcol <- (col2rgb(cols[index])/255)[3]
        }
    }
    return(list(rcol = rcol, gcol = gcol, bcol = bcol))
}

# creating bar plots in Livegraphics3D-format
`.lbar` <-
function (data, row.labels = rownames(data), col.labels = colnames(data), 
    filename = "out.m", space = 0.5, cols = rainbow(length(as.matrix(data))), 
    rcols = NULL, ccols = NULL, origin = c(0, 0, 0), scalefac = 4, autoscale = TRUE, 
    ignore_zeros = TRUE, lab.axis = c("X-axis", "Y-axis", "Z-axis"), col.axis = "black", 
    showaxis = TRUE, col.lab = "black", col.bg = "white", cex.lab = 1, 
    cex.rowlab = 1, cex.collab = 1, ambientlight = 0.5, htmlout = NULL, 
    hwidth = 1200, hheight = 800, showlegend = TRUE) 
{
    
    data <- as.matrix(data)
    
    curdir <- NULL
    
    # verify if .lbar is called within a vrmlgen-environment
    # initialized with lg3d.open()
    if (!exists(".vrmlgenEnv")) {
        
        # set the general 3D scene settings, if they are
        # not provided by a higher-level environment
        
        # write LiveGraphics3D header
        write("Graphics3D[\n{\n", file = filename, append = FALSE)
        
        # set background color
        bg_rcol <- (col2rgb(col.bg)/255)[1]
        bg_gcol <- (col2rgb(col.bg)/255)[2]
        bg_bcol <- (col2rgb(col.bg)/255)[3]
        
        write(paste("Background {\n\t skyColor [\n\t\t ", bg_rcol, 
            bg_gcol, bg_bcol, " \n\t]\n}", sep = " "), file = filename, 
            append = TRUE)
    }
    else {
    
    		# obtain the general settings from the
    		# higher level environment
    		
    		vrmlgenEnv <- get(".vrmlgenEnv",envir=.GlobalEnv)
        filename <- vrmlgenEnv$filename
        htmlout <- vrmlgenEnv$html
        hheight <- vrmlgenEnv$hheight
        hwidth <- vrmlgenEnv$hwidth
        curdir <- getwd()
        setwd(vrmlgenEnv$VRMLDir)
        
    }
    
    # scale the data to fit within the axes ranges
    if(autoscale)
    	data <- scalefac * (data/max(data))
    
    
    # pre-process row- and column-labels
    if (!is.null(col.labels) && showlegend) 
        col.labels <- sapply(strsplit(col.labels, " "), function(x) paste(x, 
            collapse = "\n"))
            
    if (!is.null(row.labels) && showlegend) 
        row.labels <- sapply(strsplit(row.labels, " "), function(x) paste(x, 
            collapse = "\n"))
            
    # calculate bar-length and -width (account for scaling and spacing)
    blength <- scalefac/(nrow(data) * (1 + space))
    bwidth <- scalefac/(ncol(data) * (1 + space))
    
    # draw the axes
    if (showaxis) {
    
    		# call the axis3d helper function to draw the axes
        axis3d(lab.axis, filename, type = "lg3d", col.lab, col.axis, 
            cex.lab, local_scale = scalefac, global_scale = 1)
            
    }
    
    # draw row-labels
    if (!is.null(row.labels) && showlegend) {
    
        for (j in 1:length(row.labels)) {
            
            # apply automatic color selection
            colsel <- .colorselection(type = "bar", j, rcols, 
                1, row.labels, "")
            
            rcol <- colsel$rcol
            gcol <- colsel$gcol
            bcol <- colsel$bcol
            
            cur_xwidth <- j/nrow(data) * scalefac - blength/2
            
            write(paste("RGBColor[", rcol, ",", gcol, ",", bcol, 
                "], Text [ \"", row.labels[j], "\", {", cur_xwidth, 
                ",", -0.4, ",", scalefac + 0.2, "}],\n", sep = ""), 
                file = filename, append = TRUE)
        }
    }
    
    # draw column-labels
    if (!is.null(col.labels) && showlegend) {
        
        for (j in 1:length(col.labels)) {
            
            # apply automatic color selection
            colsel <- .colorselection(type = "bar", j, ccols, 
                1, col.labels, "")
                
            rcol <- colsel$rcol
            gcol <- colsel$gcol
            bcol <- colsel$bcol
            
            cur_ywidth <- j/ncol(data) * scalefac - bwidth/2
            
            write(paste("RGBColor[", rcol, ",", gcol, ",", bcol, 
                "], Text [ \"", col.labels[j], "\", {", -0.4, 
                ",", cur_ywidth, ",", scalefac + 0.2, "}],\n", 
                sep = ""), file = filename, append = TRUE)
        }
    }
    
    
    rcol <- NULL
    gcol <- NULL
    bcol <- NULL
    
    # set RGB-colors
    if (length(cols) >= (ncol(data) * nrow(data))) {
        rcol <- sapply(cols, function(x) (col2rgb(x)/255)[1])
        gcol <- sapply(cols, function(x) (col2rgb(x)/255)[2])
        bcol <- sapply(cols, function(x) (col2rgb(x)/255)[3])
    }
    else {    		
        if (length(cols) >= 1) {
            rcol <- (col2rgb(cols[1])/255)[1]
            gcol <- (col2rgb(cols[1])/255)[2]
            bcol <- (col2rgb(cols[1])/255)[3]
        }
        else {
            rcol <- (col2rgb("lightblue")/255)[1]
            gcol <- (col2rgb("lightblue")/255)[2]
            bcol <- (col2rgb("lightblue")/255)[3]
        }
    }
    
    
    # main loop: iterate over data points
    
    bwidth <- bwidth/2
    for (k in 1:ncol(data)) {
    
        for (j in 1:nrow(data)) {
        
            x <- j/nrow(data) * scalefac
            y <- k/ncol(data) * scalefac
            z <- data[j, k]
            
            
        		if(ignore_zeros && (z == 0))
  							next
        
            
            # set current color
            if (!is.null(ccols)) {
                rcol <- (col2rgb(ccols[k])/255)[1]
                gcol <- (col2rgb(ccols[k])/255)[2]
                bcol <- (col2rgb(ccols[k])/255)[3]
            }
            else if (!is.null(rcols)) {
                rcol <- (col2rgb(rcols[j])/255)[1]
                gcol <- (col2rgb(rcols[j])/255)[2]
                bcol <- (col2rgb(rcols[j])/255)[3]
            }
            else {
                if(length(rcol) > 1)
                {
	                rcol <- rcol[(k - 1) * ncol(data) + 
	                  j]
	                gcol <- gcol[(k - 1) * ncol(data) + 
	                  j]
	                bcol <- bcol[(k - 1) * ncol(data) + 
	                  j]
	              }
            }
            
            # draw data point
            write(paste("SurfaceColor[RGBColor[", rcol, ",", 
                gcol, ",", bcol, "]],  Cuboid[{", origin[1] + 
                  x - bwidth, ",", origin[2] + y - bwidth, ",", 
                origin[3], "},{", origin[1] + x + bwidth, ",", 
                origin[2] + y + bwidth, ",", origin[3] + z, "}],", 
                sep = ""), file = filename, append = TRUE)
        }
    }
    
    # configure general plot settings (font size, text style, etc.)
    write(paste("\n}, Boxed -> False, Axes -> False, AmbientLight->GrayLevel[", 
        ambientlight, "], TextStyle -> {FontFamily -> \"TimesRoman\", FontSlant ->\"Italic\", FontSize -> ", 
        14 * cex.lab, "}, Lighting -> True, BoxRatios -> Automatic, PlotRange -> All ]\n", 
        sep = ""), file = filename, append = TRUE)
        
    
    # create the output files, if this is not already
    # done by a higher-level lg3d.close-function    
    if (!exists(".vrmlgenEnv")) {
    
    		# create HTML output
        if (!is.null(htmlout)) {
            cat("<HTML>", file = htmlout, append = FALSE)
            cat("<HEAD><TITLE>VRMLGen-visualization</TITLE></HEAD><BODY>", 
                file = htmlout, append = TRUE)
            cat(paste("<APPLET ARCHIVE=\"live.jar\" CODE=\"Live.class\" WIDTH=", 
                hwidth, " HEIGHT=", hheight, " ALIGN=LEFT>", 
                sep = ""), file = htmlout, append = TRUE)
            coln <- col2rgb(col.bg)
            cat(paste("<PARAM NAME=\"BGCOLOR\" VALUE=\"", rgb(red = coln[1, 
                ]/255, green = coln[2, ]/255, blue = coln[3, 
                ]/255), "\">", sep = ""), file = htmlout, append = TRUE)
            cat("<PARAM NAME=\"MAGNIFICATION\" VALUE=1.0>", file = htmlout, 
                append = TRUE)
            cat(paste("<PARAM NAME=\"INPUT_FILE\" VALUE=\"", 
                filename, "\">", sep = ""), file = htmlout, append = TRUE)
            cat("</APPLET>", file = htmlout, append = TRUE)
            cat("</HTML>", file = htmlout, append = TRUE)
        }
        
        # copy the Livegraphics3D live.jar to the output directory                
        datadir <- system.file("extdata", package = "vrmlgen")       	
       	curdir <- getwd()
       	
       	# use if file.exists
        if (data.class(result<-try(file.copy(file.path(datadir, "live.jar"), file.path(curdir, "live.jar")), TRUE))=="try-error")
        {
          warning("\nCannot copy file live.jar from vrlmgen-folder to current directory. You might need to copy the file manually.")
        }
        
        # show success message
        cat(paste("\nOutput file \"", filename, "\" was generated in folder ", 
            getwd(), ".\n", sep = ""))
                    
    }
    else {
        # return to current directory, if higher-level plotting
        # functions wrote re-directed the output to a temporary directory
        setwd(curdir)
    }
}

# creating scatter plots in Livegraphics3D-format
`.lcloud` <-
function (x, y = NULL, z = NULL, labels = rownames(data), metalabels = NULL, 
    hyperlinks = NULL, filename = "out.m", pointstyle = c("s", 
        "b", "t"), cols = rainbow(length(unique(labels))), scalefac = 4, 
    autoscale = "independent", lab.axis = c("X-axis", "Y-axis", 
        "Z-axis"), col.axis = "black", showaxis = TRUE, col.lab = "black", 
    col.metalab = "black", col.bg = "white", cex.lab = 1, ambientlight = 0.5, 
    htmlout = NULL, hwidth = 1200, hheight = 800, showlegend = TRUE) 
{
		
		# extract the coordinate vectors from the input    
    xyz_parse <- xyz.coords(x, y, z)
    
    data <- cbind(xyz_parse$x, xyz_parse$y, xyz_parse$z)
    if (ncol(data) != 3) {
        stop("Data matrix does not have 3 columns!")
    }
    
    # identify the unique set of row labels
    numlabels <- NULL
    if (length(labels)) {
    
        lab <- unique(unlist(labels))
        numlabels <- apply(as.matrix(labels), 1, function(x) match(x, 
            as.matrix(lab)))
            
    }
    
    
    
        
    
    # verify if .lcloud is called within a vrmlgen-environment
    # initialized with lg3d.open()
    curdir <- NULL
    if (!exists(".vrmlgenEnv")) {

        # set the general 3D scene settings, if they are
        # not provided by a higher-level environment

        write("Graphics3D[\n{\n", file = filename, append = FALSE)
        
        # set background color
        bg_rcol <- (col2rgb(col.bg)/255)[1]
        bg_gcol <- (col2rgb(col.bg)/255)[2]
        bg_bcol <- (col2rgb(col.bg)/255)[3]
        
        write(paste("Background {\n\t skyColor [\n\t\t ", bg_rcol, 
            bg_gcol, bg_bcol, " \n\t]\n}, \n", sep = " "), file = filename, 
            append = TRUE)
    }
    else {
        
        # obtain the general settings from the
    		# higher level environment
    		vrmlgenEnv <- get(".vrmlgenEnv",envir=.GlobalEnv)
        filename <- vrmlgenEnv$filename
        htmlout <- vrmlgenEnv$html
        hheight <- vrmlgenEnv$hheight
        hwidth <- vrmlgenEnv$hwidth
        curdir <- getwd()
        setwd(vrmlgenEnv$VRMLDir)
    }
    
    # apply the chosen scaling method
    scaledat <- function(data) {
        diff <- max(data) - min(data)
        if (diff == 0) 
            return(data)
        else return(scalefac * (data - min(data))/(diff))
    }
    
    if (autoscale == "independent") 
        data <- apply(data, 2, scaledat)
    else if (autoscale == "equidist") 
        data <- scalefac * (data/max(data))
    
    # create the plot legend
    if (length(labels) && showlegend) {
    
        cur_height <- scalefac + 1 + cex.lab/5
        
        for (j in 1:length(unique(labels))) {
        
        		# write text using the lower-level text
        		# drawing function
            .lg3dtext(0, 0, cur_height, unique(labels)[j], filename, 
                cols[unique(numlabels)[j]], 1, "Arial", "Plain", 
                NULL)
            cur_height <- cur_height + 0.4
        }
    }
    
    
    # draw the axes
    
    lab_rcol <- NULL
    lab_gcol <- NULL
    lab_bcol <- NULL
    ax_rcol <- NULL
    ax_gcol <- NULL
    ax_bcol <- NULL
    if (showaxis) {
    
    		# call the axis3d helper function to draw the axes
        axis3d(lab.axis, filename, type = "lg3d", col.lab, col.axis, 
            cex.lab, local_scale = scalefac, global_scale = 1)
            
    }
    
    # use empty metalabels if only url-labels have been defined
    if ((length(metalabels) == 0) && (length(hyperlinks) >= 1)) {
    
        metalabels <- rep("  ", length(hyperlinks))
    }
    
    # main loop: iterate over data points    
    for (j in 1:nrow(data)) {
    
        x <- data[j, 1]
        y <- data[j, 2]
        z <- data[j, 3]
        
        # apply automatic color selection
        colsel <- .colorselection(type = "cloud", j, cols, nrow(data), 
            labels, numlabels)
            
        rcol <- colsel$rcol
        gcol <- colsel$gcol
        bcol <- colsel$bcol
        
        # draw data points with given point-style
        
        if (length(pointstyle) == 1) {
        
            .lg3dpoints(x, y, z, filename, rcol, gcol, bcol, 
                pointstyle, hyperlinks = NULL, local_scale = 1)
        }        
        else if (length(pointstyle) >= length(unique(numlabels))) {
        
            if (length(labels)) {
                stylevec <- c("s", "b", "t")
                curstyle <- stylevec[numlabels[j]]
                .lg3dpoints(x, y, z, filename, rcol, gcol, bcol, 
                  curstyle, hyperlinks = NULL, local_scale = 1)
            }
            else {
            
                .lg3dpoints(x, y, z, filename, rcol, gcol, bcol, 
                  pointstyle[1], hyperlinks = NULL, local_scale = 1)
            }
        }
        else {
        
            .lg3dpoints(x, y, z, filename, rcol, gcol, bcol, 
                pointstyle[1], hyperlinks = NULL, local_scale = 1)
        }
        if (length(metalabels) >= j) {
            
            frcol <- (col2rgb(col.metalab)/255)[1]
            fgcol <- (col2rgb(col.metalab)/255)[2]
            fbcol <- (col2rgb(col.metalab)/255)[3]
            
            # add metalabels and hyperlinks
            write(paste("Text[StyleForm[\"", metalabels[j], "\",FontFamily -> \"Arial\", FontWeight->\"Bold\", FontColor -> RGBColor[", 
                frcol, ",", fgcol, ",", fbcol, "], FontSize->9,", 
                sep = ""), file = filename, append = TRUE)
                
            if (length(hyperlinks) >= j) {
                write(paste("URL -> \"", hyperlinks[j], ",target=_blank\",", 
                  sep = ""), file = filename, append = TRUE)
                  
            }
            
            write(paste("],{", x, ",", y, ",", z, "}],", sep = ""), 
                file = filename, append = TRUE)
                
        }
    }
    
    # create the output files, if this is not already
    # done by a higher-level lg3d.close-function    
    if (!exists(".vrmlgenEnv")) {
    
        write(paste("\n}, Boxed -> False, Axes -> False, AmbientLight->GrayLevel[", 
            ambientlight, "], TextStyle -> {FontFamily -> \"TimesRoman\", FontSlant ->\"Italic\", FontSize -> ", 
            14 * cex.lab, "}, Lighting -> True, BoxRatios -> Automatic, PlotRange -> All ]\n", 
            sep = ""), file = filename, append = TRUE)
            
        # create HTML output    
        if (!is.null(htmlout)) {
            cat("<HTML>", file = htmlout, append = FALSE)
            cat("<HEAD><TITLE>VRMLGen-visualization</TITLE></HEAD><BODY>", 
                file = htmlout, append = TRUE)
            cat(paste("<APPLET ARCHIVE=\"live.jar\" CODE=\"Live.class\" WIDTH=", 
                hwidth, " HEIGHT=", hheight, " ALIGN=LEFT>", 
                sep = ""), file = htmlout, append = TRUE)
            coln <- col2rgb(col.bg)
            cat(paste("<PARAM NAME=\"BGCOLOR\" VALUE=\"", rgb(red = coln[1, 
                ]/255, green = coln[2, ]/255, blue = coln[3, 
                ]/255), "\">", sep = ""), file = htmlout, append = TRUE)
            cat("<PARAM NAME=\"MAGNIFICATION\" VALUE=1.0>", file = htmlout, 
                append = TRUE)
            cat(paste("<PARAM NAME=\"INPUT_FILE\" VALUE=\"", 
                filename, "\">", sep = ""), file = htmlout, append = TRUE)
            cat("</APPLET></BODY>", file = htmlout, append = TRUE)
            cat("</HTML>", file = htmlout, append = TRUE)
        }
        
        # copy the Livegraphics3D live.jar to the output directory
				datadir <- system.file("extdata", package = "vrmlgen")    	
       	curdir <- getwd()
       	
        if (data.class(result<-try(file.copy(file.path(datadir, "live.jar"), file.path(curdir, "live.jar")), TRUE))=="try-error")
        {
          warning("\nCannot copy file live.jar from vrlmgen-folder to current directory. You might need to copy the file manually.")
        }
        
        cat(paste("\nOutput file \"", filename, "\" was generated in folder ", 
            getwd(), ".\n", sep = ""))
    }
    else {
    		# return to current directory, if higher-level plotting
        # functions wrote re-directed the output to a temporary directory
        setwd(curdir)
    }
}


# lower-level function for plotting data points in the Livegraphics3D-format
`.lg3dpoints` <-
function (x, y, z, filename, rcol, gcol, bcol, pointstyle, hyperlinks = NULL, 
    global_scale = 1, local_scale = 1) 
{
    
    # check the point-style and draw the data points at the
    # corresponding coordinates
    
    if (pointstyle == "s") {
        
        for (i in 1:length(x))
        		write(paste("PointSize[", 0.02 * 
            local_scale, "], RGBColor[", rcol, ",", 
            gcol, ",", bcol, "],  Point[{", x[i], 
            ",", y[i], ",", z[i], 
            "}],", sep = ""), file = filename, append = TRUE)
    }
    
    else if (pointstyle == "t") {
        
        for (i in 1:length(x))
        		write(paste("{ SurfaceColor[RGBColor[", 
            rcol, ",", gcol, ",", bcol, "]], Polygon[{{", local_scale * 
                x[i], ", ", local_scale * y[i], ",", local_scale * 
                (z[i] + 0.06415), "}, {", local_scale * global_scale * 
                (x[i] + 0.09072222), ",", local_scale * global_scale * 
                y[i], ",", local_scale * global_scale * (z[i] - 
                0.096225), "}, {", local_scale * global_scale * 
                x[i] - 0.09072222, ",", local_scale * global_scale * 
                y[i] + 0.15713333, ",", local_scale * global_scale * 
                z[i] - 0.096225, "}}], Polygon[{{", local_scale * 
                global_scale * x[i], ",", local_scale * global_scale * 
                y[i], ",", local_scale * global_scale * (z[i] + 
                0.06415), "}, {", local_scale * global_scale * 
                (x[i] - 0.09072222), ",", local_scale * global_scale * 
                (y[i] + 0.15713333), ",", local_scale * global_scale * 
                (z[i] - 0.096225), "}, {", local_scale * global_scale * 
                (x[i] - 0.09072222), ",", local_scale * global_scale * 
                (y[i] - 0.15713333), ",", local_scale * global_scale * 
                (z[i] - 0.096225), "}}], Polygon[{{", local_scale * 
                global_scale * x[i], ",", local_scale * global_scale * 
                y[i], ",", local_scale * global_scale * (z[i] + 
                0.06415), "}, {", local_scale * global_scale * 
                (x[i] - 0.09072222), ",", local_scale * global_scale * 
                (y[i] - 0.15713333), ",", local_scale * global_scale * 
                (z[i] - 0.096225), "}, {", local_scale * global_scale * 
                (x[i] + 0.09072222), ",", local_scale * global_scale * 
                y[i], ",", local_scale * global_scale * (z[i] - 
                0.096225), "}}], Polygon[{{", local_scale * global_scale * 
                (x[i] + 0.09072222), ",", local_scale * global_scale * 
                y[i], ",", local_scale * global_scale * (z[i] - 
                0.096225), "}, {", local_scale * global_scale * 
                (x[i] - 0.09072222), ",", local_scale * global_scale * 
                (y[i] - 0.15713333), ",", local_scale * global_scale * 
                (z[i] - 0.096225), "}, {", local_scale * global_scale * 
                (x[i] - 0.09072222), ",", local_scale * global_scale * 
                (y[i] + 0.15713333), ",", local_scale * global_scale * 
                (z[i] - 0.096225), "}}] },\n", sep = ""), file = filename, 
            append = TRUE)
    }
    else {
    
        for (i in 1:length(x))
        		write(paste("SurfaceColor[RGBColor[", 
            rcol, ",", gcol, ",", bcol, "]], Cuboid[{", local_scale * 
                (x[i] - 0.04), ",", local_scale * (y[i] - 0.04), 
            ",", local_scale * (z[i] - 0.04), "},{", local_scale * 
                global_scale * (x[i] + 0.04), ",", local_scale * 
                global_scale * (y[i] + 0.04), ",", local_scale * 
                global_scale * (z[i] + 0.04), "}],", sep = ""), 
            file = filename, append = TRUE)
    }
    
    # add hyperlinks to the datapoints
    if (!is.null(hyperlinks)) {
        
        for (i in 1:length(hyperlinks))
        		write(paste("Text[StyleForm[\" \", URL -> \"", 
            hyperlinks, "\"], {", local_scale * x[i], ",", local_scale * 
                y[i], ",", local_scale * z[i], "}],", sep = ""), 
            file = filename, append = TRUE)
    }
    
}


# lower-level function for adding text to a 3D-scene
# in the Livegraphics3D-format
`.lg3dtext` <-
function (x, y, z, text, filename, col, scale, fontfamily, fontweight, 
    hyperlink) 
{
    
    # convert the color to RGB-format
    rcol <- (col2rgb(col)/255)[1]
    gcol <- (col2rgb(col)/255)[2]
    bcol <- (col2rgb(col)/255)[3]
    
    # add hyperlinks to the text
    if (is.null(hyperlink)) {
        
        for (i in 1:length(x))
        		write(paste("RGBColor[", rcol, 
            ",", gcol, ",", bcol, "], Text [StyleForm[\"", text[i], 
            "\", FontFamily -> \"", fontfamily, "\", FontSize-> ", 
            14 * scale, ", FontWeight->\"", fontweight, "\"], {", 
            x[i], ",", y[i], ",", z[i], "}],\n", sep = ""), file = filename, 
            append = TRUE)
    }
    else {
    
        for (i in 1:length(x))
        		write(paste("RGBColor[", rcol, 
            ",", gcol, ",", bcol, "], Text [StyleForm[\"", text[i], 
            "\", FontFamily -> \"", fontfamily, "\", FontSize-> ", 
            14 * scale, ", URL -> \"", hyperlink, "\", FontWeight->\"", 
            fontweight, "\"], {", x[i], ",", y[i], ",", z[i], 
            "}],\n", sep = ""), file = filename, append = TRUE)
    }
}


# creating 3D-meshes in Livegraphics3D-format
`.lmesh` <-
function (xfun = "sin(v)*cos(u)", yfun = "sin(v)*sin(u)", zfun = "cos(v)", 
    param1 = "u", param2 = "v", range1 = c(0, 2 * pi), range2 = c(0, 
        pi), size1 = 30, size2 = 30, filename = "out.m", cols = "red", 
    scalefac = 4, autoscale = "independent", lab.axis = c("X-axis", 
        "Y-axis", "Z-axis"), col.axis = "black", showaxis = TRUE, 
    col.lab = "black", col.bg = "white", cex.lab = 1, ambientlight = 0.5, 
    htmlout = NULL, hwidth = 1200, hheight = 800) 
{
    
    # check if the required input data is available
    if (is.null(xfun) || is.null(yfun) || is.null(zfun)) {
        stop("Either the paramater infile or data or the parameters xfun, yfun and zfun have to be specified.")
    }
    
    if (is.null(param1) || is.null(param2)) {
        stop("The parameter names param1 and param2 have not been specified")
    }
    
    # initialize variables
    x <- NULL
    y <- NULL
    z <- NULL
    smin <- range1[1]
    smax <- range1[2]
    tmin <- range2[1]
    tmax <- range2[2]
    sn <- size1
    tn <- size2
    ds <- (smax - smin)/sn
    dt <- (tmax - tmin)/tn
    
    # compute the data points representing the parametric function
    for (i in seq(smin, (smax - ds/2), ds)) {
        
        for (j in seq(tmin, (tmax - dt/2), dt)) {
        
            eval(parse(text = paste(param1, " <- ", i)))
            eval(parse(text = paste(param2, " <- ", j)))
            x <- c(x, eval(parse(text = xfun)))
            y <- c(y, eval(parse(text = yfun)))
            z <- c(z, eval(parse(text = zfun)))
            eval(parse(text = paste(param1, " <- ", param1, " + ds")))
            x <- c(x, eval(parse(text = xfun)))
            y <- c(y, eval(parse(text = yfun)))
            z <- c(z, eval(parse(text = zfun)))
            eval(parse(text = paste(param2, " <- ", param2, " + dt")))
            x <- c(x, eval(parse(text = xfun)))
            y <- c(y, eval(parse(text = yfun)))
            z <- c(z, eval(parse(text = zfun)))
            eval(parse(text = paste(param1, " <- ", param1, " - ds")))
            x <- c(x, eval(parse(text = xfun)))
            y <- c(y, eval(parse(text = yfun)))
            z <- c(z, eval(parse(text = zfun)))
        }
    }
    
    # combine data
    data <- as.matrix(cbind(x, y, z))
    if (ncol(data) != 3) {
        stop("Data matrix does not have 3 columns!")
    }
    

    # verify if .lmesh is called within a vrmlgen-environment
    # initialized with lg3d.open()
    curdir <- NULL
    if (!exists(".vrmlgenEnv")) {
        
        # set the general 3D scene settings, if they are
        # not provided by a higher-level environment
        
        # write LiveGraphics3D header
        write("Graphics3D[\n{\n", file = filename, append = FALSE)
        
        # set background color
        bg_rcol <- (col2rgb(col.bg)/255)[1]
        bg_gcol <- (col2rgb(col.bg)/255)[2]
        bg_bcol <- (col2rgb(col.bg)/255)[3]
        
        write(paste("Background {\n\t skyColor [\n\t\t ", bg_rcol, 
            bg_gcol, bg_bcol, " \n\t]\n}", sep = " "), file = filename, 
            append = TRUE)
    }
    else {
    
    		# obtain the general settings from the
    		# higher level environment
        vrmlgenEnv <- get(".vrmlgenEnv",envir=.GlobalEnv)
        filename <- vrmlgenEnv$filename
        htmlout <- vrmlgenEnv$html
        hheight <- vrmlgenEnv$hheight
        hwidth <- vrmlgenEnv$hwidth
        curdir <- getwd()
        setwd(vrmlgenEnv$VRMLDir)
    }
    
    # apply the chosen scaling method
    scaledat <- function(data) {
        diff <- max(data) - min(data)
        if (diff == 0) 
            return(data)
        else return(scalefac * (data - min(data))/(diff))
    }
    
    if (autoscale == "independent") 
        data <- apply(data, 2, scaledat)
    else if (autoscale == "equidist") 
        data <- scalefac * (data/max(data))
    else if (autoscale == "equicenter") 
        data <- scalefac/2 * data/max(data) + scalefac/2

		        
    lab_rcol <- NULL
    lab_gcol <- NULL
    lab_bcol <- NULL
    ax_rcol <- NULL
    ax_gcol <- NULL
    ax_bcol <- NULL
    
    # draw the axes
    if (showaxis) {
        axis3d(lab.axis, filename, type = "lg3d", col.lab, col.axis, 
            cex.lab, local_scale = scalefac, global_scale = 1)
    }
    
    for (j in 1:nrow(data)) {
    
        x <- data[j, 1]
        y <- data[j, 2]
        z <- data[j, 3]
        
        # apply automatic color selection
        colsel <- .colorselection(type = "mesh", j, cols, nrow(data),"", "")
            
        rcol <- colsel$rcol
        gcol <- colsel$gcol
        bcol <- colsel$bcol
        
        if (j%%4 == 1) 
            write(paste(" SurfaceColor[RGBColor[", rcol, ",", 
                gcol, ",", bcol, "]], Polygon[{{", x, ",", y, 
                ",", z, "},", sep = ""), file = filename, append = TRUE)
        else if (j%%4 == 0) 
            write(paste(" {", x, ",", y, ",", z, "}}],", sep = ""), 
                file = filename, append = TRUE)
        else write(paste(" {", x, ",", y, ",", z, "},", sep = ""), 
            file = filename, append = TRUE)
            
    }
    
    # configure general plot settings (font size, text style, etc.)
    write(paste("\n}, Boxed -> False, Axes -> False, AmbientLight->GrayLevel[", 
        ambientlight, "], TextStyle -> {FontFamily -> \"TimesRoman\", FontSlant ->\"Italic\", FontSize -> ", 
        14 * cex.lab, "}, Lighting -> True, BoxRatios -> Automatic, PlotRange -> All ]\n", 
        sep = ""), file = filename, append = TRUE)

		
		# create the output files, if this is not already
    # done by a higher-level lg3d.close-function    
    if (!exists(".vrmlgenEnv")) {
        
        # create HTML output
        if (!is.null(htmlout)) {
            
            cat("<HTML>", file = htmlout, append = FALSE)
            cat("<HEAD><TITLE>VRMLGen-visualization</TITLE></HEAD><BODY>", 
                file = htmlout, append = TRUE)
            cat(paste("<APPLET ARCHIVE=\"live.jar\" CODE=\"Live.class\" WIDTH=", 
                hwidth, " HEIGHT=", hheight, " ALIGN=LEFT>", 
                sep = ""), file = htmlout, append = TRUE)
            coln <- col2rgb(col.bg)
            cat(paste("<PARAM NAME=\"BGCOLOR\" VALUE=\"", rgb(red = coln[1, 
                ]/255, green = coln[2, ]/255, blue = coln[3, 
                ]/255), "\">", sep = ""), file = htmlout, append = TRUE)
            cat("<PARAM NAME=\"MAGNIFICATION\" VALUE=1.0>", file = htmlout, 
                append = TRUE)
            cat(paste("<PARAM NAME=\"INPUT_FILE\" VALUE=\"", 
                filename, "\">", sep = ""), file = htmlout, append = TRUE)
            cat("</APPLET>", file = htmlout, append = TRUE)
            cat("</HTML>", file = htmlout, append = TRUE)
        }

				
       	# copy the Livegraphics3D live.jar to the output directory
       	datadir <- system.file("extdata", package = "vrmlgen")
       	curdir <- getwd()
       	
        if (data.class(result<-try(file.copy(file.path(datadir, "live.jar"), file.path(curdir, "live.jar")), TRUE))=="try-error")
        {
          warning("\nCannot copy file live.jar from vrlmgen-folder to current directory. You might need to copy the file manually.")
        }
        
        # show success message        
        cat(paste("\nOutput file \"", filename, "\" was generated in folder ", 
            getwd(), ".\n", sep = ""))
    }
    else {
    		# return to current directory, if higher-level plotting
        # functions wrote re-directed the output to a temporary directory
        setwd(curdir)
    }
}


# helper function to combine rotation vectors
# by quaternion multiplication
`.multquaternion` <-
function (q1, q2) 
{
    tmp <- numeric(4)
    tmp[1] = q2[4] * q1[1] + q2[1] * q1[4] + q2[2] * q1[3] - 
        q2[3] * q1[2]
    tmp[2] = q2[4] * q1[2] + q2[2] * q1[4] + q2[3] * q1[1] - 
        q2[1] * q1[3]
    tmp[3] = q2[4] * q1[3] + q2[3] * q1[4] + q2[1] * q1[2] - 
        q2[2] * q1[1]
    tmp[4] = q2[4] * q1[4] - q2[1] * q1[1] - q2[2] * q1[2] - 
        q2[3] * q1[3]
    return(tmp)
}


# helper function to represent rotation vectors
# by quaternions
`.quaternion2rot` <-
function (quat) 
{
    res <- numeric(4)
    res[4] = acos(quat[4]) * 2
    for (i in 1:3) res[i] = quat[i]/sin(res[4]/2)
    return(res)
}


# helper function to convert quaternion representations
# of rotation vectors back to the normal representation
`.rot2quaternion` <-
function (rot = c(1, 0, 0, pi)) 
{
    res <- numeric(4)
    for (i in 1:3) res[i] <- rot[i] * sin(rot[4]/2)
    res[4] <- cos(rot[4]/2)
    return(res)
}


# creating bar plots in VRML-format
`.vbar` <-
function (data, row.labels = rownames(data), col.labels = colnames(data), 
    metalabels = NULL, filename = "out.wrl", space = 0.5, cols = rainbow(length(as.matrix(data))), 
    rcols = NULL, ccols = NULL, origin = c(0, 0, 0), scalefac = 4, autoscale = TRUE, 
    ignore_zeros = TRUE, lab.axis = c("X-axis", "Y-axis", "Z-axis"), lab.vertical = FALSE, 
    col.axis = "black", showaxis = TRUE, col.lab = "black", col.bg = "white", 
    cex.lab = 1, cex.rowlab = 1, cex.collab = 1, navigation = "EXAMINE", 
    fov = 0.785, pos = rep(scalefac + 4, 3), dir = c(0.19, 0.45, 
        0.87, 2.45), transparency = 0, htmlout = NULL, hwidth = 1200, 
    hheight = 800, showlegend = TRUE) 
{
    
    data <- as.matrix(data)
    
    
    # verify if .vbar is called within a vrmlgen-environment
    # initialized with vrml.open()
    
    scale <- 1
    curdir <- NULL
    if (!exists(".vrmlgenEnv")) {
        
        # set the general 3D scene settings, if they are
        # not provided by a higher-level environment
        
        # write VRML header        
        write("#VRML V2.0 utf8\n", file = filename, append = FALSE)
        write(paste("\nViewpoint {\n\tfieldOfView", fov, "\n\tposition", 
            pos[1], pos[2], pos[3], "\n\torientation", dir[1], 
            dir[2], dir[3], dir[4], "\n\tjump TRUE\n\tdescription \"viewpoint1\"\n}\n", 
            sep = " "), file = filename, append = TRUE)
            
        write(paste("\nNavigationInfo { type \"", navigation, 
            "\" }\n", sep = ""), file = filename, append = TRUE)
        
        # set background color
        bg_rcol <- (col2rgb(col.bg)/255)[1]
        bg_gcol <- (col2rgb(col.bg)/255)[2]
        bg_bcol <- (col2rgb(col.bg)/255)[3]
        write(paste("Background {\n\t skyColor [\n\t\t ", bg_rcol, 
            bg_gcol, bg_bcol, " \n\t]\n}", sep = " "), file = filename, 
            append = TRUE)
    }
    else {
    
        vrmlgenEnv <- get(".vrmlgenEnv",envir=.GlobalEnv)
        filename <- vrmlgenEnv$filename
        htmlout <- vrmlgenEnv$html
        hheight <- vrmlgenEnv$hheight
        hwidth <- vrmlgenEnv$hwidth
        scale <- vrmlgenEnv$scale
        curdir <- getwd()
        setwd(vrmlgenEnv$VRMLDir)
    }
    
    # scale the data to fit within the axes ranges
    if(autoscale)
    	data <- scale * scalefac * (data/max(data))
    
    # configure popup-text for metalabels
    popuptext <- TRUE
    if (is.null(metalabels)) {
        popuptext <- FALSE
    }
    
    # generate popup-text
    if (popuptext) 
        write(paste("\n\nPROTO PopupText [\n\teventIn MFString showThis\n\texposedField SFNode fontStyle FontStyle {\n\t\t\t\tjustify [ \"MIDDLE\", \"MIDDLE\" ] }\n\texposedField SFNode appearance Appearance { material Material { } } \n\tfield SFString default \"\" \n\tfield MFString stringlist [ \"\" ] \n", 
            paste(paste("eventIn SFBool over", 0:(length(metalabels) - 
                1), "\n\t", sep = "", collapse = ""), sep = "", 
                collapse = ""), "\t] { \nGroup { children [ \n\tDEF POPIT Script { \n\t\tfield SFString defstring IS default \n\t\tfield MFString list IS stringlist \n\t\tfield MFString strout [ \"\" ] ", 
            paste(paste("\n eventIn SFBool over", 0:(length(metalabels) - 
                1), " IS over", 0:(length(metalabels) - 1), sep = "", 
                collapse = "")), "\n\t\teventOut MFString string_changed \n\t\turl [ \"javascript: \n\t\tfunction evtnum(num, value) \n\t\t{ \n\t\tif (value && (num < list.length)) \n\t\t\tstrout[0] = list[num]; \n\t\telse \n\t\t\tstrout[0] = defstring; \n\t\tstring_changed = strout; \n\t\t} \n \n\t", 
            paste(paste("\tfunction over", 0:(length(metalabels) - 
                1), "(v, t) { evtnum(", 0:(length(metalabels) - 
                1), ", v); } \n", sep = "", collapse = "")), 
            "\t\", \n\t\t\"Popup.class\"] \n\t} \n\t \n\tTransform { \n\ttranslation 0 0 ", 
            scale * scalefac, " \n\trotation 0.28108 0.67860 0.67860 2.59356\n\t \n\tchildren Shape { \n\t\tappearance IS appearance \n\tgeometry DEF POPUP Text { \n\t\tstring \"\" \n\t\tset_string IS showThis \n\t\tfontStyle IS fontStyle \n\t\t} \n\t} \n\t} \n] } \nROUTE POPIT.string_changed TO POPUP.set_string \n} \n \n \nGroup { \nchildren DEF POP PopupText { \n\t\t#default \"Nothing selected\" \n\t\tstringlist [ ", 
            paste(paste("\"", metalabels[1:(length(metalabels) - 
                1)], "\",", sep = "", collapse = " ")), "\"", 
            metalabels[length(metalabels)], "\" ] \n\t} \n} \n\n", 
            sep = "", collapse = ""), file = filename, append = TRUE)
    
    # draw the axes
    if (showaxis) {
        axis3d(lab.axis, filename, type = "vrml", col.lab, col.axis, 
            cex.lab, local_scale = scalefac, global_scale = 1)
    }
    
    # calculate bar-length and -width (account for scaling and spacing)
    blength <- scale * scalefac/(nrow(data) * (1 + space))
    bwidth <- scale * scalefac/(ncol(data) * (1 + space))
    
    # rotate the data labels to obtain horizontal or vertical labels
    rot_vec <- c(0, 0.70711, 0.70711, 3.14159)
    if (lab.vertical) 
        rot_vec <- c(0.57735, 0.57735, 0.57735, 4.18879)
    
    # draw row-labels
    if (!is.null(row.labels) && showlegend) {
        
        for (j in 1:length(row.labels)) {
        
        		# apply automatic color selection
            colsel <- .colorselection(type = "bar", j, rcols, 
                1, row.labels, "")
                
            rcol <- colsel$rcol
            gcol <- colsel$gcol
            bcol <- colsel$bcol
            
            cur_xwidth <- j/nrow(data) * scale * scalefac - blength/2
            write(paste("Transform {\n\ttranslation ", cur_xwidth, 
                -0.4, scale * scalefac + 0.2, "\n\tscale ", cex.rowlab * 
                  0.28, cex.rowlab * 0.28, cex.rowlab * 0.28, 
                "\n\trotation ", rot_vec[1], rot_vec[2], rot_vec[3], 
                rot_vec[4], "\n\tchildren Shape {\n\t\tappearance Appearance { material Material {diffuseColor ", 
                rcol, " ", gcol, " ", bcol, "  } }\n\t\tgeometry Text { string \"", 
                row.labels[j], "\" }\n\t}\n}", sep = " "), file = filename, 
                append = TRUE)
        }
    }
    
    # rotate the data labels to obtain horizontal or vertical labels
    rot_vec <- c(0.57735, 0.57735, 0.57735, 2.0944)
    if (lab.vertical) 
        rot_vec <- c(0.707107, 0, 0.707107, 3.14159)
        
    # draw column-labels
    if (!is.null(col.labels) && showlegend) {
    
        for (j in 1:length(col.labels)) {
        
        		# apply automatic color selection
            colsel <- .colorselection(type = "bar", j, ccols,  
                1, col.labels, "")
            rcol <- colsel$rcol
            gcol <- colsel$gcol
            bcol <- colsel$bcol
            
            cur_ywidth <- j/ncol(data) * scale * scalefac - bwidth/2
            write(paste("Transform {\n\ttranslation ", -0.4, 
                cur_ywidth, scale * scalefac + 0.2, "\n\tscale ", 
                cex.collab * 0.28, cex.collab * 0.28, cex.collab * 
                  0.28, "\n\trotation ", rot_vec[1], rot_vec[2], 
                rot_vec[3], rot_vec[4], "\n\tchildren Shape {\n\t\tappearance Appearance { material Material {diffuseColor ", 
                rcol, " ", gcol, " ", bcol, "  } }\n\t\tgeometry Text { string \"", 
                col.labels[j], "\" }\n\t}\n}", sep = " "), file = filename, 
                append = TRUE)
        }
    }
    
    # initialize color and meta-label variables
    popup_txt_str <- ""
    popup_txt_str2 <- ""
    
    rcol <- NULL
    gcol <- NULL
    bcol <- NULL
    
    # set RGB-colors
    if (length(cols) >= (ncol(data) * nrow(data))) {
        rcol <- sapply(cols, function(x) (col2rgb(x)/255)[1])
        gcol <- sapply(cols, function(x) (col2rgb(x)/255)[2])
        bcol <- sapply(cols, function(x) (col2rgb(x)/255)[3])
    }
    else {    		
        if (length(cols) >= 1) {
            rcol <- (col2rgb(cols[1])/255)[1]
            gcol <- (col2rgb(cols[1])/255)[2]
            bcol <- (col2rgb(cols[1])/255)[3]
        }
        else {
            rcol <- (col2rgb("lightblue")/255)[1]
            gcol <- (col2rgb("lightblue")/255)[2]
            bcol <- (col2rgb("lightblue")/255)[3]
        }
    }
    
    # main loop: iterate over data points
    counter <- 0
    for (j in 1:nrow(data)) {
        
        for (k in 1:ncol(data)) {
        
            x <- j/nrow(data) * scale * scalefac
            y <- k/ncol(data) * scale * scalefac
            z <- data[j, k]
            
            if(ignore_zeros && (z == 0))
  						next
                        
            counter <- counter + 1
                        
            # set current color
            if (!is.null(ccols)) {
                rcol <- (col2rgb(ccols[k])/255)[1]
                gcol <- (col2rgb(ccols[k])/255)[2]
                bcol <- (col2rgb(ccols[k])/255)[3]
            }
            else if (!is.null(rcols)) {
                rcol <- (col2rgb(rcols[j])/255)[1]
                gcol <- (col2rgb(rcols[j])/255)[2]
                bcol <- (col2rgb(rcols[j])/255)[3]
            }
            else {
                if(length(rcol) > 1)
                {
	                rcol <- rcol[(k - 1) * ncol(data) + 
	                  j]
	                gcol <- gcol[(k - 1) * ncol(data) + 
	                  j]
	                bcol <- bcol[(k - 1) * ncol(data) + 
	                  j]
	              }
            }
            
            
            # set meta-label touch sensor
            if (popuptext) {
                popup_txt_str <- paste("Group {\n  children [\n    DEF Schalter", 
                  (k - 1) * ncol(data) + j - 1, " TouchSensor { }", 
                  sep = "", collapse = "")
                popup_txt_str2 <- "\n]\n}"
            }
            
            # draw data point
            write(paste(popup_txt_str, "Transform {\n\ttranslation ", 
                origin[1] + x, origin[2] + y, origin[3] + z/2, 
                "\n\tchildren Shape {\n\t\tappearance Appearance { material Material { diffuseColor ", 
                rcol, gcol, bcol, " } }\n\tgeometry Box { size ", 
                blength, bwidth, z, "}\n\t\t}\n}", popup_txt_str2, 
                sep = " "), file = filename, append = TRUE)
        }
    }
    write("\n", file = filename, append = TRUE)
    
    
    # configure interactive popup-text
    if (popuptext) 
        write(paste("ROUTE Schalter", 0:(length(metalabels) - 
            1), ".isOver TO POP.over", 0:(length(metalabels) - 
            1), "\n", sep = "", collapse = ""), file = filename, 
            append = TRUE)
            
                
    # create the output files, if this is not already
    # done by a higher-level lg3d.close-function    
    if (!exists(".vrmlgenEnv")) {
        
     		# create HTML output            
		    if (!is.null(htmlout)) {
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
        cat(paste("\nOutput file \"", filename, "\" was generated in folder ", 
            getwd(), ".\n", sep = ""))
    }
    else {
    		# return to current directory, if higher-level plotting
        # functions wrote re-directed the output to a temporary directory
        setwd(curdir)
    }
}


# creating scatter plots in VRML-format
`.vcloud` <-
function (x, y = NULL, z = NULL, labels = rownames(data), metalabels = NULL, 
    hyperlinks = NULL, filename = "out.wrl", pointstyle = c("s", 
        "b", "c"), cols = rainbow(length(unique(labels))), showdensity = FALSE, 
    scalefac = 4, autoscale = "independent", lab.axis = c("X-axis", 
        "Y-axis", "Z-axis"), col.axis = "black", showaxis = TRUE, 
    col.lab = "black", col.bg = "white", cex.lab = 1, navigation = "EXAMINE", 
    transparency = 0, fov = 0.785, pos = rep(scalefac + 4, 3), 
    dir = c(0.19, 0.45, 0.87, 2.45), htmlout = NULL, hwidth = 1200, 
    hheight = 800, showlegend = TRUE) 
{
    
    # extract the coordinate vectors from the input        
    xyz_parse <- xyz.coords(x, y, z)
    
    data <- cbind(xyz_parse$x, xyz_parse$y, xyz_parse$z)
    if (ncol(data) != 3) {
        stop("Data matrix does not have 3 columns!")
    }
    
    # identify the unique set of row labels
    numlabels <- NULL
    if (length(labels)) {
    
        lab <- unique(unlist(labels))
        numlabels <- apply(as.matrix(labels), 1, function(x) match(x, 
            as.matrix(lab)))
    }
    
   
    # verify if .vcloud is called within a vrmlgen-environment
    # initialized with lg3d.open()
    curdir <- NULL
    scale <- 1
    if (!exists(".vrmlgenEnv")) {
    
    		# set the general 3D scene settings, if they are
        # not provided by a higher-level environment
        
        # write VRML header        
        write("#VRML V2.0 utf8\n", file = filename, append = FALSE)
        
        write(paste("\nViewpoint {\n\tfieldOfView", fov, "\n\tposition", 
            pos[1], pos[2], pos[3], "\n\torientation", dir[1], 
            dir[2], dir[3], dir[4], "\n\tjump TRUE\n\tdescription \"viewpoint1\"\n}\n", 
            sep = " "), file = filename, append = TRUE)
            
        write(paste("\nNavigationInfo { type \"", navigation, 
            "\" }\n", sep = ""), file = filename, append = TRUE)
            
        # set background color
        bg_rcol <- (col2rgb(col.bg)/255)[1]
        bg_gcol <- (col2rgb(col.bg)/255)[2]
        bg_bcol <- (col2rgb(col.bg)/255)[3]
        write(paste("Background {\n\t skyColor [\n\t\t ", bg_rcol, 
            bg_gcol, bg_bcol, " \n\t]\n}", sep = " "), file = filename, 
            append = TRUE)
    }
    else {
        
        # obtain the general settings from the
    		# higher level environment
        vrmlgenEnv <- get(".vrmlgenEnv",envir=.GlobalEnv)
        filename <- vrmlgenEnv$filename
        htmlout <- vrmlgenEnv$html
        VRMLDir <- vrmlgenEnv$VRMLDir
        hheight <- vrmlgenEnv$hheight
        hwidth <- vrmlgenEnv$hwidth
        scale <- vrmlgenEnv$scale
        curdir <- getwd()
        setwd(vrmlgenEnv$VRMLDir)
    }
    
    # apply the chosen scaling method
    scaledat <- function(data) {
        diff <- max(data) - min(data)
        if (diff == 0) 
            return(data)
        else return(scale * scalefac * (data - min(data))/(diff))
    }
    
    if (autoscale == "independent") 
        data <- apply(data, 2, scaledat)
    else if (autoscale == "equidist") 
        data <- scale * scalefac * (data/max(data))
    
        
     # configure popup-text for metalabels
    popuptext <- TRUE
    if (is.null(metalabels)) {
        popuptext <- FALSE
    }
    
    if (popuptext) 
        write(paste("\n\nPROTO PopupText [\n\teventIn MFString showThis\n\texposedField SFNode fontStyle FontStyle {\n\t\t\t\tjustify [ \"MIDDLE\", \"MIDDLE\" ] }\n\texposedField SFNode appearance Appearance { material Material { } } \n\tfield SFString default \"\" \n\tfield MFString stringlist [ \"\" ] \n", 
            paste(paste("eventIn SFBool over", 0:(nrow(data) - 
                1), "\n\t", sep = "", collapse = ""), sep = "", 
                collapse = ""), "\t] { \nGroup { children [ \n\tDEF POPIT Script { \n\t\tfield SFString defstring IS default \n\t\tfield MFString list IS stringlist \n\t\tfield MFString strout [ \"\" ] ", 
            paste(paste("\n eventIn SFBool over", 0:(nrow(data) - 
                1), " IS over", 0:(nrow(data) - 1), sep = "", 
                collapse = "")), "\n\t\teventOut MFString string_changed \n\t\turl [ \"javascript: \n\t\tfunction evtnum(num, value) \n\t\t{ \n\t\tif (value && (num < list.length)) \n\t\t\tstrout[0] = list[num]; \n\t\telse \n\t\t\tstrout[0] = defstring; \n\t\tstring_changed = strout; \n\t\t} \n \n\t", 
            paste(paste("\tfunction over", 0:(nrow(data) - 1), 
                "(v, t) { evtnum(", 0:(nrow(data) - 1), ", v); } \n", 
                sep = "", collapse = "")), "\t\", \n\t\t\"Popup.class\"] \n\t} \n\t \n\tTransform { \n\ttranslation 0 0 ", 
            scale * scalefac, " \n\trotation 0.28108 0.67860 0.67860 2.59356\n\t \n\tchildren Shape { \n\t\tappearance IS appearance \n\tgeometry DEF POPUP Text { \n\t\tstring \"\" \n\t\tset_string IS showThis \n\t\tfontStyle IS fontStyle \n\t\t} \n\t} \n\t} \n] } \nROUTE POPIT.string_changed TO POPUP.set_string \n} \n \n \nGroup { \nchildren DEF POP PopupText { \n\t\t#default \"Nothing selected\" \n\t\tstringlist [ ", 
            paste(paste("\"", metalabels[1:(nrow(data) - 1)], 
                "\",", sep = "", collapse = " ")), "\"", metalabels[nrow(data)], 
            "\" ] \n\t} \n} \n\n"), file = filename, append = TRUE)
    
    
    # draw the plot legend
    if (length(labels) && showlegend) {
        cur_height <- scale * scalefac + 1.2
        for (j in 1:length(unique(labels))) {
            .vrmltext(0, 0, cur_height, unique(labels)[j], filename, 
                cols[unique(numlabels)[j]], 1, c(0.28108, 0.6786, 
                  0.6786, 2.59356), "SANS", "NONE", NULL)
            cur_height <- cur_height + 0.4
        }
    }
    
    # draw the axes
    if (showaxis) {
    
    		# call the axis3d helper function to draw the axes
        axis3d(lab.axis, filename, type = "vrml", col.lab, col.axis, 
            cex.lab, local_scale = scalefac, global_scale = 1)
    }
    
    # initialize color and meta-label variables
    popup_txt_str <- ""
    popup_txt_str2 <- ""
    hyperlink <- NULL
    
    # main loop: iterate over data points
    for (j in 1:nrow(data)) {
    
    
        x <- data[j, 1]
        y <- data[j, 2]
        z <- data[j, 3]
                
        # apply automatic color selection
        colsel <- .colorselection(type = "cloud", j, cols, nrow(data), 
            labels, numlabels)
            
        rcol <- colsel$rcol
        gcol <- colsel$gcol
        bcol <- colsel$bcol
        
        if (popuptext) {
            popup_txt_str <- paste("Group {\n  children [\n    DEF Schalter", 
                j - 1, " TouchSensor { }", sep = "", collapse = "")
            popup_txt_str2 <- "\n]\n}"
        }
        
        if (!is.null(hyperlinks)) {
            hyperlink <- hyperlinks[j]
        }
        
        # draw data point with given point-style
        # and popup-text for meta-labels
        if (length(pointstyle) == 1) {
            
            write(popup_txt_str, file = filename, append = TRUE)
            
            # call lower level plotting function to draw datapoints
            .vrmlpoints(x, y, z, filename, rcol, gcol, bcol, 
                pointstyle, hyperlinks = hyperlink, global_scale = 1, 
                local_scale = 1, transparency)
           
            write(popup_txt_str2, file = filename, append = TRUE)
        }        
        else if (length(pointstyle) >= length(unique(numlabels))) {
        
            if (length(labels)) {
            
                stylevec <- c("s", "b", "c")
                curstyle <- stylevec[numlabels[j]]
                
                write(popup_txt_str, file = filename, append = TRUE)
                
                .vrmlpoints(x, y, z, filename, rcol, gcol, bcol, 
                  curstyle, hyperlinks = hyperlink, global_scale = 1, 
                  local_scale = 1, transparency)
                  
                write(popup_txt_str2, file = filename, append = TRUE)
            }
            else {
            
                popup_txt_str <- ""
                popup_txt_str2 <- ""
                
                write(popup_txt_str, file = filename, append = TRUE)
                
                .vrmlpoints(x, y, z, filename, rcol, gcol, bcol, 
                  pointstyle[1], hyperlinks = hyperlink, global_scale = 1, 
                  local_scale = 1, transparency)
                  
                write(popup_txt_str2, file = filename, append = TRUE)
            }
        }
        else {
        
            write(popup_txt_str, file = filename, append = TRUE)
            
            .vrmlpoints(x, y, z, filename, rcol, gcol, bcol, 
                pointstyle[1], hyperlinks = hyperlink, global_scale = 1, 
                local_scale = 1, transparency)
                
            write(popup_txt_str2, file = filename, append = TRUE)
        }
    }    
    write("\n", file = filename, append = TRUE)
    
    # configure interactive popup-text
    if (popuptext) 
        write(paste("ROUTE Schalter", 0:(nrow(data) - 1), ".isOver TO POP.over", 
            0:(nrow(data) - 1), "\n", sep = "", collapse = ""), 
            file = filename, append = TRUE)

		# draw density estimation contour surfaces            
    if (showdensity) 
        est <- .vdense(data, filename)
          
        
    # create the output files, if this is not already
    # done by a higher-level lg3d.close-function
    if (!exists(".vrmlgenEnv")) {
        	            
       # create HTML outputr       
		   if (!is.null(htmlout)) {
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
        cat(paste("\nOutput file \"", filename, "\" was generated in folder ", 
            getwd(), ".\n", sep = ""))
    }
    else {    
        # return to current directory, if higher-level plotting
        # functions wrote re-directed the output to a temporary directory        
        setwd(curdir)
    }
}


# helper function for density estimation contour surface computation
`.vdense` <-
function (x, filename) 
{
    nobs <- nrow(x)
    ndim <- ncol(x)
    weights <- rep(1, nobs)
    
    # compute number of bins (depending on number of observations)
    nbins <- round((nobs > 500) * 8 * log(nobs)/ndim)
    
    # compute number of breaks    
    if (nbins > 0) {
        breaks <- cbind(seq(min(x[, 1]), max(x[, 1]), length = nbins + 
            1), seq(min(x[, 2]), max(x[, 2]), length = nbins + 
            1))
        f1 <- cut(x[, 1], breaks = breaks[, 1])
        f2 <- cut(x[, 2], breaks = breaks[, 2])
        f3 <- cut(x[, 3], breaks = breaks[, 3])
        freq <- table(f1, f2, f3)
        dimnames(freq) <- NULL
        midpoints <- (breaks[-1, ] + breaks[-(nbins + 1), ])/2
        z1 <- midpoints[, 1]
        z2 <- midpoints[, 2]
        z3 <- midpoints[, 3]
        X <- as.matrix(expand.grid(z1, z2, z3))
        X.f <- as.vector(freq)
        id <- (X.f > 0)
        X <- X[id, ]
        dimnames(X) <- list(NULL, dimnames(x)[[2]])
        X.f <- X.f[id]
        bins <- list(x = X, x.freq = X.f, midpoints = midpoints, 
            breaks = breaks, table.freq = freq)
        x <- bins$x
        weights <- bins$x.freq
        nx <- length(bins$x.freq)
    }
    else nx <- nobs
    h <- .hnorm(x, weights)
    rawdata <- list(nbins = nbins, x = x, nobs = nobs, ndim = ndim)
    
    # calculate density estimation contours
    est <- .vdensecalc(x, h = h, weights = weights, rawdata = rawdata, 
        filename = filename)
    return(TRUE)
}


# compute density estimation contour surfaces
`.vdensecalc` <-
function (x, h = .hnorm(x, weights), eval.points = NULL, weights = rep(1, 
    length(x)), rawdata = list(), eval.type = "grid", filename = "test.wrl") 
{
    
    # intial settings
    xlim <- range(x[, 1])
    ylim <- range(x[, 2])
    zlim <- range(x[, 3])
    ngrid <- 20
    if (eval.type != "points") {
        evp <- cbind(seq(xlim[1], xlim[2], length = ngrid), seq(ylim[1], 
            ylim[2], length = ngrid), seq(zlim[1], zlim[2], length = ngrid))
        eval.points <- evp
    }
    h.weights <- rep(1, nrow(x))
    n <- nrow(x)
    nnew <- nrow(eval.points)
    hmult <- 1
    
    # compute density estimation    
    result <- list(eval.points = eval.points, h = h * hmult, 
        h.weights = h.weights, weights = weights)
    Wh <- matrix(rep(h.weights, nnew), ncol = n, byrow = TRUE)
    W1 <- matrix(rep(eval.points[, 1], rep(n, nnew)), ncol = n, 
        byrow = TRUE)
    W1 <- W1 - matrix(rep(x[, 1], nnew), ncol = n, byrow = TRUE)
    W1 <- exp(-0.5 * (W1/(hmult * h[1] * Wh))^2)/Wh
    W2 <- matrix(rep(eval.points[, 2], rep(n, nnew)), ncol = n, 
        byrow = TRUE)
    W2 <- W2 - matrix(rep(x[, 2], nnew), ncol = n, byrow = TRUE)
    W2 <- exp(-0.5 * (W2/(hmult * h[2] * Wh))^2)/Wh
    W3 <- matrix(rep(eval.points[, 3], rep(n, nnew)), ncol = n, 
        byrow = TRUE)
    W3 <- W3 - matrix(rep(x[, 3], nnew), ncol = n, byrow = TRUE)
    W3 <- exp(-0.5 * (W3/(hmult * h[3] * Wh))^2)/Wh
    if (eval.type == "points") {
        est <- as.vector(((W1 * W2 * W3) %*% weights)/(sum(weights) * 
            (2 * pi)^1.5 * h[1] * h[2] * h[3] * hmult^3))
        return(est)
    }
    est <- .vdensecalc(x, h, x, eval.type = "points", weights = weights, 
        filename = filename)
    props <- c(75, 50, 25)
    
    # set colors for different density levels
    cols <- topo.colors(length(props))
    alpha <- seq(1, 0.5, length = length(props))
    levels <- quantile(est, props/100)
    est <- apply(W3, 1, function(x) (W1 %*% (weights * x * t(W2)))/(sum(weights) * 
        (2 * pi)^1.5 * h[1] * h[2] * h[3] * hmult^3))
    est <- array(c(est), dim = rep(ngrid, 3))
    struct <- NULL
    
    # use misc3d to compute the final contours
    if (require(misc3d)) {
        struct <- contour3d(est, levels, eval.points[, 1], eval.points[, 
            2], eval.points[, 3], engine = "none")
    } else {
      stop("In order to draw contour surfaces, please first install the misc3d package!")
    }
    
    # write output to VRML file
    for (i in 1:length(props)) {
        if (length(props) > 1) 
            strct <- struct[[i]]
        else strct <- struct
        trngs.x <- c(t(cbind(strct$v1[, 1], strct$v2[, 1], strct$v3[, 
            1])))
        trngs.y <- c(t(cbind(strct$v1[, 2], strct$v2[, 2], strct$v3[, 
            2])))
        trngs.z <- c(t(cbind(strct$v1[, 3], strct$v2[, 3], strct$v3[, 
            3])))
        a <- list(x = trngs.x, y = trngs.y, z = trngs.z)
        rcol <- (col2rgb(cols[i])/255)[1]
        gcol <- (col2rgb(cols[i])/255)[2]
        bcol <- (col2rgb(cols[i])/255)[3]
        
        # write contours to VRML file
        for (j in seq(1, length(trngs.x), 3)) {
            write(paste("\n Shape { \n\tappearance Appearance {\n\t\t material Material {\n\t\tdiffuseColor ", 
                rcol, " ", gcol, " ", bcol, "\n\t\ttransparency ", 
                alpha[i], " \n\t}\n\t}\t\n \n\t geometry IndexedFaceSet {\n\t \t\tsolid TRUE \t \n\t\tcoord Coordinate {\n\t\t\t point [\n\t\t\t ", 
                trngs.x[j], trngs.y[j], trngs.z[j], " \n\t\t\t  ", 
                trngs.x[j + 1], trngs.y[j + 1], trngs.z[j + 1], 
                " \n\t\t\t  ", trngs.x[j + 2], trngs.y[j + 2], 
                trngs.z[j + 2], " \n\t\t ]\n\t }\n\t coordIndex [\t0 1 2\t] \n\t}\t \n }\n ", 
                sep = " "), file = filename, append = TRUE)
        }
    }
    return(est)
}


# helper function for density estimation contour surface computation
`.hnorm` <-
function (x, weights = NA) 
{
    if (all(is.na(weights))) 
        weights <- rep(1, nrow(x))
    ndim <- ncol(as.matrix(x))
    if (ndim != 3) 
        stop("only data with 3 dimensions is allowed.")
    n <- sum(weights)
    sd <- sqrt(apply(x, 2, function(x, w) {
        sum(w * (x - sum(x * w)/sum(w))^2)/(sum(w) - 1)
    }, w = weights))
    hh <- sd * (4/(5 * n))^(1/7)
    hh
}


# creating 3D-meshes in VRML-format
`.vmesh` <-
function (infile = NULL, x = NULL, y = NULL, z = NULL, edges = NULL, 
    xfun = "sin(v)*cos(u)", yfun = "sin(v)*sin(u)", zfun = "cos(v)", 
    param1 = "u", param2 = "v", range1 = c(0, 2 * pi), range2 = c(0, 
        pi), size1 = 30, size2 = 30, write_obj = FALSE, filename = "out.wrl", 
    cols = "red", scalefac = 4, autoscale = ifelse(is.null(infile), 
        "independent", "equicenter"), lab.axis = c("X-axis", 
        "Y-axis", "Z-axis"), col.axis = "black", showaxis = TRUE, 
    col.lab = "black", col.bg = "white", cex.lab = 1, navigation = "EXAMINE", 
    transparency = 0, fov = 0.785, pos = rep(scalefac + 4, 3), 
    dir = c(0.19, 0.45, 0.87, 2.45), htmlout = NULL, hwidth = 1200, 
    hheight = 800) 
{
    
    # write output in obj-format, if write_obj == TRUE and vertex-
    # and edge-data is available
    if (write_obj) {
    
        if (is.null(x) && is.null(edges)) {
            stop("\nThe generation of obj-files is currently only supported, when using vertex-coordinates (x, y, z-parameters) and polygonal face indices (edges-parameter) as input.")
        }
        
        # write OBJ header        
        write("# OBJ file created by VRMLgen\n#\n", file = filename, 
            append = FALSE)
        
        # extract coordinate vectors from input            
        xyz_parse <- xyz.coords(x, y, z)
        
        # generate the OBJ output file
        write(paste("v", xyz_parse$x, xyz_parse$y, xyz_parse$z, 
            sep = " "), file = filename, append = TRUE)
            
        write("\n", file = filename, append = TRUE)
        
        for (i in 1:nrow(edges)) {
            write(paste("f", paste(edges[i, ], collapse = " "), 
                sep = " "), file = filename, append = TRUE)
        }
        
        cat(paste("\nOutput file \"", filename, "\" was generated in folder ", 
            getwd(), ".\n", sep = ""))
    }
    
    # check which input is supplied:
    else {
    
    		# case 1: the user provides an obj-file as input
        if (!is.null(infile)) {
        
            # read the obj-file
            obj <- strsplit(readLines(infile), "\t")
            
            # collect face data
            faces <- which(as.numeric(lapply(obj, function(x) grep("^f", 
                x[1]))) == 1)
						
						# collect node data                
            nodes <- which(as.numeric(lapply(obj, function(x) grep("^v ", 
                x[1]))) == 1)
                
            # parse node data
            tmpmat <- lapply(obj[nodes], function(x) strsplit(x, 
                " ")[[1]])
                
            # extract edges from faces
            redmat <- t(sapply(tmpmat, function(x) as.numeric(x[2:length(x)])))
            data <- redmat[, which(apply(redmat, 2, function(x) !any(is.na(x))))]
            edges_lst <- lapply(obj[faces], function(x) strsplit(x, 
                " ")[[1]])
            edges_filt1 <- lapply(edges_lst, function(y) as.numeric(sapply(y[2:length(y)], 
                function(x) strsplit(x, "/")[[1]][1])))
            edges_filt <- lapply(edges_filt1, function(x) x[which(!is.na(x))])
            
            # split faces with more than 4 edges into smaller faces
            filter <- sapply(edges_filt, length) > 4
            filt <- edges_filt[filter]
            reduced_lst <- NULL
            if (length(filt) > 0) {
                for (j in 1:length(filt)) {
                  reduced_lst <- c(reduced_lst, list(filt[[j]][1:4]))
                  for (k in seq(4, length(filt[[j]]), 3)) {
                    if (!is.na(filt[[j]][k + 2])) {
                      reduced_lst <- c(reduced_lst, list(c(filt[[j]][k:(k + 
                        2)], filt[[j]][1])))
                    }
                    else if (!is.na(filt[[j]][k + 1])) {
                      reduced_lst <- c(reduced_lst, list(c(filt[[j]][k:(k + 
                        1)], filt[[j]][1], 0)))
                      break
                    }
                    else {
                      break
                    }
                  }
                }
            }
            
            # combine edge data
            comb_lst <- c(edges_filt[!filter], reduced_lst)
            filtless <- which(sapply(comb_lst, length) < 4)
            comb_lst[filtless] <- lapply(comb_lst[filtless], 
                function(x) c(x, 0))
            all(sapply(comb_lst, length) == 4)
            edges <- t(sapply(comb_lst, rbind)) - 1
        }
        
        # case 2: a parametric function is used as input
        else if (is.null(x)) {
        
        		# check if the required data is available
            if (is.null(xfun) || is.null(yfun) || is.null(zfun)) {
                stop("Either the paramater infile or data or the parameters xfun, yfun and zfun have to be specified.")
            }
            if (is.null(param1) || is.null(param2)) {
                stop("The parameter names param1 and param2 have not been specified")
            }
            
            # initialize variables
            x <- NULL
            y <- NULL
            z <- NULL
            
            smin <- range1[1]
            smax <- range1[2]
            tmin <- range2[1]
            tmax <- range2[2]
            
            sn <- size1
            tn <- size2
            ds <- (smax - smin)/sn
            dt <- (tmax - tmin)/tn
            
            # compute the data points representing the parametric function
            for (i in seq(smin, (smax - ds/2), ds)) {
            
                for (j in seq(tmin, (tmax - dt/2), dt)) {
                
                  eval(parse(text = paste(param1, " <- ", i)))
                  eval(parse(text = paste(param2, " <- ", j)))
                  x <- c(x, eval(parse(text = xfun)))
                  y <- c(y, eval(parse(text = yfun)))
                  z <- c(z, eval(parse(text = zfun)))
                  eval(parse(text = paste(param1, " <- ", param1, 
                    " + ds")))
                  x <- c(x, eval(parse(text = xfun)))
                  y <- c(y, eval(parse(text = yfun)))
                  z <- c(z, eval(parse(text = zfun)))
                  eval(parse(text = paste(param2, " <- ", param2, 
                    " + dt")))
                  x <- c(x, eval(parse(text = xfun)))
                  y <- c(y, eval(parse(text = yfun)))
                  z <- c(z, eval(parse(text = zfun)))
                  eval(parse(text = paste(param1, " <- ", param1, 
                    " - ds")))
                  x <- c(x, eval(parse(text = xfun)))
                  y <- c(y, eval(parse(text = yfun)))
                  z <- c(z, eval(parse(text = zfun)))
                }
            }
            
            # combine data
            data <- as.matrix(cbind(x, y, z))
        }
        else {
            
            # mesh vertex data is directly provided as input
            xyz_parse <- xyz.coords(x, y, z)
            data <- cbind(xyz_parse$x, xyz_parse$y, xyz_parse$z)
        }
                
        if (ncol(data) != 3) {
            stop("Data matrix does not have 3 columns!")
        }
        
        
        # verify if .lbar is called within a vrmlgen-environment
    		# initialized with vrml.open()
        scale <- 1
        curdir <- NULL
        if (!exists(".vrmlgenEnv")) {
        
        		# set the general 3D scene settings, if they are
        		# not provided by a higher-level environment
        
        		# write VRML header
            write("#VRML V2.0 utf8\n", file = filename, append = FALSE)
            
            write(paste("\nViewpoint {\n\tfieldOfView", fov, 
                "\n\tposition", pos[1], pos[2], pos[3], "\n\torientation", 
                dir[1], dir[2], dir[3], dir[4], "\n\tjump TRUE\n\tdescription \"viewpoint1\"\n}\n", 
                sep = " "), file = filename, append = TRUE)
                
            write(paste("\nNavigationInfo { type \"", navigation, 
                "\" }\n", sep = ""), file = filename, append = TRUE)
                
            # set background color
            bg_rcol <- (col2rgb(col.bg)/255)[1]
            bg_gcol <- (col2rgb(col.bg)/255)[2]
            bg_bcol <- (col2rgb(col.bg)/255)[3]
            write(paste("Background {\n\t skyColor [\n\t\t ", 
                bg_rcol, bg_gcol, bg_bcol, " \n\t]\n}", sep = " "), 
                file = filename, append = TRUE)
                
        }
        else {
        
        		# obtain the general settings from the
    				# higher level environment
            vrmlgenEnv <- get(".vrmlgenEnv",envir=.GlobalEnv)
            filename <- vrmlgenEnv$filename
            htmlout <- vrmlgenEnv$html
            hheight <- vrmlgenEnv$hheight
            hwidth <- vrmlgenEnv$hwidth
            scale <- vrmlgenEnv$scale
            curdir <- getwd()
            setwd(vrmlgenEnv$VRMLDir)
        }
        
                
        # apply the chosen scaling method
        scaledat <- function(data) {
            diff <- max(data) - min(data)
            if (diff == 0) 
                return(data)
            else return(scale * scalefac * (data - min(data))/(diff))
        }
        
        if (autoscale == "independent") 
            data <- apply(data, 2, scaledat)
        else if (autoscale == "equidist") 
            data <- scale * scalefac * (data/max(data))
        else if (autoscale == "equicenter") 
            data <- scale * scalefac/2 * data/max(data) + scalefac/2
            
        # draw the axes
        if (showaxis) {
        
        		# call the axis3d helper function to draw the axes
            axis3d(lab.axis, filename, type = "vrml", col.lab, 
                col.axis, cex.lab, local_scale = scalefac, global_scale = 1)
        }
                
        # write the data point coordinates
        write(paste("Transform {\nscale 1.0 1.0 1.0\n\tchildren [\n\t\tShape {\n\t\t\tappearance Appearance{ material Material{ transparency ",transparency,"}}\n", 
            sep = " "), file = filename, append = TRUE)
            
        write(paste("\n\t\t\tgeometry IndexedFaceSet {\n\t\t\tcoord DEF\n\tVertexArray Coordinate{\n\tpoint [\n", 
            sep = " "), file = filename, append = TRUE)
            
        for (j in 1:nrow(data)) {
            write(paste(paste(data[j, ], collapse = " "), ",\n", 
                sep = ""), file = filename, append = TRUE)
        }        
        write("]\n}\n", file = filename, append = TRUE)
        
        # write the edges (using coordinate indices)        
        write("coordIndex [\n", file = filename, append = TRUE)
        mat <- NULL
        
        if (!is.null(edges)) {
            mat <- edges
            while (ncol(mat) < 4) {
                mat <- cbind(mat, rep(-1, nrow(mat)))
            }
        }
        else {
            mat <- matrix(1:nrow(data), ncol = 4, byrow = TRUE) - 
                1
        }
        
        for (j in 1:nrow(mat)) {
            write(paste(paste(mat[j, ], collapse = " "), "-1", 
                sep = " "), file = filename, append = TRUE)
        }
        
        
        # write vertex colors
        write(paste("]\n solid FALSE\n colorPerVertex TRUE\n color \nDEF VertexColorArray Color{\ncolor [", 
            sep = ""), file = filename, append = TRUE)
            
        # use rainbow colors if no colors are pre-defined
        if (!length(cols) || (is.null(cols))) {
            
            cols <- rainbow(nrow(data))
            mat <- sapply(cols, function(x) (col2rgb(x)/255))
            for (j in 1:nrow(data)) {
                write(paste(mat[1, j], mat[2, j], mat[3, j], 
                  ",\n", collapse = " "), file = filename, append = TRUE)
            }
        }
                
        # use a single color if the number of colors is less than the number of data rows
        else if ((length(cols) == 1) || (length(cols) < nrow(data))) {
            
            rcol <- (col2rgb(cols[1])/255)[1]
            gcol <- (col2rgb(cols[1])/255)[2]
            bcol <- (col2rgb(cols[1])/255)[3]
            for (j in 1:nrow(data)) {
                write(paste(rcol, gcol, bcol, ",\n", collapse = " "), 
                  file = filename, append = TRUE)
            }
        }
        else {
            
            mat <- sapply(cols, function(x) (col2rgb(x)/255))            
            for (j in 1:nrow(data)) {
                write(paste(mat[1, j], mat[2, j], mat[3, j], 
                  ",\n", collapse = " "), file = filename, append = TRUE)
            }
        }        
        write(paste("]\n}\n}\n}\n]\n}", sep = ""), file = filename, 
            append = TRUE)
        
        
        # create the output files, if this is not already
    		# done by a higher-level lg3d.close-function    
        if (!exists(".vrmlgenEnv")) {
            
            # create HTML output
            if (!is.null(htmlout)) {
                cat("<HTML>", file = htmlout, append = FALSE)
                cat("<HEAD><TITLE>VRMLGen-visualization</TITLE></HEAD><BODY><br>", 
                  file = htmlout, append = FALSE)
                cat(paste("<object type=\"x-world/x-vrml\" data=\"", 
                  filename, "\" ", sep = ""), file = htmlout, 
                  append = TRUE)
                cat(paste("width=\"", hwidth, "\" height=\"", 
                  hheight, "\"><br>", sep = ""), file = htmlout, 
                  append = TRUE)
                cat(paste("<param name=\"src\" value=\"", filename, 
                  "\"><br>", sep = ""), file = htmlout, append = TRUE)
                cat("Your browser cannot display VRML files.<br>Please INSTALL a VRML-plugin or open the file in an external VRML-viewer.</object><br>", 
                  file = htmlout, append = TRUE)
                cat("</BODY></HTML>", file = htmlout, append = TRUE)
            }
            
            # show success message        
            cat(paste("\nOutput file \"", filename, "\" was generated in folder ", 
                getwd(), ".\n", sep = ""))
        }
        else {        
        		# return to current directory, if higher-level plotting
        		# functions wrote re-directed the output to a temporary directory
            setwd(curdir)
        }
    }
}


# lower-level function for plotting data points in the VRML-format
`.vrmlpoints` <-
function (x, y, z, filename, rcol, gcol, bcol, pointstyle, hyperlinks = NULL, 
    global_scale = 1, local_scale = 1, transparency = 0) 
{
    
    # add hyperlinks to the data points
    hyperlinks_start <- rep("", length(x))
    hyperlinks_end <- rep("", length(x))
    
    if (!is.null(hyperlinks)) {
    
        hyperlinks_start <- paste("\n\tchildren Anchor { url \"", 
            hyperlinks, "\"", sep = "")
        hyperlinks_end <- rep("}", length(x))
    }
    
    # check the point-style and draw the data points at the
    # corresponding coordinates
    if (pointstyle == "s") {
    
        for (i in 1:length(x))
        	  write(paste("\nTransform {\n\ttranslation ", 
            x[i] * global_scale, y[i] * global_scale, z[i] * 
                global_scale, hyperlinks_start[i], "\n\tchildren Shape {\n\t\tappearance Appearance { material Material { diffuseColor", 
            rcol, gcol, bcol, "transparency", transparency, "} }\n\tgeometry Sphere { radius ", 
            0.08 * local_scale, "}\n\t\t}", hyperlinks_end[i], 
            "\n}", sep = " "), file = filename, append = TRUE)
    }
    else if (pointstyle == "c") {
    
        for (i in 1:length(x))
        	  write(paste("\nTransform {\n\trotation 1 0 0 1.5708\n\ttranslation ", 
            x[i] * global_scale, y[i] * global_scale, z[i] * 
                global_scale, hyperlinks_start[i], "\n\tchildren Shape {\n\t\tappearance Appearance { material Material { diffuseColor", 
            rcol, gcol, bcol, "transparency", transparency, "} }\n\tgeometry Cone { bottomRadius ", 
            0.08 * local_scale, "height", 0.08 * local_scale, 
            "\n\t\tside TRUE}\n\t\t}", hyperlinks_end[i], "\n}", 
            sep = " "), file = filename, append = TRUE)
    }
    else {
        for (i in 1:length(x))
        		write(paste("\nTransform {\n\ttranslation ", 
            x[i] * global_scale, y[i] * global_scale, z[i] * 
                global_scale, hyperlinks_start[i], "\n\tchildren Shape {\n\t\tappearance Appearance { material Material { diffuseColor", 
            rcol, gcol, bcol, "transparency", transparency, "} }\n\tgeometry Box { size ", 
            0.08 * local_scale, 0.08 * local_scale, 0.08 * local_scale, 
            "}\n\t\t}", hyperlinks_end[i], "\n}", sep = " "), 
            file = filename, append = TRUE)
    }
}


# lower-level function for adding text to a 3D-scene
# in the Livegraphics3D-format
`.vrmltext` <-
function (x, y, z, text, filename, col, scale, rot, fontfamily, 
    fontweight, hyperlink) 
{
    
    # convert the color to RGB-format
    rcol <- (col2rgb(col)/255)[1]
    gcol <- (col2rgb(col)/255)[2]
    bcol <- (col2rgb(col)/255)[3]
    
    # add hyperlinks to the text
    if (is.null(hyperlink)) {
        
        for (i in 1:length(x))
        		write(paste("Transform {\n\ttranslation ", 
            x[i], y[i], z[i], "\n\tscale ", 0.28, 0.28, 0.28, 
            "\n\trotation ", rot[1], rot[2], rot[3], rot[4], 
            "\n\tchildren Shape {\n\t\tappearance Appearance { material Material {diffuseColor ", 
            rcol, gcol, bcol, "  } }\n\tgeometry Text { string \"", 
            text[i], "\"\nfontStyle FontStyle {\nsize ", scale, 
            "\nfamily [\"", fontfamily, "\"]\njustify \"MIDDLE\"\n style \"", 
            fontweight, "\" }\n }\n\t}\n}", sep = " "), file = filename, 
            append = TRUE)
    }
    else {
        
        for (i in 1:length(x))
        		write(paste("Transform {\n\ttranslation ", 
            x[i], y[i], z[i], "\n\tscale ", 0.28, 0.28, 0.28, 
            "\n\trotation ", rot[1], rot[2], rot[3], rot[4], 
            "\n\tchildren Anchor { url \"", hyperlink, "\" \n\tchildren Shape {\n\t\tappearance Appearance { material Material {diffuseColor ", 
            rcol, gcol, bcol, "  } }\n\tgeometry Text { string \"", 
            text[i], "\"\nfontStyle FontStyle {\nsize ", scale, 
            "\nfamily [\"SANS\"]\njustify \"MIDDLE\"\n style \"", 
            fontweight, "\" }\n }\n\t}\n}\n}", sep = " "), file = filename, 
            append = TRUE)
    }
}


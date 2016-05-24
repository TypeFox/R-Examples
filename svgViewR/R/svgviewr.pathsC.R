svgviewr.pathsC <- function(path, file=NULL, col=NULL, col.fill="none", col.stroke = "black", z.index = 0, 
	layer = "", label = "", lwd = 1, opacity.stroke = 1, opacity.fill = 1, index.add = 0, append = TRUE){

	# IF PATH IS VECTOR CONVERT TO LIST
	if(!is.list(path) && !is.matrix(path)) path <- list(path)

	# IF PATH IS MATRIX, CONVERT TO LIST
	if(is.matrix(path)){
		pathn <- list()
		for(i in 1:nrow(path)) pathn[[i]] <- path[i, ]
		path <- pathn
	}

	if(length(path) == 0) return(1)

	# IF COL IS SPECIFIED, OVERWRITE FILL AND STROKE
	if(!is.null(col)){col.fill <- col;col.stroke <- col}

	# SET GRAPHICAL PARAMETERS
	svg_gp <- c("col", "col.fill", "col.stroke", "label", "layer", "opacity.fill", "opacity.stroke", "lwd", "z.index")

	# CONVERT GRAPHICAL PARAMETERS TO VECTORS WITH SAME NUMBER OF ELEMENTS OF FIRST X DIMENSION
	for(gpar in svg_gp) if(length(get(gpar)) == 1) assign(gpar, rep(get(gpar), length(path)))

	# WRITE PATHS TO SVG STRUCTURE
	new_lines <- rep(NA, length(path))
	for(i in 1:length(path)){

		if(is.null(path[[i]])) next

		# IF INITIAL INDEX IS NON-ZERO ADD INDEX_I TO ALL PATH_LIST VALUES
		if(index.add > 0) path[[i]] <- path[[i]] + rep(index.add, length(path[[i]]))		

		# SAVE ROWS AS SPACE DELIMITED VECTOR FOR SVG, SUBTRACT 1 SINCE JAVASCRIPT VECTORS START AT 0
		new_lines[i] <- paste("\t<pathC z-index=\"", z.index[i], "\" layer=\"", layer[i], "\" label=\"", label[i], "\" stroke=\"", col.stroke[i], "\" stroke-width=\"", lwd[i], "\" stroke-opacity=\"", opacity.stroke[i], "\" fill=\"", col.fill[i], "\" fill-opacity=\"", opacity.fill[i], "\" d=\"", paste(path[[i]], collapse=","), "\" />", sep="")
	}

	# IF FILE IS NULL, RETURN LINES OF SVG OBJECTS
	if(is.null(file)) return(new_lines)

	# SAVE NEW LINES TO FILE
	svgviewr.write(new_lines, file, append=append)
}
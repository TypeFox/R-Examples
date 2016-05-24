drawShapes <- function(shapes, file, path.connect = NULL, connect.curve.points = TRUE, 
	window.title = NULL, animate = TRUE, animate.duration = 1, animate.reverse = FALSE, 
	animate.repeat = -1, draw='all', path.fixed = TRUE, 
	point.cex=1, point.col.fill="black", 
	point.col.stroke="black", point.lwd=2, path.col.fill='none', path.opacity.fill=1, 
	path.opacity.stroke=1, path.col.stroke="black", path.lwd = 1, 
	add = FALSE, fdir = NULL){

	#require(svgViewR)

	if(!is.list(shapes)){

		# Read shapes file(s)
		shapes <- readShapes(shapes)

		# Make sure file ends with .html
		file <- gsub('[.]txt$', '', file, ignore.case=TRUE)

		# Make sure file ends with .html
		if(!grepl('[.]html$', file, ignore.case=TRUE)) file <- paste0(file, '.html')
	}

	# Get proper shapes level
	if('shapes' %in% names(shapes)) shapes <- shapes$shapes

	points <- NULL
	curves <- NULL
	
	if(draw == 'all') draw <- c('landmarks', 'curves')

	# Get shapes
	if('landmarks.pixel' %in% names(shapes) && 'landmarks' %in% draw) points <- shapes$landmarks.pixel
	if('landmarks.scaled' %in% names(shapes) && 'landmarks' %in% draw) points <- shapes$landmarks.scaled
	if('landmarks' %in% names(shapes) && 'landmarks' %in% draw) points <- shapes$landmarks
	if('curves.pixel' %in% names(shapes) && 'curves' %in% draw) curves <- shapes$curves.pixel
	if('curves.scaled' %in% names(shapes) && 'curves' %in% draw) curves <- shapes$curves.scaled
	if('curves' %in% names(shapes) && 'curves' %in% draw) curves <- shapes$curves
	
	# Read path.connect
	if(!is.list(path.connect)){
	
		# Check that file exists
		if(!file.exists(path.connect)) stop(paste0("File '", path.connect, "' not found."))
		
		# Read lines
		read_lines <- suppressWarnings(readLines(con=path.connect))
		
		path_connect <- list()
		for(i in 1:length(read_lines)){
			if(read_lines[i] == '') next
			path_connect[[length(path_connect)+1]] <- gsub('(^[ ]+)|([ ]+$)', '', strsplit(read_lines[i], ',')[[1]])
		}

	}else{
		path_connect <- path.connect
	}

	# SET PATHS CONNECTING POINTS
	path_list <- NULL
	if(!is.null(points) && !is.null(path_connect)){

		# CREATE PATH LIST TO CONNECT POINTS
		path_list <- vector("list", length(path_connect))

		if(!is.numeric(path_connect[[1]]) && !is.null(dimnames(points)[[1]])){
			for(i in 1:length(path_connect)){
				for(j in 1:length(path_connect[[i]])){

					# FIND REVERSE MARKER
					rev_str <- FALSE
					if(grepl('^rev:', path_connect[[i]][j])) rev_str <- TRUE
					
					# REMOVE MARKER
					if(rev_str) path_connect[[i]][j] <- gsub('^rev:', '', path_connect[[i]][j])

					# PATH LABELS TO ROWNAMES OF POINT ARRAY
					if(path.fixed){
						grepl_match <- grepl(paste0('^', path_connect[[i]][j], '$'), dimnames(points)[[1]])
					}else{
						grepl_match <- grepl(path_connect[[i]][j], dimnames(points)[[1]])
					}

					# ADD TO PATH LIST
					if(sum(grepl_match) > 0){

						match <- which(grepl_match)

						# REVERSE MATCH
						if(rev_str) match <- rev(match)

						path_list[[i]] <- c(path_list[[i]], match)
					}
				}
			}
		}else{
			for(i in 1:length(path_connect)) path_list[[i]] <- path_connect[[i]]
		}
	}

	# Get window title
	if(is.null(window.title)){window_title <- 'Shapes Viewer'}else{window_title <- window.title}

	# Set whether to show control panel
	show_control <- FALSE
	start_rotate <- TRUE
	if(animate && length(dim(points)) > 2 && dim(points)[3] > 1) show_control <- TRUE
	if(dim(points)[2] == 2) start_rotate <- FALSE
	#if(length(dim(points)) == 2 || (length(dim(points)) > 2 && dim(points)[3] == 1)) start_rotate <- FALSE

	# Create new html file
	if(!add){
		svgviewr.new(file=file, window.title=window_title, animate.duration=animate.duration, 
			animate.reverse=animate.reverse, animate.repeat=animate.repeat, 
			fdir=fdir)
		#show.control=show_control, start.rotate=start_rotate, 
	}

	# Draw landmarks
	if(!is.null(points)){
		svgviewr.points(points, file=file, col.fill=point.col.fill, 
			col.stroke=point.col.stroke, cex=point.cex, lwd=point.lwd)
	}

	# SET CURVE CONNECT START
	ccstart <- 0

	# DRAW PATHS CONNECTING POINTS
	if(!is.null(path_list)){

		# CONNECT POINTS WITH PATHS
		if(animate){
			for(j in 1:length(path_list)){
				if(is.null(path_list[[j]])) next
				path_list[[j]] <- path_list[[j]] + ccstart
			}
			svgviewr.pathsC(path_list, file=file, col.fill=path.col.fill, opacity.fill=path.opacity.fill, col.stroke=path.col.stroke, opacity.stroke=path.opacity.stroke, lwd=path.lwd)
		}else{
			for(i in 1:dim(points)[3]){
				for(j in 1:length(path_list)){
					if(is.null(path_list[[j]])) next
					path_list[[j]] <- path_list[[j]] + ccstart + dim(points)[1]*(i-1)
				}
				svgviewr.pathsC(path_list, file=file, col.fill=path.col.fill, opacity.fill=path.opacity.fill, col.stroke=path.col.stroke, opacity.stroke=path.opacity.stroke, lwd=path.lwd)
			}
		}
	}
	
	ccstart <- nrow(points)

	# CREATE LIST FOR CONNECTING CURVE POINTS
	connect_curves <- list()

	# DRAW CURVE POINTS (ONLY WITH SINGLE FILE) FALSE && 
	if(!is.null(curves) && length(curves) > 0 && is.matrix(curves[[1]])){
		
		curve_points <- matrix(NA, nrow=0, ncol=ncol(curves[[1]]))
		
		for(i in 1:length(curves)){

			# GET CURVE
			curve <- curves[[names(curves)[i]]]

			# ADD POINTS TO CURVE POINT MATRIX
			curve_points <- rbind(curve_points, curve)

			# ADD VECTOR OF POINT INDICES TO CONNECT PATH LIST
			connect_curves[[i]] <- (ccstart + 1):(ccstart + nrow(curve))
			
			# UPDATE PATHC START INDEX
			ccstart <- ccstart + nrow(curve)
		}

		# ADD CURVE POINTS TO HTML FILE
		svgviewr.points(curve_points, file=file, col.fill=point.col.fill, 
			col.stroke=point.col.stroke, cex=point.cex, lwd=point.lwd)

		# ADD PATHSC FOR CURVES
		if(connect.curve.points) svgviewr.pathsC(connect_curves, file=file, 
			col.fill=path.col.fill, opacity.fill=path.opacity.fill, col.stroke=path.col.stroke, 
			opacity.stroke=path.opacity.stroke, lwd=path.lwd)

		names(connect_curves) <- names(curves)
	}
	
	return(1)
}
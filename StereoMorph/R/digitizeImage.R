digitizeImage <- function(image.file, landmarks.file=NULL, control.points.file=NULL, curve.points.file=NULL, shapes.file=NULL, landmarks.ref=NULL, 
	curves.ref=NULL, image.id=NULL, landmark.color.blur = 'blue', landmark.color.focus = 'green', curve.color.blur = 'purple', 
	control.point.color.blur = 'purple', control.point.color.focus = 'red', landmark.radius = 4, control.point.radius = 4, marker.stroke.width = 2,
	app.dir=NULL){

	if(!is.null(shapes.file) && !is.null(landmarks.file)) 
		cat("Warning: If 'shapes.file' is non-NULL, 'landmarks.file' must be NULL.\n")

	if(!is.null(shapes.file) && !is.null(control.points.file)) 
		cat("Warning: If 'shapes.file' is non-NULL, 'control.points.file' must be NULL.\n")

	if(is.null(shapes.file) && is.null(landmarks.file) && is.null(control.points.file)) 
		cat("Warning: 'shapes.file', 'landmarks.file' and 'control.points.file' are all NULL. If you would like to save landmarks or curves, please specify a location to save either the landmarks and/or curves.\n")

	# GET STEREOMORPH SHINY APP DIRECTORY
	if(is.null(app.dir)){
		app_dir <- paste0(path.package("StereoMorph"), '/extdata/apps/digitizeImage')
	}else{
		app_dir <- app.dir
	}

	# REMOVE ANY IMAGE FILES IN IMG FOLDER
	if(length(list.files(paste0(app_dir, '/www/img/'))) > 0)
		file.remove(paste0(app_dir, '/www/img/', list.files(paste0(app_dir, '/www/img/'))))

	# IF image.file IS DIRECTORY, ADD ALL IMAGE FILES FROM DIRECTORY
	if(length(list.files(image.file) > 0)){
		is_image <- grepl(pattern='[.]jpg$|[.]jpeg$|[.]tif$|[.]tiff$|[.]png$', x=list.files(image.file), ignore.case=TRUE)
		image.file <- paste0(image.file, '/', list.files(image.file)[is_image])
	}
	
	# CHECK THAT FILES EXIST
	if(sum(!file.exists(image.file)) > 0){
		stop(paste0("The following files were not found:\n\t", paste(image.file[!file.exists(image.file)], collapse="\n\t")))
	}

	# CHECK IMAGE EXTENSIONS FOR COMPATIBILITY
	check_img_type <- grepl(pattern='[.]jpg$|[.]jpeg$|[.]tif$|[.]tiff$|[.]png$', x=image.file, ignore.case=TRUE)
	if(sum(!check_img_type) > 0) stop(paste0("Only images of type JPG, JPEG, TIF, TIFF and PNG are currently supported. The following are unsupported file types:\n\t", paste(image.file[!check_img_type], collapse="\n\t")))

	# GET IMAGE NAMES
	img_file_split <- strsplit(image.file, '/')
	img_names <- rep(NA, length(image.file))
	for(i in 1:length(img_file_split)) img_names[i] <- img_file_split[[i]][length(img_file_split[[i]])]

	# SET IMAGE IDS
	if(is.null(image.id)) image.id <- gsub(" ", "_", gsub('[.][a-zA-Z]+$', '', img_names))

	# IF ANY OF THE OUTPUT PATHS ARE DIRECTORIES, MAKE NAMES SAME AS IMAGE NAMES
	if(sum(grepl('[.]txt$', shapes.file)) != length(shapes.file)) shapes.file <- paste0(shapes.file, '/', gsub('[.][a-zA-Z]+$', '.txt', img_names))
	if(sum(grepl('[.]txt$', landmarks.file)) != length(landmarks.file)) landmarks.file <- paste0(landmarks.file, '/', gsub('[.][a-zA-Z]+$', '.txt', img_names))
	if(sum(grepl('[.]txt$', control.points.file)) != length(control.points.file)) control.points.file <- paste0(control.points.file, '/', gsub('[.][a-zA-Z]+$', '.txt', img_names))
	if(sum(grepl('[.]txt$', curve.points.file)) != length(curve.points.file)) curve.points.file <- paste0(curve.points.file, '/', gsub('[.][a-zA-Z]+$', '.txt', img_names))

#	if(!is.null(landmarks.file) && length(list.files(landmarks.file) > 0)) landmarks.file <- paste0(landmarks.file, '/', gsub('.[a-zA-Z]+$', '.txt', img_names))
#	if(!is.null(control.points.file) && length(list.files(control.points.file) > 0)) control.points.file <- paste0(control.points.file, '/', gsub('.[a-zA-Z]+$', '.txt', img_names))
#	if(!is.null(curve.points.file) && length(list.files(curve.points.file) > 0)) curve.points.file <- paste0(curve.points.file, '/', gsub('.[a-zA-Z]+$', '.txt', img_names))

	# CHECK THAT NUMBER OF IMAGES MATCHES OUTPUT PATHS IF NOT NULL
	if(!is.null(shapes.file) && length(image.file) != length(shapes.file))
		stop(paste0("The number of image file paths input (", length(image.file), ") does not match the number of shape file paths input (", length(shapes.file), ")."))
	if(!is.null(landmarks.file) && length(image.file) != length(landmarks.file))
		stop(paste0("The number of image file paths input (", length(image.file), ") does not match the number of landmark file paths input (", length(landmarks.file), ")."))
	if(!is.null(control.points.file) && length(image.file) != length(control.points.file))
		stop(paste0("The number of image file paths input (", length(image.file), ") does not match the number of control point file paths input (", length(control.points.file), ")."))
	if(!is.null(curve.points.file) && length(image.file) != length(curve.points.file))
		stop(paste0("The number of image file paths input (", length(image.file), ") does not match the number of curve point file paths input (", length(curve.points.file), ")."))

	# READ IN LANDMARKS
	landmarks_ref <- NULL
	if(!is.null(landmarks.ref)){
		
		# SINGLE ELEMENT IN VECTOR AND EXISTS AS A FILE
		if(length(landmarks.ref) == 1) if(file.exists(landmarks.ref)) landmarks.ref <- as.vector(read.table(landmarks.ref, sep="\n")[,1])

		landmarks_ref <- landmarks.ref
	}

	# READ IN CURVE REF
	curves_ref <- list()
	landmarks_from_curves_ref <- rep(NA, 0)
	if(!is.null(curves.ref)){

		# NOT A MATRIX - ASSUME TO BE FILE PATH
		if(!is.matrix(curves.ref)) curves.ref <- suppressWarnings(as.matrix(read.table(curves.ref, sep="\t")))

		colnames(curves.ref) <- NULL
		for(i in 1:nrow(curves.ref)){
			curves_ref[[i]] <- c(curves.ref[i, ])
			landmarks_from_curves_ref <- c(landmarks_from_curves_ref, curves.ref[i, 2:3])
		}
	}

	# ADD LANDMARKS FROM CURVES
	if(!is.null(landmarks.file) || !is.null(shapes.file)){
		landmarks_from_curves_ref <- unique(landmarks_from_curves_ref)
		landmarks_from_curves_ref <- landmarks_from_curves_ref[!landmarks_from_curves_ref %in% landmarks_ref]
		landmarks_ref <- c(landmarks_ref, landmarks_from_curves_ref)
	}

	# SET INITIAL SETTINGS
	run_app <- list('settings'=list(
		'copy_landmarks'=FALSE,
		'copy_curves'=FALSE,
		'copy_ruler_points'=FALSE,
		'copy_corners'=FALSE,
		'copy_scaling'=FALSE
	))

	img_num <- prev_img_num <- 1
	while(img_num <= length(image.file)){

		# SET INITIAL PARAMETERS
		init_params <- list()
		init_params$app_dir <- app_dir
		init_params$prev_wd <- getwd()
		init_params$img_name <- gsub(" ", "_", img_names[img_num])
		init_params$image_id <- image.id[img_num]
		init_params$img_size <- file.info(image.file[img_num])$size
		init_params$img_file <- image.file[img_num]
		init_params$scaling_units <- 'NA'
		init_params$scaling <- 'NA'
		init_params$ruler_interval <- 'NA'
		init_params$checkerboard_nx <- 'NA'
		init_params$checkerboard_ny <- 'NA'
		
		if(!is.null(shapes.file)) 
			init_params$shapes_file <- shapes.file[img_num]
		if(!is.null(landmarks.file)) 
			init_params$landmarks_file <- landmarks.file[img_num]
		if(!is.null(control.points.file)) 
			init_params$control_points_file <- control.points.file[img_num]
		if(!is.null(curve.points.file)) 
			init_params$curve_points_file <- curve.points.file[img_num]

		# GET IMAGE DIMENSIONS
		if(grepl(pattern='[.]jpg$|[.]jpeg$', x=image.file[img_num], ignore.case=TRUE)) img_dim <- dim(readJPEG(image.file[img_num], native=TRUE))
		if(grepl(pattern='[.]png$', x=image.file[img_num], ignore.case=TRUE)) img_dim <- dim(readPNG(image.file[img_num], native=TRUE))
		if(grepl(pattern='[.]tif$|[.]tiff$', x=image.file[img_num], ignore.case=TRUE)) img_dim <- dim(readTIFF(image.file[img_num], native=TRUE))
		
		init_params$img_width <- img_dim[2]
		init_params$img_height <- img_dim[1]

		init_params$landmark_color_blur <- landmark.color.blur
		init_params$landmark_color_focus <- landmark.color.focus
		init_params$curve_color_blur <- curve.color.blur
		init_params$control_point_color_blur <- control.point.color.blur
		init_params$control_point_color_focus <- control.point.color.focus
		init_params$landmark_radius <- landmark.radius
		init_params$control_point_radius <- control.point.radius
		init_params$marker_stroke_width <- marker.stroke.width
		init_params$landmarks_ref <- landmarks_ref
		init_params$curves_ref <- curves_ref
		init_params$unsaved_landmarks <- FALSE
		init_params$unsaved_curves <- FALSE

		# SET WHETHER PREVIOUS OR CURRENT IMAGES EXIST
		init_params$prev_img <- FALSE
		init_params$next_img <- FALSE
		if(img_num < length(image.file)) init_params$next_img <- TRUE
		if(img_num > 1) init_params$prev_img <- TRUE

		# COPY IMAGE TO WWW FOLDER
		file.copy(image.file[img_num], paste0(app_dir, '/www/img/', gsub(" ", "_", img_names[img_num])))

		# READ IN SHAPES
		init_params$landmarks <- list()
		init_params$control_points <- list()
		rulers_img_num <- img_num
		landmark_img_num <- img_num
		curve_img_num <- img_num
		scaling_img_num <- img_num
		corners_img_num <- img_num

		if(!is.null(shapes.file)){

			# READ IN INITIAL FILE
			if(file.exists(shapes.file[img_num])) read_shapes_init <- readShapes(shapes.file[rulers_img_num], 
				c('ruler.pixel', 'scaling.units', 'ruler.interval', 'checkerboard.nx', 'checkerboard.ny', 'square.pixel', 'square.size', 
				'ruler.points', 'checker.pixel', 'landmarks.pixel', 'curves.control'))

			# IF COPYING OVER RULER POINTS FROM PREVIOUS FILE
			if(run_app$settings$copy_scaling && file.exists(shapes.file[prev_img_num])){
				read_shapes <- readShapes(shapes.file[prev_img_num])
				if(!is.null(read_shapes$ruler.pixel)) init_params$ruler_pixel <- read_shapes$ruler.pixel
				if(!is.null(read_shapes$scaling.units)) init_params$scaling_units <- read_shapes$scaling.units
				if(!is.null(read_shapes$ruler.interval)) init_params$ruler_interval <- read_shapes$ruler.interval
				if(!is.null(read_shapes$checkerboard.nx)) init_params$checkerboard_nx <- read_shapes$checkerboard.nx
				if(!is.null(read_shapes$checkerboard.ny)) init_params$checkerboard_ny <- read_shapes$checkerboard.ny
				if(!is.null(read_shapes$square.pixel)) init_params$checker_square_pixel <- read_shapes$square.pixel
				if(!is.null(read_shapes$square.size)) init_params$checker_square_world <- read_shapes$square.size
			}

			# GET SCALING DATA IF AVAILABLE
			if(file.exists(shapes.file[scaling_img_num])){
				if(!is.null(read_shapes_init$ruler.pixel)) init_params$ruler_pixel <- read_shapes_init$ruler.pixel
				if(!is.null(read_shapes_init$scaling.units)) init_params$scaling_units <- read_shapes_init$scaling.units
				if(!is.null(read_shapes_init$ruler.interval)) init_params$ruler_interval <- read_shapes_init$ruler.interval
				if(!is.null(read_shapes_init$checkerboard.nx)) init_params$checkerboard_nx <- read_shapes_init$checkerboard.nx
				if(!is.null(read_shapes_init$checkerboard.ny)) init_params$checkerboard_ny <- read_shapes_init$checkerboard.ny
				if(!is.null(read_shapes_init$square.pixel)) init_params$checker_square_pixel <- read_shapes_init$square.pixel
				if(!is.null(read_shapes_init$square.size)) init_params$checker_square_world <- read_shapes_init$square.size
			}
			
			# TRY READING IN RULER POINTS FROM CURRENT FILE
			ruler.points <- NULL
			if(file.exists(shapes.file[rulers_img_num])) ruler.points <- read_shapes_init$ruler.points

			# IF NO RULER POINTS AVAILABLE FROM CURRENT FILE
			if(!file.exists(shapes.file[rulers_img_num]) || (is.null(ruler.points) || nrow(ruler.points) == 0)){

				# IF COPYING OVER RULER POINTS FROM PREVIOUS FILE
				if(run_app$settings$copy_ruler_points && file.exists(shapes.file[prev_img_num])){
					rulers_img_num <- prev_img_num
					ruler.points <- readShapes(shapes.file[prev_img_num], 'ruler.points')$ruler.points
				}
			}

			# IF ANY RULER POINTS WERE FOUND
			if(!is.null(ruler.points) && nrow(ruler.points) > 0){
				for(i in 1:nrow(ruler.points)) init_params$ruler_points[[i]] <- c(rownames(ruler.points)[i], ruler.points[i, ])
				if(img_num != rulers_img_num) init_params$unsaved_ruler_points<- TRUE
			}


			# TRY READING IN CORNERS FROM CURRENT FILE
			corners <- NULL
			if(file.exists(shapes.file[corners_img_num])){
				corners <- read_shapes_init$checker.pixel
				if(!is.null(read_shapes_init$checkerboard.nx)) init_params$checkerboard_nx <- read_shapes_init$checkerboard.nx
				if(!is.null(read_shapes_init$checkerboard.ny)) init_params$checkerboard_ny <- read_shapes_init$checkerboard.ny
			}

			# IF NO CORNERS AVAILABLE FROM CURRENT FILE
			if(!file.exists(shapes.file[corners_img_num]) || (is.null(corners) || nrow(corners) == 0)){

				# IF COPYING OVER CORNERS FROM PREVIOUS FILE
				if(run_app$settings$copy_corners && file.exists(shapes.file[prev_img_num])){
					corners_img_num <- prev_img_num
					read_shapes <- readShapes(shapes.file[prev_img_num])
					corners <- read_shapes$checker.pixel
					if(!is.null(read_shapes$checkerboard.nx)) init_params$checkerboard_nx <- read_shapes$checkerboard.nx
					if(!is.null(read_shapes$checkerboard.ny)) init_params$checkerboard_ny <- read_shapes$checkerboard.ny
				}
			}

			# IF ANY CORNERS WERE FOUND
			if(!is.null(corners) && nrow(corners) > 0){
				for(i in 1:nrow(corners)) init_params$corners[[i]] <- c(rownames(corners)[i], corners[i, ])
				if(img_num != corners_img_num) init_params$unsaved_corners <- TRUE
			}


			# TRY READING IN LANDMARKS FROM CURRENT FILE
			landmarks.pixel <- NULL
			if(file.exists(shapes.file[landmark_img_num])) landmarks.pixel <- read_shapes_init$landmarks.pixel

			# IF NO LANDMARKS AVAILABLE FROM CURRENT FILE
			if(!file.exists(shapes.file[landmark_img_num]) || (is.null(landmarks.pixel) || nrow(landmarks.pixel) == 0)){

				# IF COPYING OVER LANDMARKS FROM PREVIOUS FILE
				if(run_app$settings$copy_landmarks && file.exists(shapes.file[prev_img_num])){
					landmark_img_num <- prev_img_num
					landmarks.pixel <- readShapes(shapes.file[prev_img_num], 'landmarks.pixel')$landmarks.pixel
				}
			}

			# IF ANY LANDMARKS WERE FOUND
			if(!is.null(landmarks.pixel) && nrow(landmarks.pixel) > 0){
				for(i in 1:nrow(landmarks.pixel)) init_params$landmarks[[i]] <- c(rownames(landmarks.pixel)[i], landmarks.pixel[i, ])
				if(img_num != landmark_img_num) init_params$unsaved_landmarks <- TRUE
			}


			# TRY READING IN CURVE CONTROL POINTS FROM CURRENT FILE
			curves.control <- NULL
			if(file.exists(shapes.file[curve_img_num])) curves.control <- read_shapes_init$curves.control

			# IF NO CURVES AVAILABLE FROM CURRENT FILE
			if(!file.exists(shapes.file[curve_img_num]) || (is.null(curves.control) || length(curves.control) == 0)){

				# IF COPYING OVER CURVES FROM PREVIOUS FILE
				if(run_app$settings$copy_curves && file.exists(shapes.file[prev_img_num])){
					curve_img_num <- prev_img_num
					curves.control <- readShapes(shapes.file[prev_img_num], 'curves.control')$curves.control
				}
			}

			# IF ANY CURVES WERE FOUND
			if(!is.null(curves.control) && length(curves.control) > 0){
				for(i in 1:length(curves.control)){
					init_params$control_points[[i]] <- c(names(curves.control)[i], t(curves.control[[names(curves.control)[i]]]))
				}
				if(img_num != curve_img_num) init_params$unsaved_curves <- TRUE
			}
			
		}else{

			# READ IN CURRENT LANDMARKS
			if(!is.null(landmarks.file) && file.exists(landmarks.file[img_num]) && file.info(landmarks.file[img_num])$size > 1){}else{
				if(run_app$settings$copy_landmarks) landmark_img_num <- prev_img_num
			}
			if(!is.null(landmarks.file) && file.exists(landmarks.file[landmark_img_num]) && file.info(landmarks.file[landmark_img_num])$size > 1){
				landmarks <- as.matrix(read.table(landmarks.file[landmark_img_num], row.names=1, sep="\t"))
				colnames(landmarks) <- NULL
				for(i in 1:nrow(landmarks)) init_params$landmarks[[i]] <- c(rownames(landmarks)[i], landmarks[i, ])
				if(img_num != landmark_img_num) init_params$unsaved_landmarks <- TRUE
			}

			# READ IN CURRENT CONTROL POINTS	
			if(!is.null(control.points.file) && file.exists(control.points.file[img_num]) && file.info(control.points.file[img_num])$size > 1){}else{
				if(run_app$settings$copy_curves) curve_img_num <- prev_img_num
			}
			if(!is.null(control.points.file) && file.exists(control.points.file[curve_img_num]) && file.info(control.points.file[curve_img_num])$size > 1){
				control_points <- readBezierControlPoints(control.points.file[curve_img_num])
				for(i in 1:length(control_points)){
					init_params$control_points[[i]] <- c(names(control_points)[i], c(t(control_points[[names(control_points)[i]]][[1]])))
				}
				if(img_num != curve_img_num) init_params$unsaved_curves <- TRUE
			}
		}
		
		# CONVERT PARAMETERS INTO STRING FOR R TO READ
		param_str <- 'init_params <- list()\n'
		for(param_name in names(init_params)){

			if(!is.list(init_params[[param_name]])){
				param_str <- paste0(param_str, paste0('init_params[[\'', param_name, '\']]',  ' <- ', 'c(\'', 
					paste0(init_params[[param_name]], collapse='\',\''), '\')', sep='\n'))
			}else{
				
				if(length(init_params[[param_name]]) == 0) next
				for(i in 1:length(init_params[[param_name]])){
					param_str <- paste0(param_str, paste0('init_params[[\'', param_name, '\']][[', i, ']]',  ' <- ', 'c(\'', 
						paste0(init_params[[param_name]][[i]], collapse='\',\''), '\')', sep='\n'))
				}
			}

			param_str <- paste0(param_str, '\n')
		}

		json_str <- paste0("json_string <- '", listToJSONStr(init_params), "'")

		# WRITE PARAMETERS TO FILE IN ORDER TO PASS THEM TO SHINY APP
		write(x = paste0(param_str, '\n', json_str), file = paste0(app_dir, "/initial_parameters.R"))

		# COPY HTML FILE WITHOUT IMAGE TAG
		file.copy(paste0(app_dir, "/digitize_image_pre.html"), paste0(app_dir, "/digitize_image.html"), overwrite=TRUE)

		# ADD IMAGE TAG TO HTML DOCUMENT (TO LOAD IMAGE AND GET SIZE)
		img_tag <- paste0('\n<img style="display:none;" id="img1" src="img/', gsub(" ", "_", img_names[img_num]),'" ></img>')
		write(img_tag, file=paste0(app_dir, "/digitize_image.html"), append=TRUE)

		# NOTIFICATION IN R CONSOLE
		#cat(paste0("Loading image '", img_names[img_num], "'"))

		# SAVE CURRENT IMAGE NUMBER
		prev_img_num <- img_num

		# INITIATE SHINY APP
		#run_app <- runApp(app_dir, port = NULL, host = "127.0.0.1", launch.browser = TRUE, display.mode = "auto")
		run_app <- runApp(app_dir)
		#run_app <- 'exit'
		
		# REMOVE ANY IMAGE FILES IN IMG FOLDER
		if(length(list.files(paste0(app_dir, '/www/img/'))) > 0)
			file.remove(paste0(app_dir, '/www/img/', list.files(paste0(app_dir, '/www/img/'))))

		# DETERMINE WHETHER TO CHANGE IMAGE OR QUIT
		if(run_app$next.command == 'next'){
			img_num <- img_num + 1
		}else if (run_app$next.command == 'prev'){
			img_num <- img_num - 1
		}else{			
			return(NULL)
		}
		
		# MAKE SURE IMAGE NUMBER IS SENSIBLE
		if(img_num < 1) img_num <- 1
		if(img_num > length(image.file)) img_num <- length(image.file)
	}
}
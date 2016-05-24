process_digitize_images_input <- function(image.file = image.file, 
	shapes.file = NULL, landmarks.file = NULL, control.points.file = NULL, 
	curve.points.file = NULL, cal.file = cal.file, landmarks.ref=landmarks.ref, 
	curves.ref=curves.ref){

	images_fpaths <- NULL
	shapes_fpaths <- NULL
	landmarks_fpaths <- NULL
	control_points_fpaths <- NULL
	curve_points_fpaths <- NULL
	landmarks_ref <- NULL
	curves_ref <- NULL

	if(is.null(shapes.file) && is.null(landmarks.file) && is.null(control.points.file) && is.null(curve.points.file)) 
		cat("Warning: 'shapes.file', 'landmarks.file', 'control.points.file' and 'curve.points.file' are all NULL. If you would like to save landmarks or curves, please specify a location to save either the landmarks and/or curves.\n")

	# DIRECTORY INPUT
	if(length(image.file) == 1 && !grepl('[.][a-zA-Z]+$', image.file[1])){
	
		# CHECK THAT SHAPES FILE IS ALSO A DIRECTORY
		if(!is.null(shapes.file) && grepl('[.]txt$', shapes.file, ignore.case=TRUE)) stop("If 'image.file' is a directory, 'shapes.file' should also be a directory.")
		if(!is.null(landmarks.file) && grepl('[.][a-zA-Z]+$', landmarks.file)) stop("If 'image.file' is a directory, 'landmarks.file' should also be a directory.")
		if(!is.null(control.points.file) && grepl('[.][a-zA-Z]+$', control.points.file)) stop("If 'image.file' is a directory, 'control.points.file' should also be a directory.")
		if(!is.null(curve.points.file) && grepl('[.][a-zA-Z]+$', curve.points.file)) stop("If 'image.file' is a directory, 'curve.points.file' should also be a directory.")

		# MAKE SURE IMAGE FILE FOLDER EXISTS
		if(!file.exists(image.file)[1]) stop(paste0("The 'image.file' folder, '", image.file, "', was not found."))

		# IF SHAPES.FILE DOESNT EXIST, CREATE IT
		if(!is.null(shapes.file) && !file.exists(shapes.file)) dir.create(shapes.file)

		# MAKE SURE THAT DIRECTORIES EXIST
		if(!is.null(landmarks.file) && !file.exists(landmarks.file)) stop(paste0("'landmarks.file' (", landmarks.file, ") not found."))
		if(!is.null(control.points.file) && !file.exists(control.points.file)) stop(paste0("'control.points.file' (", control.points.file, ") not found."))
		if(!is.null(curve.points.file) && !file.exists(curve.points.file)) stop(paste0("'curve.points.file' (", curve.points.file, ") not found."))

		# GET CONTENTS OF IMAGE FOLDER
		image_fdir <- list.files(image.file)
		
		# DETERMINE WHETHER THEY ARE DIRECTORIES OR FILES

		if(grepl('[.][a-zA-Z]+$', image_fdir[1])){

			## SINGLE VIEW CASE
			# CREATE MATRIX OF IMAGE FILES
			images_fpaths <- matrix(image_fdir, nrow=length(image_fdir), ncol=1)

		}else{

			## STEREO CASE
			# CHECK THAT EACH FOLDER IN THE IMAGE FOLDER HAS THE SAME NUMBER OF FILES
			number_images <- rep(NA, length(image_fdir))
			for(i in 1:length(image_fdir)) number_images[i] <- length(list.files(paste0(image.file, '/', image_fdir[i])))
			if(sd(number_images) > 0) stop(paste0("When digitizing stereo image sets each folder within 'image.file' ('", image.file, "') must have the same number of files."))

			# IF SUB-FOLDERS OF SHAPES.FILE DO NOT MATCH IMAGE_FDIR, CREATE THEM
			if(!is.null(shapes.file) && sum(!image_fdir %in% list.files(shapes.file)) > 0)
				for(i in 1:length(image_fdir)) if(!file.exists(paste0(shapes.file, '/', image_fdir[i]))) dir.create(paste0(shapes.file, '/', image_fdir[i]))

			# CHECK THAT FOLDERS IN SHAPES MATCH IMAGE FOLDERS
			if(!is.null(landmarks.file) && sum(!image_fdir %in% list.files(landmarks.file)) > 0) stop(paste0("The landmarks folder (", landmarks.file, ") does not contain the same folders as the image.file (", image.file, ")."))
			if(!is.null(control.points.file) && sum(!image_fdir %in% list.files(control.points.file)) > 0) stop(paste0("The control points folder (", control.points.file, ") does not contain the same folders as the image.file (", image.file, ")."))
			if(!is.null(curve.points.file) && sum(!image_fdir %in% list.files(curve.points.file)) > 0) stop(paste0("The curve points folder (", curve.points.file, ") does not contain the same folders as the image.file (", image.file, ")."))

			# CREATE MATRIX OF IMAGE FILES
			images_fpaths <- matrix(NA, nrow=number_images, ncol=length(image_fdir))
		
			# FILL MATRIX WITH FILE NAMES ONLY
			for(i in 1:length(image_fdir)) images_fpaths[, i] <- list.files(paste0(image.file, '/', image_fdir[i]))
		
			# CHECK THAT NAMES MATCH ACROSS ROWS (VIEWS)
			found_in_all <- matrix(FALSE, nrow=nrow(images_fpaths), ncol=ncol(images_fpaths))
			for(i in 1:nrow(images_fpaths)){
				for(j in 1:ncol(images_fpaths)){
					if(images_fpaths[i, j] %in% images_fpaths[i, (1:ncol(images_fpaths))[(1:ncol(images_fpaths)) != j]]) found_in_all[i, j] <- TRUE
				}
			}

			# REPORT ERROR IF ANY FILE NAMES DO NOT MATCH
			if(sum(!found_in_all) > 0) stop(paste0("The contents of each folder in 'image.file' (", image.file, ") are not identical. File names within each view folder must match exactly across all views."))

			# ADD FILE PATHS TO EACH FILE
			for(i in 1:ncol(images_fpaths)) images_fpaths[, i] <- paste0(image_fdir[i], '/', images_fpaths[, i])
		}

		# CREATE MATCHING FILES FOR SHAPE DATA
		if(!is.null(shapes.file)) shapes_fpaths <- gsub('.[a-zA-Z]+$', '.txt', images_fpaths)
		if(!is.null(landmarks.file)) landmarks_fpaths <- gsub('.[a-zA-Z]+$', '.txt', images_fpaths)
		if(!is.null(control.points.file)) control_points_fpaths <- gsub('.[a-zA-Z]+$', '.txt', images_fpaths)
		if(!is.null(curve.points.file)) curve_points_fpaths <- gsub('.[a-zA-Z]+$', '.txt', images_fpaths)

		# ADD DIRECTORY PREFIX
		for(i in 1:ncol(images_fpaths)) images_fpaths[, i] <- paste0(image.file, '/', images_fpaths[, i])
		if(!is.null(shapes.file)) for(i in 1:ncol(shapes_fpaths)) shapes_fpaths[, i] <- paste0(shapes.file, '/', shapes_fpaths[, i])
		if(!is.null(landmarks.file)) for(i in 1:ncol(landmarks_fpaths)) landmarks_fpaths[, i] <- paste0(landmarks.file, '/', landmarks_fpaths[, i])
		if(!is.null(control.points.file)) for(i in 1:ncol(control_points_fpaths)) control_points_fpaths[, i] <- paste0(control.points.file, '/', control_points_fpaths[, i])
		if(!is.null(curve.points.file)) for(i in 1:ncol(curve_points_fpaths)) curve_points_fpaths[, i] <- paste0(curve.points.file, '/', curve_points_fpaths[, i])
	}
	
	# VECTOR INPUT, NOT DIRECTORY
	if(is.vector(image.file) && grepl('[.][a-zA-Z]+$', image.file[1])){

		# CHECK THAT SHAPES FILE IS ALSO A VECTOR
		if(!is.null(shapes.file) && !is.vector(shapes.file)) stop("If 'image.file' is a vector, 'shapes.file' must also be a vector.")
		if(!is.null(landmarks.file) && !is.vector(landmarks.file)) stop("If 'image.file' is a vector, 'landmarks.file' must also be a vector.")
		if(!is.null(control.points.file) && !is.vector(control.points.file)) stop("If 'image.file' is a vector, 'control.points.file' must also be a vector.")
		if(!is.null(curve.points.file) && !is.vector(curve.points.file)) stop("If 'image.file' is a vector, 'curve.points.file' must also be a vector.")

		# SET VECTORS AS SINGLE COLUMN MATRICES
		images_fpaths <- matrix(image.file, ncol=1)
		shapes_fpaths <- matrix(shapes.file, ncol=1)
		if(!is.null(landmarks.file)) landmarks_fpaths <- matrix(landmarks.file, ncol=1)
		if(!is.null(control.points.file)) control_points_fpaths <- matrix(control.points.file, ncol=1)
		if(!is.null(curve.points.file)) curve_points_fpaths <- matrix(curve.points.file, ncol=1)
	}

	# MATRIX INPUT
	if(is.matrix(image.file)){

		# CHECK THAT SHAPES FILE IS ALSO A MATRIX
		if(!is.null(shapes.file) && !is.matrix(shapes.file)) stop("If 'image.file' is a matrix, 'shapes.file' must also be a matrix.")
		if(!is.null(landmarks.file) && !is.matrix(landmarks.file)) stop("If 'image.file' is a matrix, 'landmarks.file' must also be a matrix.")
		if(!is.null(control.points.file) && !is.matrix(control.points.file)) stop("If 'image.file' is a matrix, 'control.points.file' must also be a matrix.")
		if(!is.null(curve.points.file) && !is.matrix(curve.points.file)) stop("If 'image.file' is a matrix, 'curve.points.file' must also be a matrix.")

		# RENAME TO MATCH PROCESSED DIRECTORY INPUT
		images_fpaths <- image.file
		shapes_fpaths <- shapes.file
		landmarks_fpaths <- landmarks.file
		control_points_fpaths <- control.points.file
		curve_points_fpaths <- curve.points.file
	}
	
	# MAKE SURE THAT NUMBER OF ELEMENTS MATCHES ACROSS IMAGE AND SHAPE FILES
	if(!is.null(shapes_fpaths) && nrow(shapes_fpaths) != nrow(images_fpaths)) stop(paste0("The length of 'shapes.file' (", nrow(shapes_fpaths), ") does not match the length of 'image.file' (", nrow(images_fpaths), ")."))
	if(!is.null(landmarks_fpaths) && nrow(landmarks_fpaths) != nrow(images_fpaths)) stop(paste0("The length of 'landmarks.file' (", nrow(landmarks_fpaths), ") does not match the length of 'image.file' (", nrow(images_fpaths), ")."))
	if(!is.null(control_points_fpaths) && nrow(control_points_fpaths) != nrow(images_fpaths)) stop(paste0("The length of 'control.points.file' (", nrow(control_points_fpaths), ") does not match the length of 'image.file' (", nrow(images_fpaths), ")."))
	if(!is.null(curve_points_fpaths) && nrow(curve_points_fpaths) != nrow(images_fpaths)) stop(paste0("The length of 'curve.points.file' (", nrow(curve_points_fpaths), ") does not match the length of 'image.file' (", nrow(images_fpaths), ")."))
	
	# GET CALIBRATION COEFFICIENTS
	cal_coeffs <- NULL
	if(!is.null(cal.file)){
		
		cal_list <- XML4R2list(cal.file)
		
		if(is.list(cal_list)){
			if('calibration' %in% names(cal_list)) cal_coeffs <- cal_list$calibration$cal.coeff
			if('cal.coeff' %in% names(cal_list)) cal_coeffs <- cal_list$cal.coeff
			if(is.null(cal_coeffs)) stop(paste0("Calibration coefficients not found in '", cal.file, "'."))
		}else{

			# SPLIT AT TABS
			tab_split <- strsplit(cal_list, '\t')

			# READ IN AS MATRIX
			cal_coeffs <- matrix(NA, nrow=length(tab_split), ncol=length(tab_split[[1]]))
			
			# FILL MATRIX
			for(i in 1:length(tab_split)) cal_coeffs[i, ] <- as.numeric(tab_split[[i]])
		}
		
		# CHECK THAT NUMBER OF COLUMNS IN COEFFICIENT MATRIX MATCHES NUMBER IN IMAGE FPATHS
		if(ncol(cal_coeffs) != ncol(images_fpaths)) stop(paste0("The number of columns in cal_coeffs (", ncol(cal_coeffs), ") does not match the number of image views (", ncol(images_fpaths), ")."))
	}

	# CHECK IMAGE EXTENSIONS FOR COMPATIBILITY
	check_img_type <- grepl(pattern='[.]jpg$|[.]jpeg$|[.]tif$|[.]tiff$|[.]png$', x=images_fpaths, ignore.case=TRUE)
	if(sum(!check_img_type) > 0) stop(paste0("Only images of type JPG, JPEG, TIF, TIFF and PNG are currently supported. The following are unsupported file types:\n\t", paste(images_fpaths[!check_img_type], collapse="\n\t")))

	# READ IN LANDMARKS
	landmarks_ref <- NULL
	if(!is.null(landmarks.ref)){
		
		# SINGLE ELEMENT IN VECTOR AND EXISTS AS A FILE
		if(length(landmarks.ref) == 1){

			# IS FILE
			if(grepl('[.]txt$', landmarks.ref)){
				if(!file.exists(landmarks.ref)) stop(paste0("landmarks.ref ('", landmarks.ref, "') not found."))
				landmarks.ref <- as.vector(suppressWarnings(read.table(landmarks.ref, sep="\n"))[,1])
			}
		}

		landmarks_ref <- landmarks.ref
	}

	# MAKE SURE CURVE REFERENCE IS FILE OR MATRIX IF NOT NULL
	if(!is.null(curves.ref)){
		if(!grepl('[.]txt$', curves.ref[1]) && !is.matrix(curves.ref)) stop("'curves.ref' must either be a matrix or a .txt file containing a tab-delimited matrix of curve names, start points and end points.")
		if(is.matrix(curves.ref) && ncol(curves.ref) != 3) stop(paste0("If 'curves.ref' is a matrix, it must have three columns. The current curves.ref input has ", ncol(curves.ref), " column(s)."))

		# READ IN CURVE REF
		curves_ref <- matrix(NA, nrow=0, ncol=3)
		landmarks_from_curves_ref <- rep(NA, 0)
		if(!is.null(curves.ref)){

			# NOT A MATRIX - ASSUME TO BE FILE PATH
			if(!is.matrix(curves.ref)){
				if(!file.exists(curves.ref)) stop(paste0("curves.ref ('", curves.ref, "') not found."))
				curves.ref <- suppressWarnings(as.matrix(read.table(curves.ref, sep="\t")))
			}

			colnames(curves.ref) <- NULL
			for(i in 1:nrow(curves.ref)){
				curves_ref <- rbind(curves_ref, curves.ref[i, ])
				landmarks_from_curves_ref <- c(landmarks_from_curves_ref, curves.ref[i, 2:3])
			}
		}

		# ADD LANDMARKS FROM CURVES
		if(!is.null(landmarks.file) || !is.null(shapes.file)){
			landmarks_from_curves_ref <- unique(landmarks_from_curves_ref)
			landmarks_from_curves_ref <- landmarks_from_curves_ref[!landmarks_from_curves_ref %in% landmarks_ref]
			landmarks_ref <- c(landmarks_ref, landmarks_from_curves_ref)
		}
	}

	rlist <- list(
		images_fpaths = images_fpaths, 
		shapes_fpaths = shapes_fpaths,
		landmarks_fpaths = landmarks_fpaths, 
		control_points_fpaths = control_points_fpaths,
		curve_points_fpaths = curve_points_fpaths, 
		cal_coeffs = cal_coeffs,
		landmarks_ref = landmarks_ref,
		curves_ref = curves_ref
	)

	# CONVERT MATRICES TO LISTS/VECTORS IN ORDER TO PASS THROUGH JSON
	for(i in 1:length(rlist)){

		if(!is.matrix(rlist[[names(rlist)[i]]])) next

		if(ncol(rlist[[names(rlist)[i]]]) == 1){
			rlist[[names(rlist)[i]]] <- c(rlist[[names(rlist)[i]]])
			next
		}

		matrix_to_list <- list()
		for(j in 1:nrow(rlist[[names(rlist)[i]]])) matrix_to_list[[j]] <- rlist[[names(rlist)[i]]][j, ]
		rlist[[names(rlist)[i]]] <- matrix_to_list
	}
	
	rlist
}
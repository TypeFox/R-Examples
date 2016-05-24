testCalibration <- function(img.dir, cal.file, corner.dir, sq.size, nx, ny, plot.dir = NULL,
	verify.dir = NULL, print.progress = TRUE){

	# READ CALIBRATION FILE
	cal.list <- XML4R2list(file=cal.file)$calibration

	# GET CALIBRATION DATA
	cal_coeff <- cal.list$cal.coeff
	img_type <- cal.list$img.type

	# READ IMAGE DIRECTORY
	imgs_list_files <- list.files(img.dir)
	
	# SET NUMBER OF VIEWS
	num_views <- length(imgs_list_files)

	# GET FILE PATHS TO CALIBRATON IMAGE SUB-FOLDERS
	img_fpaths <- paste0(img.dir, '/', imgs_list_files)

	# GET SUB-FOLDER NAMES
	img_sub_dir <- gsub('(.mov|.avi|.mp4|.mpg)$', '', imgs_list_files, ignore.case=TRUE)

	# GET NUMBER OF SPACES TO ALIGN RIGHT OF SUB DIRECTORY NAMES
	img_sub_dir_salign <- (1 + max(nchar(img_sub_dir))) - nchar(img_sub_dir)

	# CHECK IF PLOT.DIR EXISTS
	if(!is.null(plot.dir)) if(!file.exists(plot.dir)) dir.create(path=plot.dir)

	# CHECK IF CORNER FOLDERS EXIST
	if(!file.exists(corner.dir)) dir.create(path=corner.dir)

	# CREATE SUB FOLDERS IF NOT PRESENT
	for(dir_name in img_sub_dir) if(!file.exists(paste0(corner.dir, '/', dir_name))) dir.create(paste0(corner.dir, '/', dir_name))

	# GET SQUARE SIZES AND UNITS
	sq.size.num <- as.numeric(gsub('[[:alpha:], ]', '', sq.size))
	sq.size.units <- gsub('[[:digit:]., ]', '', sq.size)

	# SET DEFAULTS
	img_fnames <- NULL
	vid_fnames <- NULL
	verify_fpaths <- NULL

	if(img_type == 'image'){

		# GET IMAGE NAMES
		img_fnames_v1 <- list.files(paste0(img.dir, '/', imgs_list_files[1]))
		img_fnames_v2 <- list.files(paste0(img.dir, '/', imgs_list_files[2]))

		# FIND COMMON CALIBRATION IMAGE FILENAMES
		img_fnames <- img_fnames_v1[img_fnames_v1 %in% img_fnames_v2]

	}else if(img_type == 'video'){
	
		# GET CALIBRATION VIDEO NAMES
		vid_fnames <- setNames(imgs_list_files, gsub('[.][A-Za-z0-9]*$', '', imgs_list_files))
	}

	# CHECK IF VERIFY FOLDERS EXIST
	if(!is.null(verify.dir)){

		# CREATE VERIFY FOLDER IF DOES NOT EXIST
		if(!file.exists(verify.dir)) dir.create(verify.dir)

		# CREATE VIEW FOLDERS IF VERIFY FOLDER IS EMPTY
		for(dir_name in img_sub_dir) if(!file.exists(paste0(verify.dir, '/', dir_name))) dir.create(paste0(verify.dir, '/', dir_name))

		# SET VERIFY FILE PATHS AND FILE NAMES
		verify_fpaths <- paste0(verify.dir, '/', list.files(verify.dir))
		if(img_type == 'image') verify_fnames <- img_fnames
	}

	if(print.progress){
		cat("calibrateCameras\n\n")
	}

	## DETECT CORNERS
	if(print.progress) cat("\tTest checkerboard corner detection...")

	# CHECK IF CORNERS ARE ALREADY FOUND
	detect_corners <- TRUE

	# CHECK IF CORNERS ARE FOUND IN ALL VIEWS
	list_files_length <- rep(NA, num_views)
	for(i in 1:num_views) list_files_length[i] <- length(list.files(paste0(corner.dir, '/', img_sub_dir[i])))

	# IF ALL FOLDERS HAVE FILES, PROMPT WHETHER TO RE-DETECT CORNERS
	if(sum(list_files_length > 0) == num_views){

		if(print.progress){
			cat(paste0("Saved calibration corners found for all views in '", corner.dir, "' folder.\n\n"))

			for(i in 1:num_views){
				num_frames_detected <- length(list.files(paste0(corner.dir, '/', img_sub_dir[i])))
				if(img_type == 'image'){
					cat(paste0("\t\t", img_sub_dir[i], paste(rep(' ', img_sub_dir_salign[i]), collapse='')))
					cat(paste0(": Corners detected in ", num_frames_detected, " aspects\n"))
				}else{
					cat(paste0("\t\t", vid_fnames[i], paste(rep(' ', img_sub_dir_salign[i]), collapse='')))
					cat(paste0(": Corners detected in ", num_frames_detected, " frames\n"))
				}
			}
		}else{
			if(print.progress) cat('\n\n')
		}

		detect_corners <- FALSE

		response <- readline(prompt="\n\t\tDo you wish to re-detect the calibration corners? (y/n) : ");
		if(print.progress) cat("\n")
		#response <- 'n'

		if(tolower(response) %in% c('yes', 'y')) detect_corners <- TRUE
	}else{
		if(print.progress) cat('\n')
	}

	## DETECT CALIBRATION CHECKERBOARD CORNERS
	if(detect_corners){

		if(img_type == 'image'){

			# CREATE ARRAY FOR CORNERS
			cal_corners <- array(NA, dim=c(nx*ny, 2, length(img_fnames), num_views))

			# FIND CHECKERBOARD CORNERS AND READ INTO ARRAY
			if(print.progress) cat("\t\tRunning automated checkerboard corner detection on calibration images...\n")
			for(i in 1:length(img_fnames)){

				# FIND IMAGE FILES IN EACH VIEW FOLDER
				if(print.progress) cat("\t\t\t", img_fnames[i], "\n", sep="")

				for(j in 1:num_views){
		
					# SET VERIFY FILEPATH IF NON-NULL
					verify_fpath <- NULL
					if(!is.null(verify_fpaths)) verify_fpath <- paste0(verify_fpaths[j], '/', verify_fnames[i])

					# SET CORNER FILEPATH
					corner_fpath <- paste0(corner.dir, '/', img_sub_dir[j], '/', gsub('[.][A-Za-z]+$', '.txt', img_fnames[i]))
				
					# SPECIFY WHETHER TO FLIP CORNER ORDER
					if(cal.list$flip.view && j == 2){flip <- TRUE}else{flip <- FALSE}

					cal_corners[, , i, j] <- findCheckerboardCorners(image.file=paste0(img_fpaths[j], '/', img_fnames[i]), 
						nx=nx, ny=ny, flip=flip, corner.file=corner_fpath, verify.file=verify_fpath, print.progress=FALSE)

					if(print.progress){
						cat("\t\t\t\tView ", j, " : ", sep="")
						if(is.na(cal_corners[1, 1, i, j])){cat("findCheckerboardCorners() unsuccessful")}else{cat(nx*ny, " corners found", sep="")}
						cat("\n")
					}
				}
			}

			# SET DIMNAMES FOR CORNER ARRAY
			dimnames(cal_corners) <- list(NULL, NULL, gsub('[.][A-Za-z0-9]*$', '', img_fnames), img_sub_dir)

			if(print.progress) cat('\n')
		}
	}

	# READ IN CORNERS IF NOT DETECTED FROM IMAGES
	if(!detect_corners || img_type == 'video'){

		# GET FRAME NAMES
		frame_names <- c()
		for(i in 1:num_views) frame_names <- c(frame_names, gsub('.txt', '', list.files(paste0(corner.dir[i], '/', img_sub_dir[i]))))

		# GET UNIQUE FRAME NAMES
		frame_names_unique <- unique(frame_names)

		# MAKE CORNER ARRAY
		cal_corners <- array(NA, dim=c(nx*ny, 2, length(frame_names_unique), num_views), 
			dimnames=list(NULL, NULL, frame_names_unique, img_sub_dir))

		# FILL CORNER ARRAY
		for(i in 1:num_views){

			# GET CORNER FILES
			corner_files <- gsub('.txt', '', list.files(paste0(corner.dir, '/', img_sub_dir[i])))
	
			# READ CORNERS INTO ARRAY FROM FILES
			for(j in 1:length(corner_files)){
				cal_corners[, , corner_files[j], i] <- 
					as.matrix(read.table(paste0(corner.dir, '/', img_sub_dir[i], '/', corner_files[j], '.txt')))
			}
		}
	}

	# FIND NUMBER OF NON-NA VIEWS FOR EACH ASPECT
	aspect_non_na <- rowSums(apply(!is.na(cal_corners), c(3, 4), 'sum') > 0)

	# FIND PAIRS WITH CORNERS IN MINIMUM NUMBER OF VIEWS
	cal_min_views_found <- aspect_non_na >= 2
	
	# FIND TOTAL NUMBER OF ASPECTS FOUND FOR CALIBRATION
	cal_views_found_num <- sum(cal_min_views_found)

	# FIND NUMBER OF NON-NA VIEWS FOR EACH ASPECT AND VIEW
	non_na_by_view <- apply(!is.na(cal_corners), c(3, 4), 'sum') > 0

	if(print.progress){
		cat(paste0("\t\tNumber of cases in which corners were found in at least ", 2, " views: ", cal_views_found_num, "\n"))
		if(img_type == 'image') cat(paste0("\t\t\tFilenames: ", paste0(gsub('.jpeg|.jpg|.tiff', '', img_fnames[cal_min_views_found], ignore.case=TRUE), collapse=", "), "\n"))
	}
	
	# REMOVE NA PAIRS
	cal_corners <- cal_corners[, , cal_min_views_found, ]

	# CREATE ERROR PLOTS
	dlt_test <- createErrorPlots(cal.coeff=cal_coeff, corners=cal_corners, nx=nx, sq.size.num=sq.size.num, 
		sq.size.units=sq.size.units, file=plot.dir)

	# PRINT ERROR SUMMARY
	if(print.progress) print(summary(dlt_test, print.tab='\t'))
}
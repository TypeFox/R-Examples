calibrateCameras <- function(img.dir, sq.size, nx, ny, cal.file, corner.dir,
	print.progress = TRUE, flip.view = FALSE, verify.dir = NULL, 
	min.views = 'max', exec.dir = '', undistort = FALSE, num.aspects.read = 'auto', 
	num.sample.est = 'auto', num.sample.sets = 'auto', num.aspects.sample = 'auto', 
	max.sample.optim = 30, nlm.calls.max = 20, fit.min.break = 1, objective.min = 1, 
	objective.min.break = 5, with.circles = FALSE, sample.est = NULL, ...){

	################################ CHECK INPUT PARAMETERS ##############################

	## READ CALIBRATION FILE, IF EXISTS
	cal.list <- list()
	if(file.exists(cal.file)) cal.list <- XML4R2list(file=cal.file)$calibration

	# FIND PATH TO FOLDER WHERE CALIBRATION FILE IS
	cal_file_str_split <- strsplit(cal.file, '/')[[1]]
	
	# SET CALIBRATION DIRECTORY BY REMOVING FILENAME
	if(length(cal_file_str_split) > 1){calib_dir <- paste0(paste(head(cal_file_str_split, -1), collapse="/"), "/")}else{calib_dir <- ""}
	
	# SET INPUT PARAMETERS TO OVERWRITE FROM CALIBRATION FILE (IF NON-NULL IN FILE)
	write_param_from_file <- c('img.dir', 'sq.size', 'nx', 'ny', 'corner.dir', 'flip.view', 'verify.dir')
	
	# OVERWRITE ANY NULL INPUT PARAMETERS WITH VALUE IN CALIBRATION FILE (IF NON-NULL IN FILE)
	for(write_param in write_param_from_file)
		if(is.null(get(write_param)) && !is.null(cal.list[[write_param]])) assign(write_param, cal.list[[write_param]])

	# GET SQUARE SIZES AND UNITS
	sq.size.num <- as.numeric(gsub('[[:alpha:], ]', '', sq.size))
	sq.size.units <- gsub('[[:digit:]., ]', '', sq.size)

	# CHECK CALIBRATION IMAGE INPUTS UNLESS CAL FILE IMGS MATCHES INPUT
	check_img_input <- TRUE
	if(!is.null(cal.list[['img.dir']])) if(cal.list[['img.dir']] == img.dir) check_img_input <- FALSE

	# CHECK CALIBRATION IMAGE INPUTS
	if(is.null(cal.list$img.list.files)){
		if(!file.exists(img.dir)) stop(paste0("Folder '", img.dir, "' not found."))
		if(length(list.files(img.dir)) == 0) stop(paste0("No files/folders found in '", img.dir, "'."))
		if(length(list.files(img.dir)) == 1) stop(paste0("Only one file/folder found in '", img.dir, "'. Two views are required for stereo calibration."))
	}
	
	# GET FILE NAMES IN CALIBRATON IMAGE FOLDER
	if(file.exists(img.dir) && length(list.files(img.dir)) > 0){
		imgs_list_files <- list.files(img.dir)
	}else{
		imgs_list_files <- cal.list$img.list.files
	}

	# GET FILE PATHS TO CALIBRATON IMAGE SUB-FOLDERS
	img_fpaths <- paste0(img.dir, '/', imgs_list_files)
	
	# SET NUMBER OF CAMERA VIEWS
	num_views <- length(img_fpaths)
	
	# SET MIN NUMBER OF VIEWS TO USE IN ESTIMATING CALIBRATION COEFFICIENTS
	if(min.views == 'max'){min_views <- max(num_views)}else{min_views <- min.views}

	# CHECK THAT CAL MIN VIEWS IS EQUAL TO OR GREATER THAN TWO
	if(min_views < 2) stop("'min.views' must be greater or equal to 2.")
	
	# CHECK IF AVI FILES
	if(sum(grepl('(.mov|.avi|.mp4|.mpg)$', img_fpaths, ignore.case=TRUE)) > 0){img_type <- 'video'}else{img_type <- 'image'}

	# GET SUB-FOLDER NAMES
	img_sub_dir <- gsub('(.mov|.avi|.mp4|.mpg)$', '', imgs_list_files, ignore.case=TRUE)

	# GET NUMBER OF SPACES TO ALIGN RIGHT OF SUB DIRECTORY NAMES
	img_sub_dir_salign <- (1 + max(nchar(img_sub_dir))) - nchar(img_sub_dir)

	# SET DEFAULTS
	img_fnames <- NULL
	vid_fnames <- NULL
	vid_nframes <- NULL
	img_size <- NULL

	if(img_type == 'image'){

		# GET IMAGE NAMES
		img_fnames_v1 <- list.files(paste0(img.dir, '/', imgs_list_files[1]))
		img_fnames_v2 <- list.files(paste0(img.dir, '/', imgs_list_files[2]))

		# FIND COMMON CALIBRATION IMAGE FILENAMES
		img_fnames <- img_fnames_v1[img_fnames_v1 %in% img_fnames_v2]

		# CHECK THAT THERE ARE AT LEAST TWO IMAGES COMMON ASPECTS
		if(length(img_fnames) < 2) stop(paste0("Number of common images between two views in '", img.dir, "' is less than two. A minimum of two images are required for stereo calibration and more than five are recommended."))

		# GET IMAGE SIZE IN EACH VIEW
		if(undistort){
			stop("Undistortion not yet available for image input.")
		}

	}else if(img_type == 'video'){
	
		# GET CALIBRATION VIDEO NAMES
		vid_fnames <- setNames(imgs_list_files, gsub('[.][A-Za-z0-9]*$', '', imgs_list_files))

		# CHECK EXEC PATH EXISTS
		if(exec.dir != '' && !file.exists(exec.dir)) stop(paste0("exec.dir '", exec.dir, "' not found."))

		# CHECK THAT EXECUTABLES ARE IN FOLDER
		if(!'get_frame_count' %in% list.files(exec.dir)) stop(paste0("'get_frame_count' not found in '", exec.dir, "'."))
		if(!'find_checkerboard_corners' %in% list.files(exec.dir)) stop(paste0("'find_checkerboard_corners' not found in '", exec.dir, "'."))

		# ADD SLASH AT END IF NOT PRESENT
		if(!grepl('[/]$', exec.dir)) exec.dir <- paste0(exec.dir, '/')

		# GET NUMBER OF FRAMES IN EACH VIDEO
		if(!is.null(cal.list$vid.nframes)){

			vid_nframes <- cal.list$vid.nframes
		}else{

			vid_nframes <- rep(NA, length(vid_fnames))

			for(i in 1:length(vid_fnames)){

				# SET COMMAND TO FIND FRAME COUNT
				command <- paste0('./', gsub(' ', '\\\\ ', exec.dir), 'get_frame_count ', 
					gsub(' ', '\\\\ ', img.dir), '/', gsub(' ', '\\\\ ', vid_fnames[i]))

				# FIND NUMBER OF FRAMES IN VIDEO
				vid_nframes[i] <- as.numeric(system(command=command, intern=TRUE))
			}
		}

		# GET IMAGE SIZE IN EACH VIEW
		if(undistort){
			if(!is.null(cal.list$img.size)){
				img_size <- cal.list$img.size
			}else{

				img_size <- matrix(NA, nrow=num_views, ncol=2, dimnames=list(img_sub_dir, c('w', 'h')))

				for(i in 1:length(vid_fnames)){

					# SET COMMAND TO FIND FRAME COUNT
					command <- paste0('./', gsub(' ', '\\\\ ', exec.dir), 'get_frame_size ', 
						gsub(' ', '\\\\ ', img.dir), '/', gsub(' ', '\\\\ ', vid_fnames[i]))

					# FIND NUMBER OF FRAMES IN VIDEO
					img_size[i, ] <- as.numeric(strsplit(x=system(command=command, intern=TRUE), split=',')[[1]])
				}
			}
		}
	}
	
	# SET DEFAULTS
	verify_fpaths <- NULL
	verify_fnames <- NULL

	# IF VIDEO, MAKE SURE CAL VERIFY IS MADE
	if(img_type == 'video' && is.null(verify.dir)) verify.dir <- 'Corner detection'

	# CHECK IF CORNERS FOLDER EXISTS
	if(!file.exists(corner.dir)) dir.create(path=corner.dir)

	# CREATE SUB FOLDERS IF NOT PRESENT
	for(dir_name in img_sub_dir) if(!file.exists(paste0(corner.dir, '/', dir_name))) dir.create(paste0(corner.dir, '/', dir_name))

	# CHECK CAL VERIFY IF NON-NULL
	if(!is.null(verify.dir)){

		# CREATE VERIFY FOLDER IF DOES NOT EXIST
		if(!file.exists(verify.dir)) dir.create(verify.dir)

		# CREATE VIEW FOLDERS IF VERIFY FOLDER IS EMPTY
		for(dir_name in img_sub_dir) if(!file.exists(paste0(verify.dir, '/', dir_name))) dir.create(paste0(verify.dir, '/', dir_name))

		# SET VERIFY FILE PATHS AND FILE NAMES
		verify_fpaths <- paste0(verify.dir, '/', list.files(verify.dir))
		if(img_type == 'image') verify_fnames <- img_fnames
	}

	# SAVE TO CAL.LIST
	cal.list[['img.type']] <- img_type
	cal.list[['img.size']] <- img_size
	cal.list[['img.fpaths']] <- img_fpaths
	cal.list[['img.sub.dir']] <- img_sub_dir
	cal.list[['img.fnames']] <- img_fnames
	cal.list[['img.list.files']] <- imgs_list_files
	cal.list[['min.views']] <- min_views
	cal.list[['vid.fnames']] <- vid_fnames
	cal.list[['vid.nframes']] <- vid_nframes
	cal.list[['vid.nframes']] <- vid_nframes

	# WRITE INPUT PARAMETERS TO CAL.LIST
	for(write_param in write_param_from_file) cal.list[[write_param]] <- get(write_param)

	# PRINT INPUT PARAMETERS
	if(print.progress){
		cat("calibrateCameras\n\n")
		cat("\tChecking input parameters...\n")
		cat(paste0("\t\tCalibration file: '", cal.file, "'\n"))
		cat(paste0("\t\tCalibration input type: ", img_type, "\n"))
		cat(paste0("\t\tNumber of camera views: ", num_views, "\n"))
		if(num_views == 2) cat(paste0("\t\tOne view upside-down?: ", flip.view, "\n"))
		
		if(img_type == 'image'){
			cat(paste0("\t\tCalibration image set:\n"))
			cat(paste0("\t\t\tNumber of common images per view: ", length(img_fnames), "\n"))
			cat(paste0("\t\t\t\tFilenames: ", paste0(gsub('.jpeg|.jpg|.tiff', '', img_fnames, ignore.case=TRUE), collapse=", "), "\n"))
		}else if(img_type == 'video'){
			cat(paste0("\t\tCalibration videos:\n"))
			cat(paste0("\t\t\tFilenames:\n"))
			for(i in 1:length(vid_fnames)){
				cat(paste0("\t\t\t\t", vid_fnames[i], paste(rep(' ', img_sub_dir_salign[i]), collapse=''), 
					" (", vid_nframes[i], " frames)\n"))
			}
		}
		cat(paste0("\t\t\tSquare size: ", sq.size, "\n"))
		cat(paste0("\t\t\tInternal corners: ", nx, " x ", ny, " (", nx*ny, " total)\n"))
		
		cat("\n")
	}

	# SAVE CALIBRATION LIST TO FILE
	list2XML4R(list('calibration' = cal.list), file=cal.file)

	
	########################## CALIBRATION CHECKERBOARD DETECTION ########################
	if(print.progress) cat("\tCalibration checkerboard corner detection...")

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
		cat("\n")
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
					if(flip.view && j == 2){flip <- TRUE}else{flip <- FALSE}

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

		}else if(img_type == 'video'){

			if(num.aspects.read == 'auto') num.aspects.read <- 60

			# CHECK THAT NUMBER OF ASPECTS IN VIDEO EXCEEDS SAMPLE NUMBER
			if(min(vid_nframes) < num.aspects.read + 40){
				cat('\n')
				stop(paste0("The number of video frames (", min(vid_nframes), ") is less than the number of aspects to be sampled for checkerboard detection (", num.aspects.read + 40, ")."))
			}
		
			min_nframes <- min(vid_nframes)

			for(i in 1:num_views){

				if(print.progress) cat(paste0("\t\tDetecting corners in '", vid_fnames[i], "'..."))

				# WRITE COMMAND
				command <- paste0(
					'./', gsub(' ', '\\\\ ', exec.dir), 'find_checkerboard_corners ', 
					gsub(' ', '\\\\ ', img.dir), '/', gsub(' ', '\\\\ ', vid_fnames[i]), ' ', 
					gsub(' ', '\\\\ ', corner.dir), '/', gsub(' ', '\\\\ ', img_sub_dir[i]), ' ', 
					gsub(' ', '\\\\ ', verify.dir), '/', gsub(' ', '\\\\ ', img_sub_dir[i]), ' ', 
					nx, ' ', ny, ' ', num.aspects.read, ' ', num.aspects.read+40, ' ', 
					min_nframes, ' 0 ', as.numeric(with.circles))
				
				# CALL COMMAND
				#cat(command, '\n')
				system(command=command)
				
				if(print.progress){

					# NUMBER OF CORNERS DETECTED
					num_frames_detected <- length(list.files(paste0(corner.dir, '/', img_sub_dir[i])))

					cat(paste(rep(' ', img_sub_dir_salign[i]), collapse=''))
					cat(paste0(" Corners detected in ", num_frames_detected, " frames\n"))
				}
			}
		}
	}

	# READ IN CORNERS IF NOT DETECTED FROM IMAGES
	if(!detect_corners || img_type == 'video'){

		# GET FRAME NAMES
		frame_names <- c()
		for(i in 1:num_views) frame_names <- c(frame_names, gsub('.txt', '', list.files(paste0(corner.dir, '/', img_sub_dir[i]))))

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
	cal_min_views_found <- aspect_non_na >= min_views
	
	# FIND TOTAL NUMBER OF ASPECTS FOUND FOR CALIBRATION
	cal_views_found_num <- sum(cal_min_views_found)

	# FIND NUMBER OF NON-NA VIEWS FOR EACH ASPECT AND VIEW
	non_na_by_view <- apply(!is.na(cal_corners), c(3, 4), 'sum') > 0

	if(print.progress){
		cat(paste0("\t\tNumber of cases in which corners were found in at least ", min_views, " views: ", cal_views_found_num, "\n"))
		if(img_type == 'image') cat(paste0("\t\t\tFilenames: ", paste0(gsub('.jpeg|.jpg|.tiff', '', img_fnames[cal_min_views_found], ignore.case=TRUE), collapse=", "), "\n"))
		cat('\n')

		# THIS WILL ONLY APPLY WITH GREATER THAN TWO VIEWS (CURRENTLY WRITTEN ONLY TO WORK WITH 3 VIEWS)
		if(num_views > 2){
			
			# VIEW COMBINATIONS
			view_combos <- list(1:2, 2:3, c(1,3), 1:3)

			cat(paste0("\t\tNumber of aspects in which corners were found for all view combinations\n"))
			for(view_combo in view_combos){
				cat(paste0("\t\t\t", paste(img_sub_dir[view_combo], collapse=', '), ": ", sum(rowSums(non_na_by_view[, view_combo]) == length(view_combo)), "\n"))
			}

			cat('\n')
		}
	}

	######################### ESTIMATE UNDISTORTION COEFFICIENTS #########################

	if(undistort){

		if(print.progress) cat("\tEstimating undistortion parameters...")

		# CHECK IF UNDISTORTION PARAMETERS ARE ALREADY FOUND
		estimate_undistortion <- TRUE
		if(!is.null(cal.list$undistort.params)){

			# READ CALIBRATION CORNERS INTO ARRAY
			undistort_params <- cal.list$undistort.params

			if(print.progress) cat("Undistortion parameters found in the calibration file.\n\n")

			estimate_undistortion <- FALSE

			response <- readline(prompt="\t\tDo you wish to re-estimate the undistortion parameters? (y/n) : ");cat("\n")
			#response <- 'n'
		
			if(tolower(response) %in% c('yes', 'y')) estimate_undistortion <- TRUE
		}else{

			if(print.progress) cat('\n')
			undistort_params <- matrix(NA, nrow=num_views, ncol=7, dimnames=list(img_sub_dir, c('cx', 'cy', 'k1', 'k2', 'k3', 'p1', 'p2')))
		}

		# SET MAXIMUM NUMBER OF ASPECTS TO USE IN UNDISTORTION
		max_undist_sample <- 15
		
		if(estimate_undistortion){
		
			# FIND DISTORTION COEFFICIENTS BY VIEW
			for(view in 1:num_views){

				# FIND NON-NA ASPECTS
				undist_sample_nona <- (1:dim(cal_corners)[3])[!is.na(cal_corners[1, 1, , view])]
				
				# SELECT AMONG NON-NA ASPECTS, SPECIFIED NUMBER BUT NO MORE THAN LENGTH
				undist_sample_nona <- undist_sample_nona[1:min(max_undist_sample, length(undist_sample_nona))]
				#undist_sample_nona <- undist_sample_nona[round(seq(1, length(undist_sample_nona), length=min(max_undist_sample, length(undist_sample_nona))))]

				# ESTIMATE DISTORTION COEFFICIENTS
				dist_params <- estimateDistortion(coor.2d=cal_corners[, , undist_sample_nona, view], nx, 
					image.size=img_size[view, ])
				
				undistort_params[view, ] <- dist_params
			}

			cal.list[['undistort.params']] <- undistort_params

			# SAVE CALIBRATION LIST
			list2XML4R(list('calibration' = cal.list), file=cal.file)
		}

		if(print.progress){

			for(view in 1:num_views){

				# FIND NON-NA ASPECTS
				undist_sample_nona <- (1:dim(cal_corners)[3])[!is.na(cal_corners[1, 1, , view])]
				
				# SELECT AMONG NON-NA ASPECTS, SPECIFIED NUMBER BUT NO MORE THAN LENGTH
				undist_sample_nona <- undist_sample_nona[1:min(max_undist_sample, length(undist_sample_nona))]
				#undist_sample_nona <- undist_sample_nona[round(seq(1, length(undist_sample_nona), length=min(max_undist_sample, length(undist_sample_nona))))]

				# FIND HOMOGRAPHY ERROR OF ORIGINAL COORDINATES
				findHo <- findHomography(cal_corners[, , undist_sample_nona, view], nx=nx)
				#findHo <- findHomography(cal_corners[, , !is.na(cal_corners[1, 1, , view]), view], nx=nx)

				# UNDISTORT CORNERS USED IN UNDISTORTION SAMPLE
				coor_2d_u <- undistort(cal_corners[, , undist_sample_nona, view], image.size=img_size[view, ],
					center=undistort_params[view, 1:2], k=undistort_params[view, 3:5], p=undistort_params[view, 6:7])
				#coor_2d_u <- undistort(cal_corners[, , !is.na(cal_corners[1, 1, , view]), view], image.size=img_size[view, ],
				#	center=undistort_params[view, 1:2], k=undistort_params[view, 3:5], p=undistort_params[view, 6:7])

				# FIND HOMOGRAPHY ERROR OF UNDISTORTED COORDINATES
				findHu <- findHomography(coor_2d_u, nx=nx)

				cat(paste0("\t\t", img_sub_dir[view], "\n"))
				cat(paste0("\t\t\tParameters (cx, cy, k1): ", paste(undistort_params[view, ], collapse=', '), "\n"))
				cat(paste0("\t\t\tMean homography fit (original, after distortion correction): ", 
					round(mean(findHo$error), 2), " px, ", round(mean(findHu$error), 2), " px\n"))
				cat(paste0("\t\t\tMax homography fit (original, after distortion correction): ", 
					round(max(findHo$error), 2), " px, ", round(max(findHu$error), 2), " px\n"))
			}
			
			cat('\n')
		}

		# UNDISTORT ALL CORNERS
		for(view in 1:num_views) cal_corners[, , , view] <- undistort(cal_corners[, , , view], 
			image.size=img_size[view, ], center=undistort_params[view, 1:2], 
			k=undistort_params[view, 3:5], p=undistort_params[view, 6:7])

	}else{

		if(print.progress){

			for(view in 1:num_views){

				# FIND HOMOGRAPHY ERROR OF ORIGINAL COORDINATES
				#findHo <- findHomography(cal_corners[, , !is.na(cal_corners[1, 1, , view]), view], nx=nx)

				#cat(paste0("\t\t", img_sub_dir[view], "\n"))
				#cat(paste0("\t\t\tMean homography fit (original): ", 
				#	round(mean(findHo$error), 2), " px\n"))
				#cat(paste0("\t\t\tMax homography fit (original): ", 
				#	round(max(findHo$error), 2), " px\n"))
			}
			
			#cat('\n')
		}
	}

	######################## ESTIMATE DLT CALIBRATION COEFFICIENTS #######################
	cat("\tCalibration coefficient estimation...")

	# IF CALIBRATION COEFFICIENTS ALREADY IN LIST, ASK USER WHETHER TO RE-ESTIMATE COEFFICIENTS
	estimate_cal_coeffs <- TRUE
	if(!is.null(cal.list$cal.coeff)){
		
		# READ FROM LIST
		cal_set_num <- cal.list$cal.set.num
		cal_coeff <- cal.list$cal.coeff
		mean_reconstruct_rmse <- cal.list$mean.reconstruct.rmse
		coefficient_rmse <- cal.list$coefficient.rmse

		if(print.progress) cat("Saved calibration DLT coefficients found in the calibration file.\n\n")

		estimate_cal_coeffs <- FALSE
		
		response <- readline(prompt="\t\tDo you wish to re-estimate the calibration coefficients? (y/n) : ");cat('\n')
		#response <- 'n'
		
		if(tolower(response) %in% c('yes', 'y')) estimate_cal_coeffs <- TRUE

	}else{
		if(print.progress) cat('\n')
	}

	# IF NAMES OF ASPECTS TO BE SAMPLED ARE EXPLICITLY DEFINED, MAKE SURE LENGTH OF NUM SAMPLE MATCHES NUMBER OF ASPECTS
	if(!is.null(sample.est)) num.sample.est <- length(sample.est)

	if(is.null(sample.est)){

		# REMOVE IMAGE PAIRS WHERE AT LEAST ONE IS NA (NUMBER OF NON-NA IMAGE PAIRS ALREADY DETERMINED)
		cal_corners_trim <- cal_corners[, , aspect_non_na >= min_views, ]

	}else{

		# sample.est IS USED FOR DE-BUGGING
		# OF THE PROVIDED SAMPLE ASPECTS FIND WHICH ARE NOT NA IN ANY VIEW
		sample_est_non_na <- sample.est[!sample.est %in% names(aspect_non_na)[which(!aspect_non_na >= min_views)]]

		# CHECK THAT NONE OF DEFINED ASPECTS ARE NA
		if(length(sample_est_non_na) < length(sample.est))
			stop(paste0("The following aspects were not found in at least one view and will not be included in the calibration estimation: ", paste(sample.est[!sample.est %in% sample_est_non_na], collapse=',')))

		# USE DEFINED ASPECTS
		cal_corners_trim <- cal_corners[, , sample_est_non_na, ]
	}

	# CHECK THAT CORNERS WERE FOUND IN AT LEAST 3 ASPECTS
	if(cal_views_found_num < 2) stop(paste0("Corners were only detected in 1 or fewer aspects. At least 2 aspects are required for calibration."))
	
	# PRINT WARNING FOR 3 ASPECTS
	if(cal_views_found_num <= 3){
		if(estimate_cal_coeffs){
			response <- readline(prompt=paste0("\n\t\tCorners were only detected in 3 or fewer aspects. More than 3 aspects are generally required for an accurate calibration. Do you still wish to continue with the calibration? (y/n) : "))
			if(tolower(response) %in% c('no', 'n')) return(1)
		}else{
			warning(paste0("Corners were only detected in 3 or fewer aspects. More than 3 aspects are generally required for an accurate calibration."))
		}
	}

	# SET DEFAULTS FOR COEFFICIENT ESTIMATION SAMPLING PARAMETERS
	if(num.sample.est == 'auto'){
		if(cal_views_found_num < 15){
			num.sample.est <- cal_views_found_num
		}else if(cal_views_found_num >= 15 && cal_views_found_num < 20){
			num.sample.est <- 10
		}else if(cal_views_found_num >= 20 && cal_views_found_num < 25){
			num.sample.est <- 15
		}else if(cal_views_found_num >= 25){
			num.sample.est <- 20
		}
	}else{

		# IF NUMBER OF FOUND ASPECTS IS LESS THAN NUMBER OF ASPECTS TO USE FOR ESTIMATION, MAKE SAMPLE NUMBER OF FOUND ASPECTS
		if(cal_views_found_num < num.sample.est) num.sample.est <- cal_views_found_num
	}

	if(num.sample.sets == 'auto'){
		if(cal_views_found_num <= 5){
			num.sample.sets <- 1
		}else if(cal_views_found_num >= 6 && cal_views_found_num <= 7){
			num.sample.sets <- 2
		}else if (cal_views_found_num >= 8){
			num.sample.sets <- 3
		}
	}		

	if(num.aspects.sample == 'auto'){
		if(cal_views_found_num <= 5){
			num.aspects.sample <- cal_views_found_num
		}else if(cal_views_found_num >= 6 && cal_views_found_num <= 7){
			num.aspects.sample <- 5
		}else if (cal_views_found_num >= 8){
			num.aspects.sample <- 6
		}
	}else{

		# CHECK THAT NUMBER OF ASPECTS TO SAMPLE DOES NOT EXCEED TOTAL NUMBER OF FOUND ASPECTS
		if(num.aspects.sample > num.sample.est) stop(paste0("'num.aspects.sample' (", num.aspects.sample, ") must be less than 'num.sample.est' (", num.sample.est, ")."))
	}

	# GET ESTIMATION SUBSAMPLE INDICES
	sample_trim_est <- floor(seq(from=1, to=dim(cal_corners_trim)[3], length=num.sample.est))
	
	# GET ESTIMATION SUBSAMPLE
	cal_corners_trim_est <- cal_corners_trim[, , sample_trim_est, ]

	# GET OPTIMIZATION SUBSAMPLE - IF AT LEAST 10 ASPECTS NOT USED IN CALIBRATION, USE THOSE FOR OPTIMIZATION
	cal_corners_trim_optim <- NULL
	optim_with_cal <- TRUE
	if(cal_views_found_num - num.sample.est >= 10){

		# INDICES
		sample_trim_optim <- (1:dim(cal_corners_trim)[3])[!(1:dim(cal_corners_trim)[3] %in% sample_trim_est)]
		
		# CAP AT MAX
		sample_trim_optim <- sample_trim_optim[floor(seq(from=1, to=length(sample_trim_optim), length=min(length(sample_trim_optim), max.sample.optim)))]

		# SUBSAMPLE
		cal_corners_trim_optim <- cal_corners_trim[, , sample_trim_optim, ]
		
		optim_with_cal <- FALSE
	}
	
	# IF MORE THAN TWO VIEWS, CHECK THAT ONE VIEW ISNT STRANDED WITHOUT COMMON CORNERS WITH OTHER VIEWS
	if(num_views > 2){
		common_by_view <- setNames(rep(NA, num_views), img_sub_dir)
		for(view in 1:num_views){
			common_by_view[view] <- sum((non_na_by_view[, view] > 0)*(rowSums(non_na_by_view[, -view]) >= min_views-1))
			if(common_by_view[view] < round((num.aspects.sample+max.sample.optim) / 2)){
				stop(paste0("There are only ", common_by_view[view], " aspects in view '", img_sub_dir[view], "' for which corners were detected in at least ", min_views-1, " other view(s). Try increasing 'num.aspects.read'."))
			}
		}
	}

	# CREATE ACCURACY CHECK FOLDERS
	if(num.sample.est != 'auto'){

		# CREATE ERROR TESTS FOLDER
		if(!file.exists(paste0(calib_dir, 'Error tests'))) dir.create(paste0(calib_dir, 'Error tests'))

		# FIND COMBINATIONS OF ALL VIEWS
		view_combo_sub_dir <- NULL
		if(num_views == 2){

			view_combos <- list(c(1,2))

		}else if(num_views == 3){

			view_combos <- list(c(1,2), c(2,3), c(1,3), c(1,2,3))

			# CREATE SUB-FOLDER NAMES
			view_combo_sub_dir <- rep(NA, length(view_combos))
			for(i in 1:length(view_combos)) view_combo_sub_dir[i] <- paste(img_sub_dir[view_combos[[i]]], collapse=', ')

			# CREATE SUB-FOLDERS IF THEY DO NOT EXIST
			for(dir_name in view_combo_sub_dir) if(!file.exists(paste0(calib_dir, 'Error tests/', dir_name))) dir.create(paste0(calib_dir, 'Error tests/', dir_name))
		}
	}

	if(num.sample.sets == 1){

		set.seed(42);
		cal_sample_sets <- list(sort(sample(1:num.sample.est, num.aspects.sample)))

	}else{

		# CREATE RANDOM SAMPLE SETS
		cal_sample_sets <- list()
		for(i in 1:max(num.sample.sets, 40)){

			# GET SAMPLE SET
			set.seed(i);
			sample_set <- sample(1:num.sample.est, num.aspects.sample);

			if(num_views > 2){

				# FIND NUMBER OF NON-NA VIEWS FOR EACH ASPECT AND VIEW
				non_na_by_view <- apply(!is.na(cal_corners_trim_est[, , sample_set, ]), c(3, 4), 'sum') > 0

				common_by_view <- setNames(rep(NA, num_views), img_sub_dir)
				for(view in 1:num_views) common_by_view[view] <- sum((non_na_by_view[, view] > 0)*(rowSums(non_na_by_view[, -view]) >= min_views-1))

				# CHECK THAT NO VIEW IS UNDERREPRESENTED
				if(sum(common_by_view < floor(num.aspects.sample / 2)) > 0) next
			}
			
			# MAKE SURE THAT SAMPLE IS NOT THE SAME AS ANY OTHER SET
			if(length(cal_sample_sets) > 0){
				go_next <- FALSE
				for(j in 1:length(cal_sample_sets)) if(sum(!sample_set %in% cal_sample_sets[[j]]) == 0){go_next <- TRUE;break}
				if(go_next) next
			}

			# SAVE SAMPLE SET
			cal_sample_sets[[length(cal_sample_sets)+1]] <- sort(sample_set)
			
			# IF NUMBER OF SAMPLE SETS HAS BEEN REACHED, BREAK
			if(length(cal_sample_sets) == num.sample.sets) break
		}
	}

	if(print.progress){
		cat(paste0("\t\tNumber of aspects from total that will be sampled for calibration coefficient estimation (num.sample.est): ", num.sample.est, "\n"))
		cat(paste0("\t\tNumber of unique sets of aspects to try (num.sample.sets): ", num.sample.sets, "\n"))
		cat(paste0("\t\tNumber of aspects to sample for each set (num.aspects.sample): ", num.aspects.sample, "\n"))
		if(num.sample.sets > 1){
			cat(paste0("\t\tAspects in each sets:\n"))
			for(i in 1:length(cal_sample_sets)) cat(paste0("\t\t\t", i, ": ", paste(dimnames(cal_corners_trim_est)[[3]][cal_sample_sets[[i]]], collapse=", "), "\n"))
		}
		cat(paste0("\t\tUse calibration aspects to determine the best calibration set?: ", optim_with_cal, "\n"))
		if(!optim_with_cal) cat(paste0("\t\tNumber of aspects in optimization set: ", dim(cal_corners_trim_optim)[3], "\n"))
	}
	
	# ESTIMATE DLT CALIBRATION COEFFICIENTS
	if(estimate_cal_coeffs){
		
		# SAVE ALL 2D COORDINATES FOR ESTIMATION
		coor_2d <- cal_corners_trim_est

		# REDUCE CORNER NUMBER FOR ALL GRIDS BY FITTING PERSPECTIVE GRID
		# SET REDUCED GRID DIMENSIONS
		rx <- ry <- 3
		
		# SET REDUCED GRID SIZES
		sx <- ((nx-1)*sq.size.num) / (rx-1)
		sy <- ((ny-1)*sq.size.num) / (ry-1)

		# EMPTY REDUCED GRID DIMENSION ARRAY
		coor_2d_red <- array(NA, dim=c(rx*ry, 2, dim(coor_2d)[3], dim(coor_2d)[4]), dimnames=list(NULL, c('x', 'y'), dimnames(coor_2d)[[3]], dimnames(coor_2d)[[4]]))

		if(print.progress) cat('\n\t\tReduce number of corners (', dim(coor_2d)[3]*dim(coor_2d)[4], ' checkerboards total)\n', sep='')

		for(i in 1:dim(coor_2d)[3]){
			for(j in 1:dim(coor_2d)[4]){

				if(print.progress) cat('\t\t\t', (i-1)*dim(coor_2d)[4] + j, ') Aspect: ', dimnames(coor_2d)[[3]][i], '; View: ', dimnames(coor_2d)[[4]][j], '; ', sep='')

				# DOWNSAMPLE THE NUMBER OF CORNERS
				coor_2d_red[, , i, j] <- resampleGridImagePoints(pts=coor_2d[, , i, j], 
					nx=nx, rx=rx, ry=ry, fit.min.break=fit.min.break, print.progress=print.progress)$pts
				
				# fit.min.break is mean error
			}
		}

		# EMPTY LIST TO STORE CALIBRATION RESULTS
		dlt_cal_cam_list <- list()

		# CALIBRATE CAMERAS WITH EACH SAMPLE SET
		cal_optim <- rep(NA, num.sample.sets)
		for(i in 1:num.sample.sets){

			if(print.progress) cat('\n\t\tAspects used in estimating coefficients: ', paste(dimnames(coor_2d_red)[[3]][cal_sample_sets[[i]]], collapse=", "), '\n', sep='')

			# ESTIMATE DLT CALIBRATION COEFFICIENTS FOR ALL VIEWS
			dlt_cal_cam <- dltCalibrateCameras(coor.2d=coor_2d_red[, , cal_sample_sets[[i]], ], nx=3, 
				grid.size=c(sx, sy), print.progress=print.progress, print.tab='\t\t', 
				reduce.grid.dim=FALSE, objective.min=objective.min, objective.min.break=objective.min.break, 
				min.views=min_views, grid.incl.min=2, nlm.calls.max=nlm.calls.max, ...)
			
			if(print.progress) cat('\n')

			dlt_cal_cam_list[[i]] <- list()

			# SET COLUMN NAMES FOR COEFFICIENTS
			colnames(dlt_cal_cam$cal.coeff) <- img_sub_dir

			# ADD COEFFICIENT ESTIMATION RESULTS
			dlt_cal_cam_list[[i]][['cal.set.num']] <- i
			dlt_cal_cam_list[[i]][['cal.coeff']] <- dlt_cal_cam$cal.coeff
			dlt_cal_cam_list[[i]][['mean.reconstruct.rmse']] <- dlt_cal_cam$mean.reconstruct.rmse
			dlt_cal_cam_list[[i]][['coefficient.rmse']] <- dlt_cal_cam$coefficient.rmse
			dlt_cal_cam_list[[i]][['cal.sample.aspects']] <- dimnames(coor_2d_red)[[3]][cal_sample_sets[[i]]]

			if(!is.null(cal_corners_trim_optim)){

				# TEST CALIBRATION ACCURACY AGAINST OPTIM SET
				test_calibration <- dltTestCalibration(dlt_cal_cam$cal.coeff, cal_corners_trim_optim, nx, sq.size.num)
			
				# SAVE IPD ERROR FOR CHOOSING BEST CALIBRATION
				cal_optim[i] <- mean(test_calibration$ipd.rmse)

				# SAVE OTHER CALIBRATION ACCURACY TEST RESULTS
				dlt_cal_cam_list[[i]][['cal.ipd.rmse']] <- test_calibration$ipd.rmse
				dlt_cal_cam_list[[i]][['cal.ipd.abs.mean']] <- mean(abs(test_calibration$ipd.error))
				dlt_cal_cam_list[[i]][['cal.ipd.abs.max']] <- max(abs(test_calibration$ipd.error))
				dlt_cal_cam_list[[i]][['cal.ipd.abs.sd']] <- sd(abs(test_calibration$ipd.error))
				dlt_cal_cam_list[[i]][['cal.ipd.mean']] <- mean(test_calibration$ipd.error)
				dlt_cal_cam_list[[i]][['cal.ee.rms']] <- test_calibration$epipolar.rmse
				dlt_cal_cam_list[[i]][['cal.ee.mean']] <- mean(test_calibration$epipolar.error)
				dlt_cal_cam_list[[i]][['cal.ee.max']] <- max(test_calibration$epipolar.error)
				dlt_cal_cam_list[[i]][['cal.ee.sd']] <- sd(test_calibration$epipolar.error)

			}else{
		
				# SAVE RECONSTRUCT ERROR FOR CHOOSING BEST CALIBRATION
				cal_optim[i] <- dlt_cal_cam$mean.reconstruct.rmse
			}
		}

		# GET SAMPLE SET WITH MINIMUM ERROR
		min_error_set <- which.min(cal_optim)
	
		# COPY ELEMENTS FROM MINIMUM ERROR SET
		for(ename in names(dlt_cal_cam_list[[min_error_set]])) cal.list[[ename]] <- dlt_cal_cam_list[[min_error_set]][[ename]]

		# SAVE CALIBRATION COEFFICIENTS
		cal_coeff <- dlt_cal_cam_list[[min_error_set]]$cal.coeff
		
		# SAVE RESULTS OF COEFFICIENT ESTIMATION
		cal_set_num <- dlt_cal_cam_list[[min_error_set]]$cal.set.num
		coefficient_rmse <- dlt_cal_cam_list[[min_error_set]]$coefficient.rmse
		mean_reconstruct_rmse <- dlt_cal_cam_list[[min_error_set]]$mean.reconstruct.rmse
		cal_sample_aspects <- dlt_cal_cam_list[[min_error_set]]$cal.sample.aspects

		# SAVE CALIBRATION LIST
		list2XML4R(list('calibration' = cal.list), file=cal.file)

		if(print.progress) cat("\n")
	}

	if(print.progress){

		cat('\n')
		if(num.sample.sets > 1){
			cat(paste0("\t\tSelected aspect set (lowest reconstruction error): ", cal_set_num, "\n"))
			cat(paste0("\t\t\tAspects in set: ", paste(dimnames(cal_corners_trim_est)[[3]][cal_sample_sets[[cal_set_num]]], collapse=", "), "\n"))
		}
		cat(paste0("\t\tMean reconstruction RMS Error: ", round(mean_reconstruct_rmse, 4), " px\n"))
		cat(paste0("\t\tDLT Coefficient RMS Error:\n"))
		for(i in 1:length(coefficient_rmse)){
			cat(paste0("\t\t\t", img_sub_dir[i], paste(rep(' ', img_sub_dir_salign[i]), collapse=''), 
				" : ", paste0(round(coefficient_rmse[i], 4), collapse=', '), " px\n"))
		}
	}

	# DEFAULT NULL
	cal_corners_test <- NULL
	if(is.null(cal_corners_trim_optim)){

		if(num_views == 2){

			# TEST ACCURACY USING ALL ASPECTS
			cal_corners_test <- cal_corners_trim

		}else{

			# TEST ACCURACY USING ALL ASPECTS
			#cal_corners_test <- cal_corners_trim
			cal_corners_test <- cal_corners[, , aspect_non_na >= 2, ]
		}

		if(print.progress){
			cat(paste0("\n\t\tTesting accuracy using all aspects (including those used in calibration):\n"))
			cat(paste0("\t\t\t", paste(dimnames(cal_corners_test)[[3]], collapse=", "), "\n"))
		}

	}else{
	
		# SET OPTIM SET FOR TESTING ACCURACY
		cal_corners_test <- cal_corners_trim_optim

		if(print.progress) cat("\n\t\tTesting accuracy using optimization set\n")
		cat(paste0("\t\t\t", paste(dimnames(cal_corners_test)[[3]], collapse=", "), "\n"))
	}

	dlt_test <- NULL
	if(!is.null(cal_corners_test)){

		for(i in 1:length(view_combos)){

			# SET SAVE AS FOLDER
			if(length(view_combos) == 1){
				save_to <- paste0(calib_dir, 'Error tests')
			}else{
				save_to <- paste0(calib_dir, 'Error tests/', view_combo_sub_dir[i])
			}

			# FIND NON-NA ASPECTS AMONG VIEWS
			aspects_non_na <- rowSums(is.na(cal_corners_test[1,1, , view_combos[[i]]])) == 0

			# CREATE ERROR PLOTS
			dlt_test <- createErrorPlots(cal.coeff=cal_coeff[, view_combos[[i]]], 
				corners=cal_corners_test[, , aspects_non_na, view_combos[[i]]], nx=nx, 
				sq.size.num=sq.size.num, sq.size.units=sq.size.units, 
				file=save_to)

			if(print.progress) if(i == length(view_combos)){print(summary(dlt_test, print.tab='\t\t'))}
		}
	}
	
	if(print.progress) cat('\n')

	rlist <- list(
		cal.coeff = cal_coeff, 
		mean.reconstruct.rmse = cal.list$mean.reconstruct.rmse, 
		coefficient.rmse = cal.list$coefficient.rmse
	)

	if(!is.null(dlt_test)) rlist$test.ipd.error = dlt_test$ipd.error

	class(rlist) <- 'calibrateCameras'

	rlist
}

print.calibrateCameras <- function(x, ...){

	class(x) <- ''

	cat('$cal.coeff\n')
	print(x$cal.coeff)
	cat('\n')

	cat('$mean.reconstruct.rmse\n')
	print(x$mean.reconstruct.rmse)
	cat('\n')

	cat('$coefficient.rmse\n')
	print(x$coefficient.rmse)
	cat('\n')

	if(!is.null(x$test.ipd.error)){
		cat('$test.ipd.error\n')
		cat('[1] ', paste(format(x$test.ipd.error[1:min(length(x$test.ipd.error), 5)]), collapse=' '), ' ... and ', length(x$test.ipd.error)-min(length(x$test.ipd.error), 5), ' more\n', sep='')
		cat('\n')
	}
}
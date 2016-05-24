shinyServer(function(input, output) {

	require(jpeg)
	require(tiff)
	require(png)
	
	output$text_output <- renderText({

		# Set print.progress
		print.progress <- FALSE

		#
		json_list <- fromJSON(input$text_input)

		#cat('-------------------------------------------------\n')
		#cat('-----------------------', json_list$submit_ct, '-----------------------\n')
		#cat('-------------------------------------------------\n')
		#print(json_list)

		# Ignore input copied from digitizeImages to server
		if(is.null(json_list$fromBrowser)) return(NULL)

		# Detect exit
		if(!is.null(json_list$exit)) stopApp()

		# Set output list - initialize with inputs
		out <- json_list
		
		# Start update status string
		out$update_status <- ''

		# Default NA
		scaling_wpp <- NA
		ruler_interval_pixels <- NA

		# Change image
		if(!is.null(json_list$change_image)){

			# Get image full file path
			image_full_fpath <- paste0(json_list$prev_wd, '/', json_list$image_fpath)
			
			# Remove images currently in img folder
			if(length(list.files('www/img/')) > 0) file.remove(paste0('www/img/', list.files('www/img/')))

			# Insert view number into image for stereo images - makes sure that image name differs for different views to clear the image cache
			if(!is.null(json_list$add_view)){
				str <- json_list$image_fname
				regexpr_r <- regexpr(pattern='[.][A-Za-z]*$', str)
				json_list$image_fname <- paste0(substr(str,1, regexpr_r[1]-1), '_view_', json_list$add_view, substr(str,regexpr_r[1],nchar(str)))
			}

			# Copy image into www/img folder
			if(print.progress) cat(paste0('Submit ', json_list$submit_ct, ': Changing to image "', json_list$image_fpath, '"\n'))
			file.copy(from=image_full_fpath, to=paste0('www/img/', json_list$image_fname))
			
			# Get image size
			out$image_size <- file.info(image_full_fpath)$size
			out$image_fname <- json_list$image_fname
			
			# Get image dimensions
			if(grepl(pattern='[.]jpg$|[.]jpeg$', x=image_full_fpath, ignore.case=TRUE)) img_dim <- dim(readJPEG(image_full_fpath, native=TRUE))
			if(grepl(pattern='[.]png$', x=image_full_fpath, ignore.case=TRUE)) img_dim <- dim(readPNG(image_full_fpath, native=TRUE))
			if(grepl(pattern='[.]tif$|[.]tiff$', x=image_full_fpath, ignore.case=TRUE)) img_dim <- dim(readTIFF(image_full_fpath, native=TRUE))
		
			# Set image dimensions
			out$image_width <- img_dim[2]
			out$image_height <- img_dim[1]
			
			# Set initial shape lists
			out$landmarks <- list()
			out$curves <- list()
			out$rulers <- list()
			out$corners <- list()
			
			# Check for shape filepath
			shapes <- list()
			if(!is.null(json_list$shape_fpath)){
				
				# Check if file exists
				full_shape_fpath <- paste0(json_list$prev_wd, '/', json_list$shape_fpath)

				if(file.exists(full_shape_fpath)){

					# Add shapes to shapes list
					shapes_list <- XML4R2list(full_shape_fpath)$shapes

					if('landmarks.pixel' %in% names(shapes_list)) shapes$landmarks <- shapes_list$landmarks.pixel
					if('curves.control' %in% names(shapes_list)){
						shapes$curves <- shapes_list$curves.control
						if(!is.null(shapes$curves) && length(shapes$curves) > 0){
							for(i in 1:length(shapes$curves)) out$curves[[i]] <- c(names(shapes$curves)[i], t(shapes$curves[[names(shapes$curves)[i]]]))
						}
					}
					if('ruler.points' %in% names(shapes_list)){
						shapes$rulers <- shapes_list$ruler.points
						if(!is.null(shapes$rulers) && nrow(shapes$rulers) > 0){
							for(i in 1:nrow(shapes$rulers)) out$rulers[[i]] <- c(rownames(shapes$rulers)[i], shapes$rulers[i, ])
						}
					}
					if('checker.pixel' %in% names(shapes_list)){
						shapes$corners <- shapes_list$checker.pixel
						if(!is.null(shapes$corners) && nrow(shapes$corners) > 0){
							for(i in 1:nrow(shapes$corners)) out$corners[[i]] <- shapes$corners[i, ]
						}
					}

					# Add metadata from shape file
					if('ruler.pixel' %in% names(shapes_list)) out$ruler_pixel <- shapes_list$ruler.pixel
					if('ruler.interval' %in% names(shapes_list)) out$ruler_interval <- shapes_list$ruler.interval
					if('checkerboard.nx' %in% names(shapes_list)) out$checkerboard_nx <- shapes_list$checkerboard.nx
					if('checkerboard.ny' %in% names(shapes_list)) out$checkerboard_ny <- shapes_list$checkerboard.ny
					if('square.pixel' %in% names(shapes_list)) out$checker_square_pixel <- shapes_list$square.pixel
					if('square.size' %in% names(shapes_list)) out$checker_square_world <- shapes_list$square.size
					if('scaling.units' %in% names(shapes_list)) out$scaling_units <- shapes_list$scaling.units
					if('scaling' %in% names(shapes_list)) out$scaling <- shapes_list$scaling
				}
			}

			# If no landmarks already added and landmarks filepath is specified, read in landmarks
			if(is.null(shapes$landmarks) && !is.null(json_list$landmark_fpath)){
				full_shape_fpath <- paste0(json_list$prev_wd, '/', json_list$landmark_fpath)
				if(file.exists(full_shape_fpath) && file.info(full_shape_fpath)$size > 1) shapes$landmarks <- as.matrix(read.table(full_shape_fpath, row.names=1, sep="\t"))
			}

			# If no curves already added and curve control point filepath is specified, read in curves
			if(is.null(shapes$curves) && !is.null(json_list$control_points_fpath)){
				full_shape_fpath <- paste0(json_list$prev_wd, '/', json_list$control_points_fpath)
				if(file.exists(full_shape_fpath) && file.info(full_shape_fpath)$size > 1){
					control_points <- readBezierControlPoints(full_shape_fpath)
					if(length(control_points) > 0){
						for(i in 1:length(control_points))
							out$curves[[i]] <- c(names(control_points)[i], c(t(control_points[[names(control_points)[i]]][[1]])))
					}
				}
			}

			# Format shape data for proper interpretation from JSON string
			if(!is.null(shapes$landmarks) && nrow(shapes$landmarks) > 0)
				for(i in 1:nrow(shapes$landmarks)) out$landmarks[[i]] <- c(rownames(shapes$landmarks)[i], shapes$landmarks[i, ])

			if(!is.null(out$landmarks) || !is.null(out$curves) || !is.null(out$rulers) || !is.null(out$corners)) out$load_shapes_from_file <- 1

			#print(out)
		}

		# Find ruler interval
		if(!is.null(json_list$ruler_points) && (!is.null(json_list$save_shapes) || !is.null(json_list$find_ruler_interval))){

			if(print.progress) cat(paste0('Estimating the ruler interval\n'))

			# READ RULER POINTS TO MATRIX
			ruler_points <- matrix(cbind(sapply(json_list$ruler_points, "[[", 1), sapply(json_list$ruler_points, "[[", 2)), 
				nrow=length(json_list$ruler_points), ncol=2)

			# CONVERT "-" TO NA
			ruler_points[ruler_points == '-'] <- NA
			
			# MAKE NUMERIC
			ruler_points <- matrix(as.numeric(ruler_points), nrow(ruler_points), ncol(ruler_points), dimnames=list(paste0("Ruler point ", 1:nrow(ruler_points)), NULL))

			# AT LEAST TWO POINTS
			if(sum(!is.na(ruler_points[, 1])) >= 2){

				# FIND CONSECUTIVE STRING OF RULER POINTS
				start <- FALSE
				ruler_points_c <- matrix(NA, 0, ncol(ruler_points))
				for(i in 1:nrow(ruler_points)){
					if(is.na(ruler_points[i, 1]) && start == TRUE) break
					if(is.na(ruler_points[i, 1]) && start == FALSE) next
					if(!is.na(ruler_points[i, 1])) start <- TRUE
				
					ruler_points_c <- rbind(ruler_points_c, ruler_points[i, ])
				}
				
				ruler_points_c <- ruler_points_c
				
				if(nrow(ruler_points_c) < 2){
					#print('too few points')
				}else if(nrow(ruler_points_c) == 2){
					
					# FIND DISTANCE BETWEEN TWO POINTS
					ruler_interval_pixels <- sqrt(sum((ruler_points_c[2, ] - ruler_points_c[1, ])^2))

				}else if(nrow(ruler_points_c) > 2){

					# FIT A REGULARLY SPACED POINTS MODEL TO EACH DIMENSION OF THE POINTS MATRIX
					fit_x <- nlminb(start=c(ruler_points_c[1, 1], ruler_points_c[2, 1]-ruler_points_c[1, 1]),
							objective=gridPointsFitError, nx=nrow(ruler_points_c), points=ruler_points_c[, 1])
					fit_y <- nlminb(start=c(ruler_points_c[1, 2], ruler_points_c[2, 2]-ruler_points_c[1, 2]),
						objective=gridPointsFitError, nx=nrow(ruler_points_c), points=ruler_points_c[, 2])

					# FIND DISTANCE BETWEEN TWO POINTS
					ruler_interval_pixels <- sqrt(fit_x$par[2]^2 + fit_y$par[2]^2)
				}
			}

			if(!is.null(json_list$find_ruler_interval)){
				out$ruler_pixel <- ruler_interval_pixels
				out$load_new_ruler_interval <- 1
			}
		}		

		# Find checkerboard corners
		if(!is.null(json_list$find_checkerboard_corners)){
			
			if(print.progress) cat(paste0('Finding checkerboard corners\n'))

			# GET PATH TO IMAGE FILE
			img_fpath <- paste0(getwd(), '/www/img/', json_list$image_fname)

			# FIND CHECKERBOARD CORNERS
			corners <- findCheckerboardCorners(img_fpath, nx=json_list$checkerboard_nx, ny=json_list$checkerboard_ny)

			#shapes_list <- XML4R2list('/Users/aaron/Documents/Research/R Package Tests/StereoMorph/Digitizing App/aStereo/Shapes/mug_001.txt')$shapes
			#corners <- shapes_list$checker.pixel
			#img_name_txt <- gsub('.jpg$|.jpeg$|.tif$|.tiff$', '.txt', json_list$image_name, ignore.case=TRUE)
			#corners <- as.matrix(read.table(file=paste0('/Users/aaron/Documents/Research/R Package Tests/StereoMorph/Digitizing App/Corners/', img_name_txt), sep="\t"))

			# IF CORNERS WERE FOUND, FIND MEAN SQUARE SIZE
			if(!is.na(corners[1, 1])){

				# FIT MINIMAL PARAMETER CHECKERBOARD MODEL TO POINTS
				nlminb_fit_x <- nlminb(start=c(corners[1, 1], corners[2, 1]-corners[1, 1], 0), 
					objective=gridPointsFitError, points=corners[, 1], nx=json_list$checkerboard_nx, ny=json_list$checkerboard_ny)
				nlminb_fit_y <- nlminb(start=c(corners[1, 2], corners[2, 2]-corners[1, 2], 0), 
					objective=gridPointsFitError, points=corners[, 2], nx=json_list$checkerboard_nx, ny=json_list$checkerboard_ny)

				# GET SQUARE SIZE FROM MODEL FIT
				out$checker_square_pixel <- sqrt(nlminb_fit_x$par[2]^2 + nlminb_fit_y$par[2]^2)

				out$update_status <- paste0(out$update_status, json_list$checkerboard_nx*json_list$checkerboard_ny, " corners found successfully.")

			}else{
				out$update_status <- paste0(out$update_status, "Checkerboard corner detection unsuccessful.")
			}

			out$corners <- list()
			for(i in 1:nrow(corners)) out$corners[[i]] <- corners[i, ]

			out$load_new_corners <- 1
		}
		
		# Save shapes
		if(!is.null(json_list$save_shapes)){

			if(print.progress) cat(paste0('Saving shape data...\n'))

			shapes <- list()
			curve_points <- matrix(NA, nrow=0, ncol=2)

			# $save_as_shapes
			# $save_as_landmarks
			# $save_as_control_points
			# $save_as_curve_points

			# Check if scaling is available from checkerboard corners
			if(!is.null(json_list$checker_square_world) && !is.null(json_list$checker_square_pixel))
				scaling_wpp <- json_list$checker_square_world / json_list$checker_square_pixel

			# Check if scaling is available from ruler points
			if(!is.null(json_list$ruler_interval_world) && !is.null(json_list$ruler_interval_pixels))
				scaling_wpp <- json_list$ruler_interval_world / json_list$ruler_interval_pixels
			
			# Check if scaling is available from just submitted ruler points
			if(!is.null(json_list$ruler_interval_world) && !is.na(ruler_interval_pixels))
				scaling_wpp <- json_list$ruler_interval_world / ruler_interval_pixels

			# If there are no original numbers from which to calculate scaling but it is saved in output, use that scaling
			if(is.na(scaling_wpp) && !is.null(json_list$save_as_shapes) && !is.null(json_list$scaling) && json_list$scaling != "NA") scaling_wpp <- json_list$scaling

			# Add metadata for shapes file
			if(!is.null(json_list$save_as_shapes)){
				if(!is.null(json_list$checker_square_world) && json_list$checker_square_world != "NA") shapes$square.size <- json_list$checker_square_world
				if(!is.null(json_list$checker_square_pixel) && json_list$checker_square_pixel != "NA") shapes$square.pixel <- json_list$checker_square_pixel
				if(!is.null(json_list$checkerboard_ny) && json_list$checkerboard_ny != "NA") shapes$checkerboard.ny <- json_list$checkerboard_ny
				if(!is.null(json_list$checkerboard_nx) && json_list$checkerboard_nx != "NA") shapes$checkerboard.nx <- json_list$checkerboard_nx
				if(!is.null(json_list$ruler_interval_pixels) && json_list$ruler_interval_pixels != "NA") shapes$ruler.pixel <- json_list$ruler_interval_pixels
				if(!is.null(json_list$ruler_interval_world) && json_list$ruler_interval_world != "NA") shapes$ruler.interval <- json_list$ruler_interval_world
				if(!is.null(json_list$scaling_units) && json_list$scaling_units != "NA") shapes$scaling.units <- json_list$scaling_units
				shapes$scaling <- scaling_wpp
				shapes$image.name <- json_list$image_fname
			}
			
			# Parse landmarks for saving
			if(!is.null(json_list$landmarks) && (!is.null(json_list$save_as_shapes) || !is.null(json_list$save_as_landmarks))){

				# Convert landmarks into matrix
				landmarks_pixel <- matrix(cbind(as.numeric(sapply(json_list$landmarks, "[[", 2)), as.numeric(sapply(json_list$landmarks, "[[", 3))), 
					nrow=length(json_list$landmarks), ncol=2,
					dimnames=list(sapply(json_list$landmarks, "[[", 1), NULL))
				
				if(nrow(landmarks_pixel) > 0){

					shapes$landmarks.pixel <- landmarks_pixel

					# Get scaled landmarks
					if(!is.na(scaling_wpp)){

						# Scale landmarks
						shapes$landmarks.scaled <- shapes$landmarks.pixel*scaling_wpp
						
						# Flip about x-axis
						shapes$landmarks.scaled <- cbind(shapes$landmarks.scaled[, 1], -shapes$landmarks.scaled[, 2])
					}

					if(nrow(shapes$landmarks.pixel) == 1){plural <- ''}else{plural <- 's'}
					out$update_status <- paste0(out$update_status, nrow(shapes$landmarks.pixel), " landmark", plural, " saved. ")
				}
			}

			# Parse control points for saving
			if(!is.null(json_list$control_points)){

				curve_save_status <- '(control points only)'
				curves_saved <- FALSE

				# Convert control points to string for saving in control point file (old version)
				curve_string <- c()
				if(length(json_list$control_points) > 0 && !is.null(json_list$save_as_control_points)){
					for(i in 1:length(json_list$control_points)){
						curve_string <- paste0(curve_string, paste(json_list$control_points[[i]], collapse='\t'))
						if(i < length(json_list$control_points)) curve_string <- paste0(curve_string, '\n')
					}
					curves_saved <- TRUE
				}

				# Convert control points to list of matrices for saving to shapes file and/or making curve points
				if(length(json_list$control_points) > 0 && (!is.null(json_list$save_as_shapes) || !is.null(json_list$save_as_curve_points))){

					shapes$curves.control <- list()
					curve_names <- rep(NA, length(json_list$control_points))
					for(i in 1:length(json_list$control_points)){
						shapes$curves.control[[i]] <- matrix(unlist(json_list$control_points[[i]][2:length(json_list$control_points[[i]])]), ncol=2, byrow=TRUE)
						curve_names[i] <- json_list$control_points[[i]][1]
					}
					names(shapes$curves.control) <- curve_names
					curves_saved <- TRUE
				}

				# Get curve points as single matrix and for saving as list of matrices
				if(length(json_list$control_points) > 0 && (!is.null(json_list$save_as_shapes) || !is.null(json_list$save_as_curve_points)) && json_list$save_curve_points){

					# Get curve points as single matrix
					shapes$curves.pixel <- list()
					num_curve_sets <- 0
					for(i in 1:length(shapes$curves.control)){

						if(nrow(shapes$curves.control[[i]]) <= 2) next

						points_on_bezier <- pointsOnBezier(p=shapes$curves.control[[i]], method='adjoining', deg=2)

						# CIRCUMVENT BUG IN POINTSONBEZIER WHERE LAST POINT OVERSHOOTS BY ONE
						# IF SECOND TO LAST POINT IS THE SAME AS THE LAST POINT OF M, GET RID OF LAST POINT
						if(sum(abs(points_on_bezier$points[nrow(points_on_bezier$points)-1, ] - shapes$curves.control[[i]][nrow(shapes$curves.control[[i]]), ])) == 0)
							points_on_bezier$points <- points_on_bezier$points[-nrow(points_on_bezier$points), ]

						shapes$curves.pixel[[names(shapes$curves.control)[i]]] <- points_on_bezier$points

						rownames(points_on_bezier$points) <- paste0(
							names(shapes$curves.control)[i], 
							formatC(1:nrow(points_on_bezier$points), width=4, format="d", flag="0"))

						curve_points <- rbind(points_on_bezier$points, curve_points)
						num_curve_sets <- num_curve_sets + 1
					}
					
					# Get scaled curve points
					if(!is.na(scaling_wpp) && length(shapes$curves.pixel) > 0){

						# Start with pixel coordinates
						shapes$curves.scaled <- shapes$curves.pixel

						for(i in 1:length(shapes$curves.scaled)){

							# Scale curve points
							shapes$curves.scaled[[i]] <- shapes$curves.scaled[[i]]*scaling_wpp

							# Flip about x-axis
							shapes$curves.scaled[[i]] <- cbind(shapes$curves.scaled[[i]][, 1], -shapes$curves.scaled[[i]][, 2])
						}
					}

					curve_save_status <- '(control and curve points)'
					curves_saved <- TRUE
				}

				if(curves_saved){
					if(length(json_list$control_points) == 1){plural <- ''}else{plural <- 's'}
					out$update_status <- paste0(out$update_status, length(json_list$control_points), " curve", plural, " saved ", curve_save_status, ". ")
				}
			}

			if(!is.null(json_list$ruler_points)){

				# REMOVE NA ROWS
				ruler_points <- matrix(ruler_points[!is.na(ruler_points[, 1]), ], ncol=2, dimnames=list(rownames(ruler_points)[!is.na(ruler_points[, 1])], NULL))

				if(nrow(ruler_points) > 0){
					shapes$ruler.points <- ruler_points
					if(nrow(shapes$ruler.points) == 1){plural <- ''}else{plural <- 's'}
					out$update_status <- paste0(out$update_status, nrow(shapes$ruler.points), " ruler point", plural, " saved. ")
				}
			}

			if(!is.null(json_list$corners)){

				# Read corners to matrix
				corners <- matrix(as.numeric(cbind(sapply(json_list$corners, "[[", 1), sapply(json_list$corners, "[[", 2))), 
					nrow=length(json_list$corners), ncol=2)
				
				if(nrow(corners) > 0){
					shapes$checker.pixel <- corners
					if(nrow(shapes$checker.pixel) == 1){plural <- ''}else{plural <- 's'}
					out$update_status <- paste0(out$update_status, nrow(shapes$checker.pixel), " checkerboard corner", plural, " saved. ")
				}
			}

			# Save to shapes file
			if(!is.null(json_list$save_as_shapes)){
				list2XML4R(list=list('shapes'=shapes), file=paste0(json_list$prev_wd, '/', json_list$save_as_shapes))
			}

			# Save to landmarks file
			if(!is.null(json_list$save_as_landmarks)){
				write.table(x=shapes$landmarks.pixel, file=paste0(json_list$prev_wd, '/', json_list$save_as_landmarks),
					quote=FALSE, sep="\t", col.names=FALSE)
			}

			# Save to curve control points file
			if(!is.null(json_list$save_as_control_points)){
				write(curve_string, file=paste0(json_list$prev_wd, '/', json_list$save_as_control_points))
			}

			# Save to curve points file
			if(!is.null(json_list$save_as_curve_points) && nrow(curve_points) > 0){
				write.table(curve_points, file=paste0(json_list$prev_wd, '/', json_list$save_as_curve_points), 
				quote=F, sep="\t", col.names=F, row.names=T)
			}
		}

		# Find epipolar line
		if(!is.null(json_list$find_epipolar_line)){

			slopes <- c()
			intercepts <- c()

			for(i in 1:length(json_list$all_views)){
			
				# Skip current view
				if(i == json_list$current_view) next

				# Get file path
				full_shape_fpath <- paste0(json_list$prev_wd, '/', json_list$all_views[i])

				# Skip non-existant or empty files
				if(!file.exists(full_shape_fpath) || file.info(full_shape_fpath)$size == 0) next
				
				# Get shapes
				shapes_list <- XML4R2list(full_shape_fpath)$shapes

				# Skip if required shape data not found
				if(!json_list$shape_type %in% names(shapes_list)) next

				# Landmark from view other than current
				landmark  <- rep(NA, 2)
				
				if(json_list$shape_type == 'landmarks.pixel'){
				
					# Skip if landmark is not in matrix
					if(!json_list$landmark_name %in% rownames(shapes_list$landmarks.pixel)) next

					# Get landmark
					landmark <- shapes_list$landmarks.pixel[json_list$landmark_name, ]
				}

				# Skip if landmark is NA
				if(is.na(landmark[1])) next
				
				# Read calibration coefficients to matrix
				cal_coeffs <- matrix(as.numeric(unlist(json_list$cal_coeffs)), nrow=11, byrow=TRUE)

				# Find epipolar line
				dlt_epipolar_line <- dltEpipolarLine(landmark, cal_coeffs[, i], cal_coeffs[, json_list$current_view])
				
				slopes <- c(slopes, dlt_epipolar_line$m)
				intercepts <- c(intercepts, dlt_epipolar_line$b)
			}
			
			if(length(slopes) == 0){
				out$slope <- "NA"
				out$intercept <- "NA"
			}else{
				out$slope <- mean(slopes)
				out$intercept <- mean(intercepts)
			}
		}

		out_str <- listToJSONStr(out, direct=TRUE)

		#cat(out_str);cat("\n")

		out_str
	})
})

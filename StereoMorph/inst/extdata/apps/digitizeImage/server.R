# SET CONSTANTS
landmarks_file <- NULL
control_points_file <- NULL
curve_points_file <- NULL
shapes_file <- NULL

if(!is.null(init_params$landmarks_file) && !grepl(init_params$prev_wd, init_params$landmarks_file, fixed=TRUE)){
	landmarks_file <- paste0(init_params$prev_wd, '/', init_params$landmarks_file)
}else{
	landmarks_file <- init_params$landmarks_file
}
if(!is.null(init_params$control_points_file) && !grepl(init_params$prev_wd, init_params$control_points_file, fixed=TRUE)){
	control_points_file <- paste0(init_params$prev_wd, '/', init_params$control_points_file)
}else{
	control_points_file <- init_params$control_points_file
}
if(!is.null(init_params$curve_points_file) && !grepl(init_params$prev_wd, init_params$curve_points_file, fixed=TRUE)){
	curve_points_file <- paste0(init_params$prev_wd, '/', init_params$curve_points_file)
}else{
	curve_points_file <- init_params$curve_points_file
}
if(!is.null(init_params$shapes_file) && !grepl(init_params$prev_wd, init_params$shapes_file, fixed=TRUE)){
	shapes_file <- paste0(init_params$prev_wd, '/', init_params$shapes_file)
}else{
	shapes_file <- init_params$shapes_file
}

shinyServer(function(input, output) {

	output$text_output <- renderText({

		# PARSE JSON
		json_list <- fromJSON(input$text_input)

		# START OUTPUT LIST
		out <- list()
		out$update_status <- ''
		out$submit_ct <- json_list$submit_ct

		# INPUT IS FROM BROWSER
		if(!is.null(json_list$fromBrowser)){

			# SET DEFAULTS
			save_shapes <- FALSE
			apply_scaling <- FALSE

			# START LINES FOR SHAPE FILE
			shapes_file_write <- c()

			# CHECK IF SCALING IS PROVIDED (WILL BE OVERWRITTEN BY RULER POINTS IF AVAILABLE)
			if(!is.null(json_list$checker_square_world) && !is.null(json_list$checker_square_pixel)){
			
				# FIND SCALING FROM PIXELS TO WORLD UNITS
				scaling_wpp <- json_list$checker_square_world / json_list$checker_square_pixel

				if(!is.na(scaling_wpp)) apply_scaling <- TRUE
			}

			# CHECK IF SCALING IS PROVIDED (WILL BE OVERWRITTEN BY RULER POINTS IF AVAILABLE)
			if(!is.null(json_list$ruler_interval_world) && !is.null(json_list$ruler_interval_pixels)){
			
				# FIND SCALING FROM PIXELS TO WORLD UNITS
				scaling_wpp <- json_list$ruler_interval_world / json_list$ruler_interval_pixels

				if(!is.na(scaling_wpp)) apply_scaling <- TRUE
			}
			
			# FIND CHECKERBOARD CORNERS
			if(!is.null(json_list$find_checkerboard_corners)){
				
				# SET DEFAULT
				checker_square_pixel <- "NA"
				
				# GET PATH TO IMAGE FILE
				img_fpath <- paste0(getwd(), '/www/img/', json_list$image_name)

				# FIND CHECKERBOARD CORNERS
				corners <- findCheckerboardCorners(img_fpath, nx=json_list$checkerboard_nx, ny=json_list$checkerboard_ny)
				
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
					checker_square_pixel <- sqrt(nlminb_fit_x$par[2]^2 + nlminb_fit_y$par[2]^2)

					out$update_status <- paste0(out$update_status, json_list$checkerboard_nx*json_list$checkerboard_ny, " corners found successfully.")

				}else{
					out$update_status <- paste0(out$update_status, "Checkerboard corner detection unsuccessful.")
				}
			}

			# GET RULER POINTS
			if(!is.null(json_list$ruler_points) && (!is.null(json_list$submit_shapes) || !is.null(json_list$find_ruler_interval))){

				# DEFAULT NA
				ruler_interval_pixels <- NA

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
				}else{
					#print('too few points')
				}
				
				# CHECK IF RULER INTERVAL IN WORLD UNITS IS PROVIDED
				if(!is.null(json_list$ruler_interval_world) && !is.na(ruler_interval_pixels)){
				
					# FIND SCALING FROM PIXELS TO WORLD UNITS
					scaling_wpp <- json_list$ruler_interval_world / ruler_interval_pixels

					if(!is.na(scaling_wpp)) apply_scaling <- TRUE
				}
			}

			# DETECT LANDMARK SUBMISSION
			if(!is.null(json_list$submit_landmarks)){
				landmarks <- matrix(cbind(as.numeric(sapply(json_list$landmarks, "[[", 2)), as.numeric(sapply(json_list$landmarks, "[[", 3))), 
					nrow=length(json_list$landmarks), ncol=2,
					dimnames=list(sapply(json_list$landmarks, "[[", 1), NULL))

				# SAVE LANDMARKS TO FILE
				if(!is.null(landmarks_file)) write.table(x=landmarks, file=landmarks_file, quote=FALSE, sep="\t", col.names=FALSE)
				if(!is.null(shapes_file)){

					shapes_file_write <- c(shapes_file_write, '\n\t<landmarks.pixel type=matrix rownames=TRUE as.numeric=TRUE>')
					if(nrow(landmarks) > 0){
						for(i in 1:nrow(landmarks)) shapes_file_write <- c(shapes_file_write, paste0('\t\t', paste(c(rownames(landmarks)[i], landmarks[i, ]), collapse='\t')))
					}
					shapes_file_write <- c(shapes_file_write, '\t</landmarks.pixel>')
					
					# SAVE SCALED LANDMARKS
					if(apply_scaling){
					
						# SCALE LANDMARKS
						landmarks_scaled <- landmarks*scaling_wpp
						
						# FLIP ABOUT X-AXIS
						landmarks_scaled <- cbind(landmarks_scaled[, 1], -landmarks_scaled[, 2])

						shapes_file_write <- c(shapes_file_write, '\n\t<landmarks.scaled type=matrix rownames=TRUE as.numeric=TRUE>')
						if(nrow(landmarks) > 0){
							for(i in 1:nrow(landmarks)) shapes_file_write <- c(shapes_file_write, paste0('\t\t', paste(c(rownames(landmarks)[i], landmarks_scaled[i, ]), collapse='\t')))
						}
						shapes_file_write <- c(shapes_file_write, '\t</landmarks.scaled>')
					}

					save_shapes <- TRUE
				}

				# NOTIFICATION IN R CONSOLE
				#cat(paste0(nrow(landmarks), " landmark", plural, " saved to '", landmarks_file, "'\n"))
				
				if(nrow(landmarks) == 1){plural <- ''}else{plural <- 's'}
				out$update_status <- paste0(out$update_status, nrow(landmarks), " landmark", plural, " saved. ")
			}

			# DETECT CURVE SUBMISSION
			if(!is.null(json_list$submit_curves)){
				
				curve_save_status <- '(control points only)'
				
				# SAVE CURVE CONTROL POINTS
				curve_string <- c()
				if(length(json_list$control_points) > 0){
					for(i in 1:length(json_list$control_points)){
						curve_string <- paste0(curve_string, paste(json_list$control_points[[i]], collapse='\t'))
						if(i < length(json_list$control_points)) curve_string <- paste0(curve_string, '\n')
					}
				}

				if(!is.null(control_points_file)) write(curve_string, control_points_file)
				if(!is.null(shapes_file) && length(json_list$control_points) > 0){
					shapes_file_write <- c(shapes_file_write, '\n\t<curves.control type=list names=TRUE as.numeric=TRUE>')
					for(i in 1:length(json_list$control_points)){
						shapes_file_write <- c(shapes_file_write, paste0('\t\t<', json_list$control_points[[i]][1], ' type=matrix names=FALSE as.numeric=TRUE>'))
						control_point_matrix <- matrix(json_list$control_points[[i]][2:length(json_list$control_points[[i]])], ncol=2, byrow=TRUE)
						for(j in 1:nrow(control_point_matrix)){
							shapes_file_write <- c(shapes_file_write, paste0('\t\t\t', paste(control_point_matrix[j, ], collapse='\t')))
						}
						shapes_file_write <- c(shapes_file_write, paste0('\t\t</', json_list$control_points[[i]][1], '>'))
					}					
					shapes_file_write <- c(shapes_file_write, '\t</curves.control>')
					save_shapes <- TRUE
				}

				# SAVE CURVE POINTS (IF FILEPATH GIVEN)
				if(!is.null(curve_points_file) || !is.null(shapes_file)){
					
					curve_points <- matrix(NA, nrow=0, ncol=2)
					curves_pixels <- list()
					num_curve_sets <- 0

					if(length(json_list$control_points) > 0){
						for(i in 1:length(json_list$control_points)){

							m <- matrix(as.numeric(json_list$control_points[[i]][2:length(json_list$control_points[[i]])]), nrow=(length(json_list$control_points[[i]])-1)/2, ncol=2, byrow=TRUE)

							if(nrow(m) <= 2) next

							points_on_bezier <- pointsOnBezier(p=m, method='adjoining', deg=2)

							# CIRCUMVENT BUG IN POINTSONBEZIER WHERE LAST POINT OVERSHOOTS BY ONE
							# IF SECOND TO LAST POINT IS THE SAME AS THE LAST POINT OF M, GET RID OF LAST POINT
							if(sum(abs(points_on_bezier$points[nrow(points_on_bezier$points)-1, ] - m[nrow(m), ])) == 0)
								points_on_bezier$points <- points_on_bezier$points[-nrow(points_on_bezier$points), ]

							curves_pixels[[json_list$control_points[[i]][1]]] <- points_on_bezier$points

							rownames(points_on_bezier$points) <- paste0(
								json_list$control_points[[i]][1], 
								formatC(1:nrow(points_on_bezier$points), width=4, format="d", flag="0"))

							curve_points <- rbind(points_on_bezier$points, curve_points)
							num_curve_sets <- num_curve_sets + 1
						}
					}

					if(!is.null(curve_points_file)) write.table(curve_points, file = curve_points_file, quote=F, sep="\t", col.names=F, row.names=T)
					if(!is.null(shapes_file) && nrow(curve_points) > 0){

						shapes_file_write <- c(shapes_file_write, '\n\t<curves.pixel type=matrix rownames=TRUE as.numeric=TRUE>')
						for(name in names(curves_pixels)){
							shapes_file_write <- c(shapes_file_write, paste0('\t\t<', name, ' type=matrix rownames=FALSE as.numeric=TRUE>'))
							for(j in 1:nrow(curves_pixels[[name]])){
								shapes_file_write <- c(shapes_file_write, paste0('\t\t\t', paste(curves_pixels[[name]][j, ], collapse='\t')))
							}
							shapes_file_write <- c(shapes_file_write, paste0('\t\t</', name, '>'))
						}
						shapes_file_write <- c(shapes_file_write, '\t</curves.pixel>')

						save_shapes <- TRUE

						# SAVE SCALED CURVES
						if(apply_scaling){
							shapes_file_write <- c(shapes_file_write, '\n\t<curves.scaled type=matrix rownames=TRUE as.numeric=TRUE>')
							for(name in names(curves_pixels)){
								shapes_file_write <- c(shapes_file_write, paste0('\t\t<', name, ' type=matrix rownames=FALSE as.numeric=TRUE>'))
								for(j in 1:nrow(curves_pixels[[name]])){
									shapes_file_write <- c(shapes_file_write, paste0('\t\t\t', paste0(curves_pixels[[name]][j, 1]*scaling_wpp, '\t', -curves_pixels[[name]][j, 2]*scaling_wpp)))
								}
								shapes_file_write <- c(shapes_file_write, paste0('\t\t</', name, '>'))
							}
							shapes_file_write <- c(shapes_file_write, '\t</curves.scaled>')
						}
					}

					curve_save_status <- '(control and curve points)'

					# NOTIFICATION IN R CONSOLE
					#if(length(num_curve_sets) == 1){plural <- ''}else{plural <- 's'}
					#cat(paste0(num_curve_sets, " set", plural, " of curve point", plural, " saved to '", control_points_file, "'\n"))
				}				

				# NOTIFICATION IN R CONSOLE
				#cat(paste0(length(json_list$control_points), " set", plural, " of control points saved to '", control_points_file, "'\n"))

				if(length(json_list$control_points) == 1){plural <- ''}else{plural <- 's'}
				out$update_status <- paste0(out$update_status, length(json_list$control_points), " curve", plural, " saved ", curve_save_status, ". ")
			}
			
			# RULER POINT SUBMISSION
			if(!is.null(json_list$ruler_points) && !is.null(json_list$submit_shapes)){
				
				if(sum(!is.na(ruler_points[, 1])) > 0){
				
					# REMOVE NA ROWS
					ruler_points <- matrix(ruler_points[!is.na(ruler_points[, 1]), ], ncol=2, dimnames=list(rownames(ruler_points)[!is.na(ruler_points[, 1])], NULL))

					if(!is.null(shapes_file)){
						shapes_file_write <- c(shapes_file_write, '\n\t<ruler.points type=matrix rownames=TRUE as.numeric=TRUE>')
						for(i in 1:nrow(ruler_points)) shapes_file_write <- c(shapes_file_write, paste0('\t\t', paste(c(rownames(ruler_points)[i], ruler_points[i, ]), collapse='\t')))
						shapes_file_write <- c(shapes_file_write, '\t</ruler.points>')
						save_shapes <- TRUE
					}

					if(nrow(ruler_points) == 1){plural <- ''}else{plural <- 's'}
					out$update_status <- paste0(out$update_status, nrow(ruler_points), " ruler point", plural, " saved. ")
				}
			}

			# CHECKERBOARD CORNER SUBMISSION
			if(!is.null(json_list$corners) && !is.null(json_list$submit_shapes)){

				# READ CORNERS TO MATRIX
				corners <- matrix(as.numeric(cbind(sapply(json_list$corners, "[[", 1), sapply(json_list$corners, "[[", 2))), 
					nrow=length(json_list$corners), ncol=2)

				shapes_file_write <- c(shapes_file_write, '\n\t<checker.pixel type=matrix rownames=FALSE as.numeric=TRUE>')
				if(nrow(corners) > 0) for(i in 1:nrow(corners)) shapes_file_write <- c(shapes_file_write, paste0('\t\t', paste(corners[i, ], collapse='\t')))
				shapes_file_write <- c(shapes_file_write, '\t</checker.pixel>')

				if(nrow(corners) == 1){plural <- ''}else{plural <- 's'}
				out$update_status <- paste0(out$update_status, nrow(corners), " checkerboard corner", plural, " saved. ")
			}

			# SAVE TO SHAPES FILE
			if(!is.null(shapes_file) && save_shapes){

				# ADD EXTRA INFORMATION
				if(!is.null(json_list$checker_square_world) && json_list$checker_square_world != "NA") 
					shapes_file_write <- c(paste0("\t<square.size asnumeric=TRUE >", json_list$checker_square_world, "</square.size>"), shapes_file_write)

				if(!is.null(json_list$checker_square_pixel) && json_list$checker_square_pixel != "NA") 
					shapes_file_write <- c(paste0("\t<square.pixel asnumeric=TRUE >", json_list$checker_square_pixel, "</square.pixel>"), shapes_file_write)

				if(!is.null(json_list$checkerboard_ny) && json_list$checkerboard_ny != "NA") 
					shapes_file_write <- c(paste0("\t<checkerboard.ny asnumeric=TRUE >", json_list$checkerboard_ny, "</checkerboard.ny>"), shapes_file_write)

				if(!is.null(json_list$checkerboard_nx) && json_list$checkerboard_nx != "NA") 
					shapes_file_write <- c(paste0("\t<checkerboard.nx asnumeric=TRUE >", json_list$checkerboard_nx, "</checkerboard.nx>"), shapes_file_write)

				if(!is.null(json_list$ruler_interval_pixels) && json_list$ruler_interval_pixels != "NA") 
					shapes_file_write <- c(paste0("\t<ruler.pixel asnumeric=TRUE >", json_list$ruler_interval_pixels, "</ruler.pixel>"), shapes_file_write)

				if(!is.null(json_list$ruler_interval_world) && json_list$ruler_interval_world != "NA") 
					shapes_file_write <- c(paste0("\t<ruler.interval asnumeric=TRUE >", json_list$ruler_interval_world, "</ruler.interval>"), shapes_file_write)

				if(!is.null(json_list$scaling_units) && json_list$scaling_units != "NA") 
					shapes_file_write <- c(paste0("\t<scaling.units>", json_list$scaling_units, "</scaling.units>"), shapes_file_write)

				if(apply_scaling) shapes_file_write <- c(paste0("\t<scaling asnumeric=TRUE >", scaling_wpp, "</scaling>"), shapes_file_write)

				shapes_file_write <- c(
					'<shapes>',
					paste0("\t<image.name>", json_list$image_name, "</image.name>"), 
					paste0("\t<image.id>", json_list$image_id, "</image.id>"), 
					"\t<class>shapes</class>", 
					shapes_file_write,
					'</shapes>')
				write(shapes_file_write, file=shapes_file)
				#cat('Write file.\n')
				#cat(out$update_status, '\n')
			}

			# SAVE SETTINGS
			settings <- list(
				'copy_landmarks' = json_list$copy_landmarks,
				'copy_curves' = json_list$copy_curves,
				'copy_ruler_points' = json_list$copy_ruler_points,
				'copy_corners' = json_list$copy_corners,
				'copy_scaling' = json_list$copy_scaling
			)

			# DETECT IMAGE CHANGE
			if(!is.null(json_list$change_image)){
				if(json_list$change_image == 1) stopApp(list('next.command'='next', 'settings'=settings))
				if(json_list$change_image == -1) stopApp(list('next.command'='prev', 'settings'=settings))
			}

			# DETECT EXIT
			if(!is.null(json_list$exit)) stopApp(list('next.command'='exit', 'settings'=settings))

		}

		# ADD OUTPUTS
		if(!is.null(json_list$find_ruler_interval)) out$ruler_interval_pixels <- ruler_interval_pixels
		if(!is.null(json_list$find_checkerboard_corners)){

			out$corners <- list()
			for(i in 1:nrow(corners)) out$corners[[i]] <- corners[i, ]

			out$checker_square_pixel <- checker_square_pixel
		}

		#if(out$update_status != "" || !is.null(json_list$find_ruler_interval)){
		#	#print(out)
		#	#out_str <- '{"update_status":"8 landmarks saved. 2 curves saved (control and curve points). "}'
		#
		#}else{
		#	out_str <- ""
		#}
		out_str <- listToJSONStr(out, direct=TRUE)

		#cat(out_str)
		#cat("\n")

		out_str
	})
})

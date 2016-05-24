reconstructStereoSets <- function(shapes.2d, shapes.3d, cal.file, 
	set.names = NULL, min.common = 3, unify = TRUE, reconstruct.curves = TRUE, 
	even.spacing = NULL, print.progress = TRUE, verbose = FALSE, update.only = FALSE, 
	min.direct.tangency = 25, min.fill.tangency = 10, epi.err.weight = 0, 
	rec.err.weight = 1, curves.as.landmarks = FALSE, curve.name.width = 5){

	# RECONSTRUCTS 2D SHAPE DATA FROM MULTIPLE VIEWS AND ASPECTS INTO 3D SHAPES

	# CHECK THAT SHAPE FILES EXIST
	if(!file.exists(shapes.2d)) stop("File '", shapes.2d, "' not found.")

	# IF SHAPES.FILE DOESNT EXIST, CREATE IT
	if(!file.exists(shapes.3d)) dir.create(shapes.3d)

	# FIND SHAPE FILE DIRECTORIES
	shape_fdir <- list.files(shapes.2d)
	shape_fdir <- shape_fdir[!grepl('[.]txt$', shape_fdir)]

	# FIND SHAPE FILES FOR EACH VIEW
	files_2d <- list()
	for(view in shape_fdir){
		files_2d[[view]] <- list.files(paste0(shapes.2d, '/', view))
	}
	
	# GET LIST OF ALL UNIQUE FILENAMES
	files_2d_unique <- unique(unlist(files_2d))

	# CREATE MATRIX WITH FILE IN EACH VIEW
	files_in_view <- matrix(FALSE, nrow=length(files_2d_unique), ncol=length(shape_fdir), 
		dimnames=list(files_2d_unique, shape_fdir))
	
	# FILL MATRIX
	for(view in names(files_2d)) files_in_view[files_2d[[view]], view] <- TRUE

	# REPORT WARNING IF FILES NOT FOUND IN AT LEAST TWO VIEWS
	#files_not_two <- rownames(files_in_view)[rowSums(files_in_view) < 2]
	#warning(paste0("The following shape files were found in less than 2 views:\n\t", paste(files_not_two, collapse="\n\t"), ""))

	# USE ONLY MATCHED FILES FOR RECONSTRUCTION
	files_2d_matched <- rownames(files_in_view)[rowSums(files_in_view) > 1]

	# REMOVE .TXT EXTENSION
	files_2d_wo_ext <- gsub('[.]txt$', '', files_2d_matched, ignore.case=TRUE)

	# REMOVE ASPECT, IF PRESENT
	files_2d_wo_asp <- gsub('[_]?a[0-9]+$', '', files_2d_wo_ext, ignore.case=TRUE)

	# IF set.names IS SPECIFIED, REMOVE FILES NOT IN set.names
	if(!is.null(set.names)){

		# REMOVE EXTENSION AND ASPECT IF PRESENT
		set.names <- gsub('[.]txt$', '', set.names, ignore.case=TRUE)
		set.names <- gsub('[_]?a[0-9]+$', '', set.names, ignore.case=TRUE)

		# CHECK THAT ALL NAMES ARE FOUND
		not_found <- set.names[!set.names %in% files_2d_wo_asp]
		if(length(not_found) > 0) stop(paste0("The following names in 'set.names' were not found in '", shapes.2d, "': ", paste(not_found, collapse=", "), ""))

		files_2d_matched <- files_2d_matched[files_2d_wo_asp %in% set.names]
		files_2d_wo_ext <- files_2d_wo_ext[files_2d_wo_asp %in% set.names]
		files_2d_wo_asp <- files_2d_wo_asp[files_2d_wo_asp %in% set.names]
	}

	#print(files_2d_matched);print(files_2d_wo_ext);print(files_2d_wo_asp)
	
	# IF update.only, REMOVE CHILD FILES THAT ARE NOT MORE RECENT THAN THEIR PARENT FILES
	if(update.only){
		
		# GET PARENT FOLDERS
		files_3d <- gsub('[.]txt$', '', list.files(shapes.3d), ignore.case=TRUE)
		
		# FIND CHILD FILES FOR EACH
		if(length(files_3d) > 0){
			for(i in 1:length(files_3d)){
			
				# REMOVE ASPECT, IF PRESENT
				file_3d_wo_asp <- gsub('[_]?a[0-9]+$', '', files_3d[i], ignore.case=TRUE)

				# FIND CHILDREN FILES
				match_children <- file_3d_wo_asp == files_2d_wo_asp

				# SKIP IF NO CHILDREN FILES
				if(sum(match_children) == 0) next

				# GET CHILDREN FILE PATHS
				children_fpaths <- c()
				for(view in names(files_2d)) children_fpaths <- c(children_fpaths, paste0(view, '/', files_2d_matched[match_children]))

				# FIND MODIFIED DATE OF PARENT FILE
				mtime_parent <- file.info(paste0(shapes.3d, '/', files_3d[i], '.txt'))$mtime

				# FIND MODIFIED DATES OF CHILDREN FILES
				mtime_children <- file.info(paste0(shapes.2d, '/', children_fpaths))$mtime

				# SKIP (RETAIN) IF ANY CHILDREN ARE OLDER THAN PARENT
				if(sum(mtime_children > mtime_parent, na.rm=TRUE) > 0) next

				# REMOVE SINCE PARENT MUST BE OLDER THAN ALL CHILDREN
				files_2d_matched <- files_2d_matched[!match_children]
				files_2d_wo_ext <- files_2d_wo_ext[!match_children]
				files_2d_wo_asp <- files_2d_wo_asp[!match_children]
			}
		}		
	}
	
	if(print.progress) cat('reconstructStereoSets\n')

	# IF NO FILES TO 
	if(length(files_2d_matched) == 0){
		if(print.progress) cat('\tAll files up to date.\n')
		return(NULL)
	}
	
	# GET CALIBRATION COEFFICIENTS
	if(!file.exists(cal.file)) stop(paste0("The calibration file ('", cal.file, "') was not found."))
	cal_coeff <- XML4R2list(cal.file)$calibration$cal.coeff
	
	# IF NO COLUMN NAMES, USE 2D DIRECTORY NAMES
	if(is.null(colnames(cal_coeff))) colnames(cal_coeff) <- names(files_2d)

	# VECTOR OF WRITTEN FILES
	files.written <- c()
	
	# GET EVEN.SPACING INPUT
	if(is.list(even.spacing)){
		
		# CONVERT TO NAMED VECTOR
		even.spacing <- setNames(unlist(even.spacing), names(even.spacing))
		
	}else if(length(even.spacing) == 1 && grepl('[.]txt$', even.spacing[1])){
		
		# READ IN MATRIX FILE
		read_table <- as.matrix(read.table(file=even.spacing, sep='\t', quote="", row.names=1))
		even.spacing <- setNames(c(read_table), rownames(read_table))
	}

	if(is.null(names(even.spacing)) && length(even.spacing) > 1) stop("If 'even.spacing' is vector of length greater than 1 the vector must be named.")

	#print(files_2d_matched);print(files_2d_wo_ext);print(files_2d_wo_asp)

	# GET 3D FILENAMES
	files_2d_wo_asp_unique <- unique(files_2d_wo_asp)

	# RECONSTRUCT AND, IF POSSIBLE, UNIFY
	for(fname in files_2d_wo_asp_unique){

		if(print.progress) cat('\t', fname, '\n', sep='')

		# FIND MATCHING ASPECTS
		match_aspects <- files_2d_matched[files_2d_wo_asp == fname]
		
		# CREATE LIST FOR 3D COORDINATES FROM EACH ASPECT
		all_aspects <- list()
		aspect_names <- c()
		all_curve_pts <- list()
		all_landmark_names <- c()
		unified_lm <- list()

		# CHECK WHETHER MULTIPLE ASPECTS
		mult_asp <- TRUE
		tabs <- '\t\t\t'
		if(length(match_aspects) == 1 && files_2d_wo_asp[files_2d_wo_asp == fname] == files_2d_wo_ext[files_2d_wo_asp == fname]){
			mult_asp <- FALSE
			tabs <- '\t\t'
		}

		for(aspect in 1:length(match_aspects)){

			if(mult_asp && print.progress) cat('\t\tAspect ', aspect, '\n', sep='')
			
			# GET FILE PATHS TO ALL VIEWS
			view_fpaths <- paste0(shapes.2d, "/", shape_fdir, "/", match_aspects[aspect])
			
			# SET VIEWS THAT EXIST
			views <- shape_fdir[file.exists(view_fpaths)]

			# 3D COORDINATE MATRIX
			coor_3d <- matrix(NA, nrow=0, ncol=3)

			# READ SHAPES
			read_shapes <- readShapes(view_fpaths[file.exists(view_fpaths)])
			#print(names(read_shapes))

			if(!is.null(read_shapes$landmarks.pixel)){

				# RECONSTRUCT LANDMARKS
				dlt_recon <- dltReconstruct(cal_coeff[, dimnames(read_shapes$landmarks.pixel)[[3]]], read_shapes$landmarks.pixel)

				if(print.progress){
					cat(tabs, 'Landmark Reconstruction RMS Error: ', mean(dlt_recon$rmse, na.rm=TRUE), ' px\n', sep='')
					if(verbose) cat(tabs, '\t', paste(rownames(dlt_recon$coor.3d), '\t', paste0(format(dlt_recon$rmse), ' px'), collapse=paste0('\n\t', tabs)), '\n', sep='')
				}
				
				all_landmark_names <- c(all_landmark_names, rownames(dlt_recon$coor.3d))
				coor_3d <- rbind(coor_3d, dlt_recon$coor.3d)
			}

			# CREATE VECTOR OF CURVE POINT ROWNAMES TO SEPARATE LATER
			curve_pt_rownames <- c()

			if(reconstruct.curves && !is.null(read_shapes$curves.pixel)){
			
				if(print.progress) cat(tabs, 'Curve point matching and reconstruction\n', sep='')

				# GET ALL CURVE NAMES
				curve_names <- c()
				for(view in names(read_shapes$curves.pixel)) curve_names <- c(curve_names, names(read_shapes$curves.pixel[[view]]))
				curve_names <- sort(unique(curve_names))
				
				# GET VIEW NAMES
				view_names <- names(read_shapes$curves.pixel)
				
				# CREATE MATRIX OF CURVES IN VIEW
				curves_in_view <- matrix(FALSE, nrow=length(curve_names), ncol=length(read_shapes$curves.pixel), 
					dimnames=list(curve_names, view_names))
				
				# FILL MATRIX
				for(view in names(read_shapes$curves.pixel)) curves_in_view[names(read_shapes$curves.pixel[[view]]), view] <- TRUE

				# RECONSTRUCT CURVES FOR EACH VIEW
				for(i in 1:nrow(curves_in_view)){
					
					# SKIP IF IN LESS THAN 2 VIEWS
					if(sum(curves_in_view[curve_names[i], ]) < 2) next
					
					# FIND FIRST TWO TRUE VIEWS
					in_views <- view_names[curves_in_view[curve_names[i], ] == TRUE][1:2]

					# CREATE LIST OF CURVES
					lm.list <- list(read_shapes$curves.pixel[[in_views[1]]][[curve_names[i]]], read_shapes$curves.pixel[[in_views[2]]][[curve_names[i]]])

					# FIND CORRESPONDING POINTS ON CURVES BETWEEN VIEWS
					#if(print.progress) cat(tabs, '\t', curve_names[i], ' (', in_views[1], '/', in_views[2], ')', sep='')
					if(print.progress) cat(tabs, '\t', curve_names[i], '\n', sep='')
					
					# FIND MATCHING CURVE POINTS
					dlt_mcp <- dltMatchCurvePoints(lm.list=lm.list, cal.coeff=cal_coeff[, in_views], 
						min.direct.tangency=min.direct.tangency, min.fill.tangency=min.fill.tangency,
						epi.err.weight=epi.err.weight, rec.err.weight=rec.err.weight)
				
					if(print.progress){
						dlt_mcp_summary <- summary(dlt_mcp, print.tab=paste0(tabs, '\t'))

						if(verbose){
							cat(c(dlt_mcp_summary[4:length(dlt_mcp_summary)]), sep='')
						}else{
							#cat(c(dlt_mcp_summary[4:6]), sep='')
						}
					}

					# CREATE ROWNAMES
					row_names <- paste0(curve_names[i], '_', formatC(1:nrow(dlt_mcp$match.lm.list[[1]]), width=curve.name.width, format="d", flag="0"))
				
					# ADD ROWNAMES TO VECTOR
					curve_pt_rownames <- c(curve_pt_rownames, row_names)

					# ADD ROWNAMES TO MATRICES
					rownames(dlt_mcp$match.lm.list[[1]]) <- row_names
					rownames(dlt_mcp$match.lm.list[[2]]) <- row_names

					# CREATE ARRAY FOR RECONSTRUCTION
					cp_array <- array(NA, dim=c(length(row_names), 2, length(shape_fdir)), dimnames=list(row_names, NULL, shape_fdir))

					for(k in 1:length(dlt_mcp$match.lm.list)) cp_array[, , k] <- dlt_mcp$match.lm.list[[k]]
				
					# RECONSTRUCT CURVE POINTS
					dlt_recon <- dltReconstruct(cal_coeff[, in_views], cp_array)

					if(print.progress){
						if(!verbose) cat(tabs, '\t\tMax Reconstruction error: ', max(dlt_recon$rmse, na.rm=TRUE), ' px\n', sep='')
						if(verbose) cat(tabs, '\t\tCurve Point Reconstruction RMS Error, Mean: ', mean(dlt_recon$rmse, na.rm=TRUE), ' px; Min: ', min(dlt_recon$rmse, na.rm=TRUE), ' px; Max: ', max(dlt_recon$rmse, na.rm=TRUE), ' px\n', sep='')
					}

					# GET CURVE POINTS
					curve_points <- dlt_recon$coor.3d
			
					# FIND EVENLY SPACED POINTS ALONG CURVE IF SPECIFIED
					if(!is.null(even.spacing)){
						if(is.null(names(even.spacing))){

							curve_points <- pointsAtEvenSpacing(x=curve_points, n=even.spacing[1])
						}else{

							if(!curve_names[i] %in% names(even.spacing)) stop(paste0("Number of curve points for curve '", curve_names[i], "' not found in 'even.spacing'."))

							curve_points <- pointsAtEvenSpacing(x=curve_points, n=even.spacing[curve_names[i]])
						}
					}

					if(curves.as.landmarks){

						# IF ADDING CURVE POINTS TO LANDMARK MATRIX, ADD CURVE NAMES AS LANDMARK NAMES
						all_landmark_names <- c(all_landmark_names, rownames(curve_points))

					}else{

						# SAVE CURVE POINT NAMES
						all_curve_pts[[curve_names[i]]] <- rownames(curve_points)
					}

					coor_3d <- rbind(coor_3d, curve_points)
					
					#return(1)
				}
			}

			# DONT CREATE LIST ELEMENT IF NO 3D COORDINATES RECONSTRUCTED
			if(nrow(coor_3d) == 0) next
			
			# SAVE 3D COORDINATES TO LIST
			all_aspects[[length(all_aspects)+1]] <- coor_3d
		}

		if(length(all_aspects) == 0){

			if(print.progress) cat('\t\tNo shapes found.\n', sep='')
			next

		}else if(length(all_aspects) == 1){

			# SAVE UNIFIED 3D SET
			unified_lm[[1]] <- all_aspects[[1]]
			shapes_3d_fname <- fname

			if(print.progress && mult_asp && unify) cat('\t\tOnly one aspect, no unification performed.\n', sep='')

		}else{

			unify_lm <- NULL
			
			if(unify){

				# GET ALL UNIQUE LANDMARK NAMES
				all_landmarks <- unique(unlist(lapply(all_aspects, 'rownames')))

				# CREATE ARRAY FOR 3D LANDMARKS FROM ALL ASPECTS
				all_aspects_arr <- array(NA, dim=c(length(all_landmarks), 3, length(all_aspects)), dimnames=list(all_landmarks, NULL, match_aspects))
		
				# FILL ARRAY
				for(aspect in 1:length(all_aspects)) all_aspects_arr[rownames(all_aspects[[aspect]]), , aspect] <- all_aspects[[aspect]]

				if(print.progress) cat('\t\tUnify landmarks', sep='')

				# UNIFY ASPECTS
				unify_lm <- unifyLandmarks(all_aspects_arr, min.common=min.common, return.on.error=TRUE)

				if(is.null(unify_lm)) cat(paste0("\n\t\t\tNumber of common points less than min.common (", min.common, "). No unification performed.\n"))
			}
			
			if(is.null(unify_lm)){

				# SAVE UNIFIED 3D SET
				unified_lm <- all_aspects
				shapes_3d_fname <- gsub('[.]txt$', '', match_aspects, ignore.case=TRUE)

			}else{

				if(print.progress){
					unify_summary <- summary(unify_lm, print.tab='\t\t', verbose=verbose)
					cat(c('\n', unify_summary[4:(length(unify_summary)-1)]), sep='')
					cat('\n')
				}

				# SAVE UNIFIED 3D SET
				unified_lm[[1]] <- unify_lm$lm.matrix
				shapes_3d_fname <- fname
			}
		}

		# GET 3D SHAPE FILE FILE PATH
		shapes_3d_fpath <- paste0(shapes.3d, '/', shapes_3d_fname, '.txt')
		
		# GET UNIQUE LANDMARK NAMES
		all_landmark_names_unique <- unique(all_landmark_names)

		for(i in 1:length(unified_lm)){

			# CREATE LIST STRUCTURE FOR CURVE POINTS
			curves_r <- list()
			if(length(all_curve_pts) > 0){
				for(j in 1:length(all_curve_pts)){
			
					# GET CURVE NAME
					curve_name <- names(all_curve_pts)[j]
					
					# IF CURVE POINTS NOT IN MATRIX, SKIP
					if(sum(!all_curve_pts[[j]] %in% rownames(unified_lm[[i]])) > 0) next

					# GET CURVE POINTS
					curve_points <- unified_lm[[i]][all_curve_pts[[j]], ]
					rownames(curve_points) <- NULL

					# ADD TO LIST
					curves_r[[curve_name]] <- curve_points
				}
			}
		
			# ADD FILE PATH TO WRITTEN FILES
			files.written <- c(files.written, shapes_3d_fpath[i])
			
			landmarks_in <- all_landmark_names_unique[all_landmark_names_unique %in% rownames(unified_lm[[i]])]
		
			# SAVE 3D COORDINATES TO SHAPE FILE
			list2XML4R(list=
				list('shapes'=
					list(
						'landmarks'=unified_lm[[i]][sort(landmarks_in), ],
						'curves'=curves_r
					)
				), file=shapes_3d_fpath[i])
		}
	}

	#rlist <- list(
	#	files.written=files.written
	#)
	#class(rlist) <- 'reconstructStereoSets'
	
	return(NULL)
}
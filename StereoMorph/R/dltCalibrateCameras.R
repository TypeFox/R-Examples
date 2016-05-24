dltCalibrateCameras <- function(coor.2d, nx, grid.size, c.run = FALSE, reduce.grid.dim = 3, 
	fit.min.break = 1, nlm.iter.max.init = 100, objective.min.init = 10, 
	nlm.eval.max = 350, nlm.iter.max = 250, nlm.calls.max = 100, min.views = 'max', 
	objective.min = 1, grid.incl.min = 2, objective.min.break = NULL, start.param = NULL, 
	sx = NULL, sy = NULL, print.progress = FALSE, print.tab = ''){

	# objective.min.break OBJECTIVE AT WHICH TO STOP ADDING CHECKERBOARDS AND START OVER WITH NEW SET

	# SET MIN NUMBER OF VIEWS TO USE IN ESTIMATING CALIBRATION COEFFICIENTS
	if(min.views == 'max'){min_views <- dim(coor.2d)[4]}else{min_views <- min.views}

	# FIND NUMBER OF NON-NA VIEWS FOR EACH ASPECT
	aspect_non_na <- rowSums(apply(!is.na(coor.2d), c(3, 4), 'sum') > 0)

	# REMOVE SETS WITH LESS THAN THE MINIMUM NUMBER OF VIEWS
	coor.2d <- coor.2d[, , aspect_non_na >= min_views, ]

	# FIND SECOND GRID DIMENSION
	ny <- dim(coor.2d)[1]/nx

	# REDUCE GRID DIMENSIONS
	if(reduce.grid.dim){

		# CHECK THAT reduce.grid.dim IS GREATER THAN TWO
		if(reduce.grid.dim < 3) stop(paste0("reduce.grid.dim (", reduce.grid.dim, ") must be greater than 2."))
		
		# SET REDUCED GRID DIMENSIONS
		rx <- reduce.grid.dim
		ry <- reduce.grid.dim
		
		# SET REDUCED GRID SIZES
		sx <- ((nx-1)*grid.size[1]) / (rx-1)
		sy <- ((ny-1)*grid.size[1]) / (ry-1)

		# EMPTY REDUCED GRID DIMENSION ARRAY
		coor_2d_red <- array(NA, dim=c(rx*ry, 2, dim(coor.2d)[3], dim(coor.2d)[4]))

		if(print.progress) cat('\n', print.tab, 'Reduce grid point density (', dim(coor.2d)[3]*dim(coor.2d)[4], ' total)\n', sep='')

		for(i in 1:dim(coor.2d)[3]){
			for(j in 1:dim(coor.2d)[4]){
				if(print.progress) cat(print.tab, '\t', (i-1)*dim(coor.2d)[4] + j, ') Aspect ', i, ', View ', j, '; ', sep='')
				coor_2d_red[, , i, j] <- resampleGridImagePoints(pts=coor.2d[, , i, j], nx=nx, rx=rx, ry=ry, fit.min.break=fit.min.break, print.progress=print.progress)$pts
			}
		}
	}else{
		
		# MAINTAIN FULL GRID DIMENSIONS
		rx <- nx;ry <- ny;
		
		#
		if(is.null(sx)) sx <- grid.size[1]
		if(is.null(sy) && length(grid.size) > 1) sy <- grid.size[2]
		
		# COPY POINTS TO REDUCED 2D COOR ARRAY
		coor_2d_red <- coor.2d
	}

	# SCALE INITIAL TRANSLATE PARAMETER TO REAL-WORLD UNITS (APPROX THIRD THE MAX DIMENSION OF GRID)
	t_init <- (max(nx, ny)*grid.size[1])/3
	
	# SET INITIAL TIME POINT
	ptm <- proc.time()

	# RUN NLM FUNCTION TO FIND TRANSFORMATION PARAMETERS THAT MINIMIZE INTERNAL RMS CALIBRATION ERROR
	if(print.progress) cat('\n', print.tab, 'Full Transform RMSE Minimization\n', print.tab,'Number of parameters: ', (dim(coor.2d)[3]-1)*6, '\n', print.tab, 'Number of points: ', rx*ry*dim(coor.2d)[3], '\n', sep='')

	# SET FIXED START PARAMETERS, FIX FIRST GRID AT ORIGIN
	p_fix <- c()

	for(num_grid in grid.incl.min:dim(coor.2d)[3]){

		if(print.progress) cat('\n', print.tab, 'Running minimization with ', num_grid, ' grids...', sep='')

		# ESTIMATE STARTING PARAMETERS BASED ON 2D VIEW
		estimate_start_params <- dltCCEstimateStartParams(coor.2d, num_grid, nx, ny, grid.size[1], p_fix)
		
		start_param_array <- estimate_start_params$start_param_array
		max_disp_rwu <- estimate_start_params$max_disp_rwu

		# GET 2D COORDINATE SUBSET MATRIX FOR TRANSFORMATION OPTIMIZATION
		coor_2d_t <- apply(coor_2d_red[, , 1:num_grid, ], c(2, 4), matrix, byrow=FALSE)

		lower <- rep(c(rep(-pi, 3), rep(-7*max(nx, ny)*grid.size[1], 3)), num_grid)

		# SET NLM CALL VARIABLES
		nlm_min <- rep(NA, nlm.calls.max)
		nlm_calls <- list()

		# SAVE PROCESSING TIME
		nlm_start <- proc.time()

		for(nlm_n in 1:nlm.calls.max){

			# SET STARTING PARAMETERS
			if(is.null(start.param)){

				# VARY EACH ITERATION
				start <- c(t(start_param_array[nlm_n, ((length(p_fix)/6)+1):(num_grid-1), ]))
				if(length(p_fix)) start <- c(p_fix, start)
			}else{
				start <- start.param[1:((num_grid-1)*6)]
				if(length(p_fix)) start[1:length(p_fix)] <- p_fix
				
				if(((num_grid-1)*6) > length(start.param)){
					start <- c(t(start_param_array[nlm_n, ((length(p_fix)/6)+1):(num_grid-1), ]))
					if(length(p_fix)) start <- c(p_fix, start)
				}

				# VARY LAST SIX PARAMETERS AFTER FIRST ITERATION
				if(nlm_n > 1){
					start[(length(start)-6+1):length(start)] <- sample(seq(0.5, 1.5, length=30), size=6)*start[(length(start)-6+1):length(start)]
				}
			}

			if(print.progress) cat(' ', nlm_n, sep='')

			#if(!is.null(sink.trace)) sink(paste0(sink.trace, "/", "num grids ", num_grid, " nlm iter ", nlm_n, ".txt"))

			start_full <- start
			stopped <- FALSE

			if(!is.null(nlm.iter.max.init)){
				nlm_fit <- tryCatch(
					expr={
						nlminb(start=start, objective=dltTransformationParameterRMSError, 
							control=list(iter.max = nlm.iter.max.init, rel.tol=1e-6),
							lower=lower, upper=-lower, coor.2d=coor_2d_t, nx=rx, ny=ry, 
							sx=sx, sy=sy, p.fixed=rep(0, 6))
					},
					error=function(cond){if(print.progress) cat('F');return(NULL)},
					warning=function(cond) return(NULL)
				)

				if(is.null(nlm_fit)) next

				#cat(" nlm_fit$objective: ", nlm_fit$objective, ", ", nlm_fit$message, ", ", sep="")
				
			}else{
				nlm_fit <- list('objective' = 0, 'message' = '')
			}

			#if(!is.null(sink.trace)) sink()

			if(nlm_fit$objective > objective.min.init && nlm_fit$message == "relative convergence (4)"){
				stopped <- TRUE
			}else{

				start_full <- nlm_fit$par

				nlm_fit <- tryCatch(
					expr={
						nlminb(start=start_full, objective=dltTransformationParameterRMSError, 
							control=list(eval.max = nlm.eval.max, iter.max = nlm.iter.max), 
							lower=lower, upper=-lower, coor.2d=coor_2d_t, nx=rx, ny=ry, 
							sx=sx, sy=sy, p.fixed=rep(0, 6))
					},
					error=function(cond){if(print.progress) cat('N');return(NULL)},
					warning=function(cond) return(NULL)
				)
			}

			#
			if(!is.null(nlm_fit)){

				# GET MINIMUM
				if(print.progress){
					if(nlm_fit$message == 'function evaluation limit reached without convergence (9)' || nlm_fit$message == 'iteration limit reached without convergence (10)'){
						cat('L')
					}else if(nlm_fit$message == 'false convergence (8)'){
						cat('F')
					}else if(nlm_fit$message == 'relative convergence (4)'){
						cat('R')
					}else{
						cat('C')
					}
					if(stopped) cat('-S')
					cat('(', round(nlm_fit$objective, 2), ')', sep='')
				}
				nlm_calls[[nlm_n]] <- nlm_fit
				nlm_min[nlm_n] <- nlm_calls[[nlm_n]]$objective

				# IF OPTIMIZATION MINIMUM IS LESS THAN objective.min, STOP ITERATING
				if(nlm_calls[[nlm_n]]$convergence == 0 && nlm_calls[[nlm_n]]$objective < objective.min) break
				if(nlm_fit$message == 'iteration limit reached without convergence (10)' && nlm_calls[[nlm_n]]$objective < objective.min) break
				if(nlm_fit$message == 'function evaluation limit reached without convergence (9)' && nlm_calls[[nlm_n]]$objective < objective.min) break
				if(nlm_fit$message == 'false convergence (8)' && nlm_calls[[nlm_n]]$objective < objective.min) break
			}
		}

		# SAVE PROCESSING TIME
		nlm_elapsed <- proc.time() - nlm_start

		# SAVE RUN WITH MINIMUM
		nlm_res_t <- nlm_calls[[which.min(nlm_min)]]

		if(print.progress){
			cat('\n', print.tab, '\tNumber of nlminb() calls: ', nlm_n, ' (', nlm.calls.max, ' max)', 
				'\n', print.tab, '\tTermination message: ', nlm_res_t$message, 
				'\n', print.tab, '\tMinimum: ', nlm_res_t$objective, 
				'\n', print.tab, '\tIterations: ', nlm_res_t$iterations, ' (', nlm.iter.max, ' max)', 
				'\n', print.tab, '\tFunction evaluations: ', nlm_res_t$evaluations['function'], ' (', nlm.eval.max, ' max)', 
				'\n', print.tab, '\tRun-time: ', nlm_elapsed[1], ' sec', 
				'\n', print.tab, '\tStarting parameters: ', sep='')
				cat(round(start, 5), sep=', ')
				cat('\n', print.tab, '\tFinal estimated parameters: ', sep='')
				cat(round(nlm_res_t$par, 5), sep=', ')
				cat('\n', print.tab, '\tSum of absolute differences between initial and final parameters: ', sum(abs(start - nlm_res_t$par)), sep='')
			cat('\n')
		}

		# SET FIXED PARAMETERS FOR NEXT ITERATION
		p_fix <- nlm_res_t$par
		
		if(!is.null(objective.min.break)) if(nlm_res_t$objective > objective.min.break) break
	}

	# SAVE PROCESSING TIME
	run_time <- proc.time() - ptm

	if(print.progress) cat('\n', print.tab, 'Total estimation processing time: ', run_time[1], ' sec', sep='')

	# SAVE OPTIMIZED PARAMETERS
	p_init <- matrix(c(rep(0, 6), nlm_res_t$par), nrow=6)

	# GET 3D COORDINATES BASED ON OPTIMIZED PARAMETERS
	if(length(grid.size) == 1){
		coor_3d_coeff <- transformPlanarCalibrationCoordinates(tpar=c(p_init), nx=nx, ny=ny, sx=grid.size[1])
	}else{
		coor_3d_coeff <- transformPlanarCalibrationCoordinates(tpar=c(p_init), nx=nx, ny=ny, sx=grid.size[1], sy=grid.size[2])
	}

	# GET 2D INPUT COORDINATES
	coor_2d_coeff <- apply(coor.2d[, , 1:(length(p_init)/6), ], c(2, 4), matrix, byrow=FALSE)

	# GET DLT CALIBRATION COEFFICIENTS FROM OPTIMIZED 3D COORDINATE SUBSET (ALL GRID POINTS - SO RMSE CAN DIFFER FROM NLM MINIMUM)
	dlt_coefficients_t <- dltCoefficients(coor.3d=coor_3d_coeff, coor.2d=coor_2d_coeff)

	if(c.run){

		# GET 2D COORDINATE SUBSET MATRIX FOR COEFFICIENT OPTIMIZATION
		coor_2d_c <- apply(coor_2d_red, c(2, 4), matrix, byrow=FALSE)

		# SET INITIAL TIME POINT AND NUMBER OF ITERATIONS
		ptm <- proc.time();#c_iter <<- 0

		# RUN NLM FUNCTION TO FIND COEFFICIENTS THAT MINIMIZE INTERNAL RMS CALIBRATION ERROR
		if(print.progress) cat('\n', print.tab, 'Coefficient RMSE Minimization\n', print.tab, 'Number of parameters: ', length(c(dlt_coefficients_t$cal.coeff)), '\n', print.tab, 'Number of points: ', rx*ry*dim(coor.2d)[3], '\n', print.tab, 'Running minimization...\n', sep='')
		nlm_res_c <- nlm(dltCoefficientRMSError, p=c(dlt_coefficients_t$cal.coeff), coor.2d=coor_2d_c)

		if(print.progress) cat('\n', print.tab, 'Termination code: ', nlm_res_c$code, '\n', print.tab, 'Iterations: ', nlm_res_c$iterations, '\n', sep='')

		# SAVE PROCESSING TIME
		run_time_c <- proc.time() - ptm

		# USE DLT COEFFICIENTS TO RECONSTRUCT ALL 2D COORDINATES
		dlt_reconstruct <- dltReconstruct(cal.coeff=matrix(nlm_res_c$estimate, nrow=11, ncol=dim(coor.2d)[4]), coor.2d=apply(coor.2d, c(2, 4), matrix, byrow=FALSE))

		# USE DLT COEFFICIENTS TO RECONSTRUCT ALL 2D COORDINATES
		dlt_coefficients_c <- dltCoefficients(coor.3d=dlt_reconstruct$coor.3d, coor.2d=apply(coor.2d, c(2, 4), matrix, byrow=FALSE))

		mean_reconstruct_rmse <- mean(dlt_reconstruct$rmse)
		coor_3d <- dlt_reconstruct$coor.3d
		cal_coeff <- dlt_coefficients_c$cal.coeff
		coefficient_rmse <- dlt_coefficients_c$rmse
	}else{

		# USE DLT COEFFICIENTS TO RECONSTRUCT ALL 2D COORDINATES
		dlt_reconstruct <- dltReconstruct(cal.coeff=dlt_coefficients_t$cal.coeff, coor.2d=apply(coor.2d, c(2, 4), matrix, byrow=FALSE))

		mean_reconstruct_rmse <- mean(dlt_reconstruct$rmse)
		coor_3d <- dlt_reconstruct$coor.3d
		cal_coeff <- dlt_coefficients_t$cal.coeff
		coefficient_rmse <- dlt_coefficients_t$rmse

		#dlt_coefficients_c <- dltCoefficients2(coor.3d=coor_3d, coor.2d=apply(coor.2d, c(2, 4), matrix, byrow=FALSE))
		#hist(dlt_coefficients_c$rmse)
		
		nlm_res_c <- list('estimate'=NA, 'minimum'=NA)
		run_time_c <- NA
	}

	l <- list(
		cal.coeff=cal_coeff, 
		coor.3d=coor_3d, 
		mean.reconstruct.rmse=mean_reconstruct_rmse, 
		coefficient.rmse=coefficient_rmse,
		t.param.final=p_init,
		t.min=nlm_res_t$objective,
		t.runtime=run_time[1],
		c.param.init=c(dlt_coefficients_t$cal.coeff),
		c.param.final=nlm_res_c$estimate,
		c.min=nlm_res_c$minimum,
		c.iter=nlm_res_c$iterations,
		c.runtime=run_time_c[1]
		)
	class(l) <- 'dltCalibrateCameras'
	l
}

summary.dltCalibrateCameras <- function(object, ...){
	r <- ''

	r <- c(r, '\ndltCalibrateCameras Summary\n')

	r <- c(r, '\tMinimize RMS Error by transformation\n')
	r <- c(r, '\t\tTotal number of parameters estimated: ', length(object$t.param.final)-6, '\n')
	r <- c(r, '\t\tFinal Minimum Mean RMS Error: ', round(object$t.min, 3), '\n')
	r <- c(r, '\t\tTotal Run-time: ', round(object$t.runtime, 2), ' sec\n')

	if(!is.null(object$c.iter)){
		r <- c(r, '\tMinimize RMSE by coefficients\n')
		r <- c(r, '\t\tNumber of parameters: ', length(object$c.param.init), '\n')
		r <- c(r, '\t\tSum of absolute differences between initial and final parameters: ', format(sum(abs(object$c.param.init - object$c.param.final))), '\n')
		r <- c(r, '\t\tTotal function calls: ', object$c.iter, '\n')
		r <- c(r, '\t\tFinal Mean RMS Error: ', format(object$c.min), '\n')
		r <- c(r, '\t\tRun-time: ', round(object$c.runtime, 2), ' sec\n')
	}

	r <- c(r, '\n')

	r <- c(r, '\tMean Reconstruct RMSE: ', round(object$mean.reconstruct.rmse, 3), '\n')

	if(length(object$coefficient.rmse) == 1){
		r <- c(r, '\tMean Calibration Coefficient RMS Error: ', round(object$coefficient.rmse, 3), '\n')
	}else{
		r <- c(r, '\tCalibration Coefficient RMS Error:\n')
		for(i in 1:length(object$coefficient.rmse)){
			r <- c(r, '\t\t', i, ': ', round(object$coefficient.rmse[i], 3), '\n')
		}
	}

	r <- c(r, '\t3D Coordinate Matrix dimensions: ', dim(object$coor.3d)[1], ' x ', dim(object$coor.3d)[2], '\n')

	class(r) <- "summary.dltCalibrateCameras"
	r
}

print.summary.dltCalibrateCameras <- function(x, ...) cat(x, sep='')


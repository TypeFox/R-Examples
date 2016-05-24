bezierCurveFit <- function(m, min.control.points = 3, max.control.points = 20, fix.start.end = FALSE, max.rse = NULL, max.rse.percent.change = 0.01, na.fill = FALSE){
	# START AND END MUST BE NON-NA

	# MAKE SURE THAT MAX IS GREATER THAN MIN CONTROL POINTS
	if(min.control.points > max.control.points) stop(paste0("min.control.points (", min.control.points, ") must be less than or equal to max.control.points (", max.control.points, ")"))

	if(na.fill){
		# FIND NA VALUES BASED
		is_na <- is.na(m[, 1])
	
		# IF NA FILL
		i <- 1
		while(i < nrow(m)){
	
			# SKIP IF NOT NA
			if(!is_na[i]){i <- i + 1;next}
	
			# FIND NEXT NON-NA VALUE
			next_non_na <- i + which(!is_na[(i+1):nrow(m)])[1]
	
			# FIND POINTS ON LINE BETWEEN NON-NA VALUES
			m_fill <- matrix(m[next_non_na, ] - m[i-1, ], nrow=next_non_na-i+2, ncol=ncol(m), byrow=TRUE)*seq(0, 1, length=next_non_na-i+2) + matrix(m[i-1, ], nrow=next_non_na-i+2, ncol=ncol(m), byrow=TRUE)

			# ENTER VALUES INTO M MATRIX
			m[(i-1):next_non_na, ] <- m_fill
	
			# WHERE TO START SEARCH FOR NEXT NA
			i <- next_non_na
		}
	}

	# CONVERT INPUT PARAMETER POINTS TO MATRIX
	if(!is.matrix(m)) m <- matrix(m)

	# GET T PARAMETER, SAME SIZE AS MATRIX
	t <- seq(0, 1, length = nrow(m))
	
	# MODEL FIT PARAMETER LIST
	p <- list()
	
	# VECTOR OF FINAL RESIDUAL STANDARD ERRORS OF FIT FOR DIMENSION
	p.rse <- rep(NA, ncol(m))

	# REASON FOR STOPPING FIT ITERATIONS
	fit.stopped.by <- rep(NA, ncol(m))

	# FIT EACH DIMENSION SEPARATELY -- NLS() DOESNT WORK WITH MATRIX INPUT PARAMETER
	for(i in 1:ncol(m)){
		#cat(i, '\n')

		# INITIAL PARAMETER VALUES
		if(fix.start.end){
			init_param_values <- rep(m[1, 1], min.control.points - 2)
		}else{
			init_param_values <- rep(m[1, 1], min.control.points)			
		}
		
		# VECTOR FOR SAVING RESIDUAL STANDARD ERROR
		rse <- rep(NA, max.control.points - 2)
		
		# DEFAULT REASON
		fit.stopped.by[i] <- paste0('Maximum number of control points (',  max.control.points, ') reached')

		# EXPONENTIAL FIT START PARAMETERS
		exp_start <- list(b1=1, b2=-0.5, b3=0.1)

		# TRY FIT WITH INCREASING NUMBER OF PARAMETERS
		for(j in (min.control.points - 2):(max.control.points - 2)){
			
			# MAKE INITIAL PARAMETER VECTOR
			init_param <- init_param_values

			# RUN NON-LINEAR MODEL FIT ON BEZIER CONTROL POINTS (CONSTRAINING START AND END POINTS)
			if(fix.start.end){
				model_r <- nls(m[, i] ~ bezier(t=t, p=p, start=m[1, i], end=m[nrow(m), i]), start = list(p = init_param), trace = F, control = nls.control(maxiter = 100, warnOnly = TRUE, minFactor = 1/2048))	#nls.control(maxiter = 100, warnOnly = TRUE)
			}else{
				model_r <- nls(m[, i] ~ bezier(t=t, p=p), start = list(p = init_param), trace = F, control = nls.control(maxiter = 100, warnOnly = TRUE, minFactor = 1/2048))
			}

			# SAVE RESIDUAL STANDARD ERROR
			rse[j] <- summary(model_r)$sigma
			#cat('\t', j, ') RSE: ', rse[j], '\n', sep='')

			# SAVE PARAMETER ESTIMATES
			parameter_estimate <- summary(model_r)$parameters[, 'Estimate']

			# STOP IF MINIMUM RSE IS REACHED
			if(!is.null(max.rse) && rse[j] < max.rse){fit.stopped.by[i] <- paste0('RSE (', round(rse[j], 5), ') is less than maximum RSE (',  max.rse, ')');break}

			# USE PARAMETER ESTIMATES IN SUBSEQUENT MODEL FIT
			init_param_values <- c(parameter_estimate, init_param_values[1])

			# TEST IF MAX.RSE.PERCENT.CHANGE HAS BEEN REACHED USING A LINEAR REGRESSION EXCLUDING FIRST POINT
			if(!is.null(max.rse.percent.change) && sum(!is.na(rse[2:length(rse)])) >= 3 && sum(!is.na(rse[2:length(rse)])) <= 6){

				# FIND SLOPE OF LAST TEN POINTS
				rse_deriv <- summary(lm(rse ~ num_param, data=data.frame(num_param=1:length(na.omit(rse[2:length(rse)])), rse=na.omit(rse[2:length(rse)]))))$coefficients[2, 1]

				# STOP IF MINIMUM PERCENT CHANGE IN RSE IS REACHED
				#print(abs(rse_deriv/rse[j]))
				if(abs(rse_deriv/rse[j]) < max.rse.percent.change){fit.stopped.by[i] <- paste0('Percent change in RSE (', round(abs(rse_deriv/rse[j]), 5), ') is less than maximum percent change in RSE (',  max.rse.percent.change, ')');break}
			}

			# TEST IF MAX.RSE.PERCENT.CHANGE HAS BEEN REACHED USING AN EXPONENTIAL REGRESSION
			if(!is.null(max.rse.percent.change) && sum(!is.na(rse)) > 7){
			
				# FIT AN EXPONENTIAL FUNCTION TO RSE VALUES
				model <- nls(rse ~ b1*exp(b2*num_param) + b3, data=data.frame(num_param=1:length(na.omit(rse)), rse=na.omit(rse)), start=exp_start, control = nls.control(maxiter = 100, warnOnly = TRUE, minFactor = 1/2048))
				#print(data.frame(num_param=1:length(na.omit(rse)), rse=na.omit(rse)))

				# GET PARAMETER ESTIMATES
				b <- summary(model)$parameters[, 'Estimate']
				#print(b)

				# SET EXPONENTIAL FIT PARAMETERS AS START PARAMETERS FOR NEXT FIT
				exp_start <- list(b1=b[1], b2=b[2], b3=b[3])

				# GET DERIVATIVE OF FUNCTION AT CURRENT ITERATION
				rse_deriv <- b[2]*b[1]*exp(b[2]*j)

				#cat('\t', 'Percent rate of change in RSE: ', rse_deriv/rse[j], '\n', sep='')

				# STOP IF MINIMUM PERCENT CHANGE IN RSE IS REACHED
				#print(abs(rse_deriv/rse[j]))
				if(abs(rse_deriv/rse[j]) < max.rse.percent.change){fit.stopped.by[i] <- paste0('Percent change in RSE (', round(abs(rse_deriv/rse[j]), 5), ') is less than maximum percent change in RSE (',  max.rse.percent.change, ')');break}
			}			
		}
		
		# SAVE VALUES TO OUTPUT VECTORS
		if(fix.start.end){
			p[[i]] <- c(m[1, i], parameter_estimate, m[nrow(m), i])
		}else{
			p[[i]] <- parameter_estimate
		}
		p.rse[i] <- rse[j]
	}
	
	l <- list(p=p, rse=p.rse, fit.stopped.by=fit.stopped.by)
	class(l) <- 'bezierCurveFit'
	l
}
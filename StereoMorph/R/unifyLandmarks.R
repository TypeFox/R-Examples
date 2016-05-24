unifyLandmarks <- function(lm.array, min.common = dim(lm.array)[2], return.on.error = FALSE){
	# Modified from R function unifyVD() written by Annat Haber
	# Rohlf 1990. Proceedings of the Michigan Morphometrics Workshop

	# CREATE ARRAY WITH ONE EXTRA DIMENSION FOR UNIFIED MATRICES
	lm_array <- array(NA, dim=c(dim(lm.array)[1], dim(lm.array)[2], dim(lm.array)[3]+1), dimnames=list(dimnames(lm.array)[[1]], NULL, NULL))

	# COPY INPUT LANDMARK ARRAY INTO ARRAY
	lm_array[, , 1:dim(lm.array)[3]] <- lm.array

	# FIND NUMBER OF SETS
	num_sets <- dim(lm.array)[3]

	# RETURN ARRAY AS MATRIX IF ONLY ONE SET
	if(num_sets == 1){	
		l <- list(lm.matrix = lm_array[, , 1], unify.seq = 1, unify.error = NA, unify.rmse = NA)
		class(l) <- 'unifyLandmarks'
		return(l)
	}

	# UNIFICATION ERROR BY LANDMARK MATRIX
	unify.error <- matrix(NA, nrow=dim(lm.array)[1], ncol=num_sets-1, dimnames=list(dimnames(lm.array)[[1]]))

	# UNIFICATION SEQUENCE
	unify.seq <- rep(NA, num_sets)

	# INITIAL PAIR SEQUENCE ELEMENTS
	n1_include <- n2_include <- c(1:num_sets, NA)

	for(i in 1:(num_sets - 1)){

		# GENERATE PAIRS OF MATRICES TO TEST FOR UNIFICATION ERROR
		pairs <- matrix(NA, nrow=0, ncol=2)		
		for(j in na.omit(n1_include)){
			for(k in na.omit(n2_include)){
	
				# CHECK IF SAME NUMBERS IN REVERSE ARE ALREADY PRESENT IN MATRIX, IF SPECIFIED
				if(k %in% pairs[, 1]) if(j %in% pairs[k == pairs[, 1], 2]) next
	
				# SKIP IDENTICAL NUMBERS, IF SPECIFIED
				if(k == j) next
	
				# ADD PAIRS TO MATRIX
				pairs <- rbind(pairs, c(j, k))
			}
		}

		# CREATE ARRAY TO SAVE EACH UNIFIED MATRIX COMBINATION
		unify_array <- array(NA, dim=c(dim(lm.array)[1], dim(lm.array)[2], nrow(pairs)), dimnames=list(dimnames(lm.array)[[1]], dimnames(lm.array)[[2]], NULL))

		# UNIFICATION ROOT MEAN SQUARE ERROR VECTOR
		error_matrix <- matrix(NA, nrow=dim(lm.array)[1], ncol=nrow(pairs), dimnames=list(dimnames(lm.array)[[1]], NULL))

		# UNIFICATION ROOT MEAN SQUARE ERROR VECTOR
		rmse_error <- rep(NA, nrow(pairs))
		
		for(j in 1:nrow(pairs)){			

			# FILL INITIAL LANDMARK MATRICES
			all_m1 <- common_m1 <- lm_array[, , pairs[j, 1]]
			all_m2 <- lm_array[, , pairs[j, 2]]

			# REPLACE NON-COMMON LANDMARKS BETWEEN TWO MATRICES WITH NA
			common_m1[which(is.na(all_m2))] <- NA
			
			# SKIP IF NUMBER OF COMMON LANDMARKS IS LESS THAN MIN.COMMON
			if(sum(!is.na(common_m1[, 1])) < min.common) next

			# TRANSLATE AND ROTATE M2 TO MINIMIZE ALIGNMENT ERROR BETWEEN M1 AND M2
			all_m2 <- findOptimalPointAlignment(all_m1, all_m2)

			# GET SUMMED DIFFERENCES BETWEEN COMMON POINT SETS AFTER UNIFICATION
			error_matrix[, j] <- rowSums((all_m1 - all_m2)^2)

			# SAVE ROOT MEAN SQUARE ERROR FOR UNIFICATION OF EACH SET PAIR
			rmse_error[j] <- sqrt(sum(na.omit(error_matrix[, j]))/length(na.omit(error_matrix[, j])))

			# TRANSFER LANDMARKS MISSING FROM ONE MATRIX TO THE OTHER
			all_m1[which(is.na(all_m1))] <- all_m2[which(is.na(all_m1))]
			all_m2[which(is.na(all_m2))] <- all_m1[which(is.na(all_m2))]

			# AVERAGE MATRICES AND SAVE TO UNIFICATION ARRAY
			unify_array[, , j] <- (all_m1 + all_m2)/2
		}

		# NOT ENOUGH POINTS TO UNIFY SETS
		if(sum(!is.na(rmse_error)) == 0){
			
			if(return.on.error) return(NULL)

			stop(paste0("The number of common points between the two aspects to unify is less than the specified minimum (", min.common, ")."))
		}

		# FIND MINIMUM UNIFICATION ERROR
		min_error <- which.min(rmse_error)

		# SAVE ERROR BETWEEN SETS FOR PARTICULAR LANDMARKS
		unify.error[, i] <- error_matrix[, min_error]

		# SAVE UNIFICATION SEQUENCE
		if(i == 1){unify.seq[1:2] <- pairs[min_error, ]}else{unify.seq[i+1] <- pairs[min_error, 2]}
		
		# APPEND CORRESPONDING MATRIX TO ARRAY
		lm_array[, , num_sets+1] <- unify_array[, , min_error]

		# BREAK WHEN ONLY ONE MATRIX REMAINS
		if(dim(unify_array)[3] == 1) break

		# REMOVE ALREADY ALIGNED MATRICES FROM POTENTIAL FUTURE PAIRS
		n2_include[c(pairs[min_error, 1], pairs[min_error, 2])] <- NA

		# AFTER FIRST RUN SET FIRST MATRIX TO FIRST UNIFIED MATRIX APPENDED TO THE END OF THE LANDMARK ARRAY
		if(i == 1) n1_include <- num_sets+1
	}

	# SAVE TO LANDMARK MATRIX WITH PROPER DIMENSION NAMES
	lm.matrix <- matrix(lm_array[, , num_sets+1], nrow=dim(lm_array)[1], ncol=dim(lm_array)[2], dimnames=list(dimnames(lm_array)[[1]], c('x', 'y', 'z')[1:dim(lm_array)[2]]))

	l <- list(
		lm.matrix=lm.matrix,
		unify.seq = unify.seq,
		unify.error = unify.error,
		unify.rmse = sqrt(colSums(unify.error, na.rm = TRUE) / colSums(!is.na(unify.error))))
	class(l) <- 'unifyLandmarks'
	l
}

summary.unifyLandmarks <- function(object, print.tab = '', verbose = TRUE, ...){

	r <- '\n'

	r <- c(r, print.tab, 'unifyLandmarks Summary\n')
	if(verbose) r <- c(r, print.tab, '\tUnification sequence: ', paste0(object$unify.seq, collapse=', '))
	if(verbose) r <- c(r, '\n')

	r <- c(r, print.tab, '\tUnification RMS Error:')
	if(is.na(object$unify.rmse[1])){
		r <- c(r, print.tab, ' NA')
	}else{
		#r <- c(r, '\n\t\t')
		#r <- c(r, print.tab, paste(format(object$unify.rmse), collapse=paste0('\n\t\t', print.tab)))
		r <- c(r, paste0(' ', paste(format(object$unify.rmse), collapse=', ')))
	}
	r <- c(r, '\n')

	if(verbose){
		r <- c(r, print.tab, '\tUnification landmark errors by sequence:')
		if(!is.matrix(object$unify.error)){
			r <- c(r, print.tab, '\tNA\n')
		}else{
			r <- c(r, '\n')
			for(i in 1:ncol(object$unify.error)){
				r <- c(r, print.tab, '\t\t[[', i, ']]\n')
			
				m <- as.matrix(object$unify.error[!is.na(object$unify.error[, i]), i])
				maxnchar <- max(nchar(rownames(m))) + 3
			
				for(j in 1:nrow(m)){
					r <- c(r, print.tab, paste0('\t\t\t', rownames(m)[j], paste(rep(' ', length=maxnchar - nchar(rownames(m)[j])), collapse=''), format(m[j, ]), sep='\n'))
				}
				r <- c(r, '\n')
			}
		}
	}

	#r <- c(r, '\n')
	class(r) <- "summary.unifyLandmarks"
	r
}

print.summary.unifyLandmarks <- function(x, ...) cat(x, sep='')
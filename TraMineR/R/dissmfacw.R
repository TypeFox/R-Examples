###Based on the program written for scipy (Python) by
###Ondrej Libiger and Matt Zapala
###
###Based on some of the methods presented in:
###McArdle, B.H. and Anderson, M.J. (2001).
###Fitting multivariate models to community data:
###a comment on distance-based redundancy analysis. Ecology, 82, 290-297.



gower_matrix <- function(diss, squared=TRUE, weights=NULL) {
	n <- nrow(diss)
	if (squared) {
		diss <- diss^2
	}
	One <- rep(1, n)
	if(is.null(weights)){
		weights <- One
	}
	s <- weights/sum(weights)
	
	return((diag(n) - tcrossprod(One,s)) %*%(diss/(-2))%*% (diag(n) - tcrossprod(s,One)))
	# Creating centered (Gower's) matrix from the distance matrix:
	
}

dissmfacw <- function(formula, data, R=1000, gower=FALSE,
					squared=FALSE, weights=NULL) {
	
		#To ensure dissmatrix can be a dist object as well
		dissmatrix <- as.matrix(eval(formula[[2]], data, parent.frame()))
		formula[[2]] <- NULL
		#Model matrix from forumla
		predictor_frame <- model.frame(formula, data, na.action=na.pass, drop.unused.levels = TRUE)
		seqok <- complete.cases(predictor_frame)
		#To make sure unused levels after removing NA values are actually removed
		predictor_frame <- model.frame(formula, data[seqok, ], drop.unused.levels = TRUE)
		
		predictor_terms <- attr(predictor_frame, "terms")
		var_complete_name <- c(attr(predictor_terms, "term.labels"), "Total")
		predictor <- model.matrix(formula, predictor_frame)
		if (is.null(weights)) {
			weights <- rep(1, sum(seqok))
		} else {
			weights <- weights[seqok]
		}
		
		dissmatrix <- dissmatrix[seqok, seqok]

		if (!gower) {
			g_matrix <- gower_matrix(dissmatrix, squared, weights=weights)

		} else {
			g_matrix <- dissmatrix
		}
		# n <- nrow(g_matrix)
		totalweight <- sum(weights)
		var_list <- attr(predictor, "assign")
		W_sqrt <- sqrt(weights)
		W_sqrt_mat <- tcrossprod(W_sqrt)
		# hat_matrix_qr <- function(pred) {
			# qr_matrix <- qr(W_sqrt*pred)
			# q_matrix <- qr.Q(qr_matrix)
			# hat_matrix <- tcrossprod(q_matrix)
		# }
		hatw_matrix_qr <- function(pred) {
			qr_matrix <- qr(W_sqrt*pred)
			q_matrix <- qr.Q(qr_matrix)
			hat_matrix <- tcrossprod(q_matrix)
			return(W_sqrt_mat*hat_matrix)
		}
		var_list_index <- (1:length(var_list))
		#List of variable
		var_names <- unique(var_list)
		#Number of variable (minus cte)
		nterms <- length(var_names) - 1
		SCtot <- sum(weights * diag(g_matrix))
		
		p_list <- numeric(nterms)
		SSbv <- numeric(nterms)
		#Compute all  "backward" SCexp based on QR decomposition
		for (var in 1:(nterms)) {
			pred <- predictor[, c(var_list_index[var_list!=var])]
			p_list[var] <- length(var_list_index[var_list==var])
			hwm <- hatw_matrix_qr(pred)
			##No need to transpose, G is symmetric
			SSbv[var]<- sum(hwm*g_matrix)
		}

		#complete_model <- nterms+1
		#If we permute dissmatrix, we can build all hat_matrix in advance
		internalbootmatrixregression2 <- function(predictor, ind) {
			perm_predictor <- predictor[ind, ]
			complete_model <- nterms+1
			R2_list <- numeric(complete_model)
			F_list <- numeric(complete_model)
			m <- ncol(perm_predictor)
			hatw_matrix <- hatw_matrix_qr(perm_predictor)
			SCexpC <- sum(hatw_matrix * g_matrix)
			SCresC <- SCtot-SCexpC
			R2_list[complete_model] <- (SCexpC/SCtot)
			F_list[complete_model] <- (SCexpC/(m-1)) / ((SCresC)/(totalweight-m))
			if (nterms==1) {
				R2_list[1] <- R2_list[complete_model]
				F_list[1] <- F_list[complete_model]
				return(c(F_list, R2_list))
			}
			for (var in 1:nterms) {
				pred <- predictor 
				pred[, c(var_list_index[var_list==var])] <- perm_predictor[, c(var_list_index[var_list==var])]
				hwm <- hatw_matrix_qr(pred)
				SCexp <- sum(hwm * g_matrix)
				SCres <- SCtot-SCexp
				F_list[[var]] <- ((SCexp-SSbv[var])/(p_list[var])) / ((SCres)/(totalweight-m))
				#calculate proportion of variance explained
				R2_list[[var]] <- (SCexp-SSbv[var]) / SCtot
				# Printing single regression analysis output
			}
			return(c(F_list, R2_list))
		}
		
			#If we permute the models, then hat_matrix is in function...
			gc()
			bts <- TraMineR.permutation(predictor, R=R, statistic=internalbootmatrixregression2)
		
		#Computing Pvalue based on permutations tests
		
		#Results
		mfac <- data.frame(Variable=var_complete_name,
				PseudoF=bts$t0[1:(nterms+1)],
				PseudoR2=bts$t0[(nterms+2):length(bts$t0)],
				p_value=bts$pval[1:(nterms+1)]
			)
		#Var names
		dmf <- list(mfac = mfac, call = match.call(), perms=bts)
		class(dmf) <- "dissmultifactor"
		return(dmf)

}
print.dissmultifactor <- function(x, pvalue.confint=0.95, digits = NULL, ...) {
	mfac <- x$mfac
	if(is.null(digits)) {
		digits = getOption("digits")
	}
	print(mfac,digits=digits, ...)
	if(!is.null(pvalue.confint)){
		## Two sided
		cat("\nInconclusive intervals: \n")
		pvalue.confint <- (1-(1-pvalue.confint)/2)
		pvalues <- c(0.01, 0.05)
		confinter <- qnorm(pvalue.confint)*sqrt(pvalues*(1-pvalues)/x$perms$R)
		cat(paste(format(pmax(0,(pvalues-confinter)), digits=3), " < ",pvalues, " < ",format(pmin(1,(pvalues+confinter)),digits=3), collapse="\n"), "\n")
	}
	invisible(x)
}
print.dissregression <- function(x, ...) {
	warning("dissreg function is deprecated. It has been renamed dissmfac.")
	print(x$mreg)
	invisible(x)
}


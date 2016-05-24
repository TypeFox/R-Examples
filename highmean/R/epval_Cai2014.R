epval_Cai2014 <- function(sam1, sam2, eq.cov = TRUE, n.iter = 1000, cov1.est, cov2.est, bandwidth1, bandwidth2, cv.fold = 5, norm = "F", seeds){
	if(missing(seeds)) seeds <- NULL
	if(length(seeds) != n.iter){
		seeds <- NULL
		cat("The length of seeds does not match the specified n.iter.\n")
		cat("Seeds for each permutation/resampling iteration are assigned randomly.\n")
	}
	if(eq.cov){
		out <- epval_Cai2014_samecov(sam1, sam2, n.iter, seeds)
	}else{
		sam.cov1 <- cov(sam1)
		sam.cov2 <- cov(sam2)
		p <- dim(sam1)[2]

		if(cv.fold != floor(cv.fold) | cv.fold <= 1) stop("cv.fold must be an integer greater than 1.")
		if(missing(bandwidth1)) bandwidth1 <- seq(from = 0, to = p, by = floor(p/50))
		if(missing(bandwidth2)) bandwidth2 <- seq(from = 0, to = p, by = floor(p/50))
		if(any(bandwidth1 < 0)){
			cat("Negative values specified in bandwidth1 are removed.\n")
			bandwidth1 <- bandwidth1[bandwidth1 < 0]
		}
		if(any(bandwidth2 < 0)){
			cat("Negative values specified in bandwidth2 are removed.\n")
			bandwidth2 <- bandwidth2[bandwidth2 < 0]
		}
		if(any(bandwidth1 != floor(bandwidth1))){
			cat("Non-integers specified in bandwidth1 are converted to their integer parts.")
			bandwidth1 <- floor(bandwidth1)
		}
		if(any(bandwidth2 != floor(bandwidth2))){
			cat("Non-integers specified in bandwidth2 are converted to their integer parts.")
			bandwidth2 <- floor(bandwidth2)
		}

		if(missing(cov1.est)){
			output.opt.bw1 <- TRUE
			if(length(bandwidth1) > 1){
				optim.bandwidth1 <- best.band(sam1, bandwidth1, cv.fold, norm)
			}
			if(length(bandwidth1) == 1){
				optim.bandwidth1 <- bandwidth1
			}
			if(optim.bandwidth1 > 0){
				cov1.est <- sam.cov1
				cov1.est[abs(row(cov1.est) - col(cov1.est)) > optim.bandwidth1] <- 0
			}
			if(optim.bandwidth1 == 0){
				cov1.est <- diag(diag(sam.cov1))
			}
			eigen1 <- eigen(cov1.est)
			eigen1.vectors <- eigen1$vectors
			eigen1.values <- eigen1$values
			eigen1.values[eigen1.values <= 0] <- 0.001
			cov1.est <- eigen1.vectors %*% diag(eigen1.values) %*% t(eigen1.vectors)
		}else{
			output.opt.bw1 <- FALSE
		}

		if(missing(cov2.est)){
			output.opt.bw2 <- TRUE
			if(length(bandwidth2) > 1){
				optim.bandwidth2 <- best.band(sam2, bandwidth2, cv.fold, norm)
			}
			if(length(bandwidth2) == 1){
				optim.bandwidth2 <- bandwidth2
			}
			if(optim.bandwidth2 > 0){
				cov2.est <- sam.cov2
				cov2.est[abs(row(cov2.est) - col(cov2.est)) > optim.bandwidth2] <- 0
			}
			if(optim.bandwidth2 == 0){
				cov2.est <- diag(diag(sam.cov2))
			}
			eigen2 <- eigen(cov2.est)
			eigen2.vectors <- eigen2$vectors
			eigen2.values <- eigen2$values
			eigen2.values[eigen2.values <= 0] <- 0.001
			cov2.est <- eigen2.vectors %*% diag(eigen2.values) %*% t(eigen2.vectors)
		}else{
			output.opt.bw2 <- FALSE
		}

		out <- epval_Cai2014_diffcov(sam1, sam2, n.iter, sam.cov1, sam.cov2, cov1.est, cov2.est, cv.fold, norm, seeds, optim.bandwidth1, optim.bandwidth2, output.opt.bw1, output.opt.bw2)
	}

	return(out)
}
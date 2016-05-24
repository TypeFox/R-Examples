apval_aSPU <- function(sam1, sam2, pow = c(1:6, Inf), eq.cov = TRUE, cov.est, cov1.est, cov2.est, bandwidth, bandwidth1, bandwidth2, cv.fold = 5, norm = "F"){
	if(eq.cov){
		n1 <- dim(sam1)[1]
		n2 <- dim(sam2)[1]
		p <- dim(sam1)[2]
		sam.cov <- ((n1 - 1)*cov(sam1) + (n2 - 1)*cov(sam2))/(n1 + n2 - 2)

		if(missing(bandwidth)) bandwidth <- seq(from = 0, to = p, by = floor(p/50))
		if(any(bandwidth < 0)){
			cat("Negative values specified in bandwidth are removed.\n")
			bandwidth <- bandwidth[bandwidth < 0]
		}
		if(any(bandwidth != floor(bandwidth))){
			cat("Non-integers specified in bandwidth are converted to their integer parts.")
			bandwidth <- floor(bandwidth)
		}

		if(missing(cov.est)){
			output.opt.bw <- TRUE
			if(length(bandwidth) > 1){
				optim.bandwidth <- best.band(rbind(sam1 - matrix(colMeans(sam1), byrow = TRUE, nrow = n1, ncol = p), sam2 - matrix(colMeans(sam2), byrow = TRUE, nrow = n2, ncol = p)), bandwidth, cv.fold, norm)
			}
			if(length(bandwidth) == 1){
				optim.bandwidth <- bandwidth
			}
			if(optim.bandwidth > 0){
				cov.est <- sam.cov
				cov.est[abs(row(cov.est) - col(cov.est)) > optim.bandwidth] <- 0
			}
			if(optim.bandwidth == 0){
				cov.est <- diag(diag(sam.cov))
			}
		}else{
			output.opt.bw <- FALSE
		}

		out <- apval_aSPU_samecov(sam1, sam2, pow, sam.cov, cov.est, cv.fold, norm, optim.bandwidth, output.opt.bw)
	}else{
		sam.cov1 <- cov(sam1)
		sam.cov2 <- cov(sam2)
		p <- dim(sam1)[2]

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
		}else{
			output.opt.bw2 <- FALSE
		}

		out <- apval_aSPU_diffcov(sam1, sam2, pow, sam.cov1, sam.cov2, cov1.est, cov2.est, cv.fold, norm, optim.bandwidth1, optim.bandwidth2, output.opt.bw1, output.opt.bw2)
	}

	return(out)
}
logistic.boot <- function(part1.list, seed = new.seed(), alpha = 0.05, paired = FALSE, N.iter = 1000, cores = 1, p.adj = "oleson", test.spots = NULL, time.test = NULL, test.params = FALSE) {
  data       <- part1.list$data
  col        <- part1.list$col
  rho.0      <- part1.list$rho.0
  N.time     <- part1.list$N.time
  coef.id1   <- part1.list$coef.id1
  coef.id2   <- part1.list$coef.id2
  coef.id3   <- part1.list$coef.id3
  coef.id4   <- part1.list$coef.id4
  sdev.id1   <- part1.list$sdev.id1
  sdev.id2   <- part1.list$sdev.id2
  sdev.id3   <- part1.list$sdev.id3
  sdev.id4   <- part1.list$sdev.id4
  sigma.id1  <- part1.list$sigma.id1
  sigma.id2  <- part1.list$sigma.id2
  sigma.id3  <- part1.list$sigma.id3
  sigma.id4  <- part1.list$sigma.id4
  id.nums.g1 <- part1.list$id.nums.g1
  id.nums.g2 <- part1.list$id.nums.g2
  groups     <- part1.list$groups
  time.all   <- part1.list$time.all
  N.g1       <- part1.list$N.g1
  N.g2       <- part1.list$N.g2
  diffs      <- part1.list$diffs
	
	if(!is.null(test.spots)) time.all <- test.spots
	N.tests <- length(time.all)
	N.time <- length(time.all)
	
	if(diffs) {
		group1.bad <- is.na(coef.id1[,1]) | is.na(coef.id3[,1])
		group2.bad <- is.na(coef.id2[,1]) | is.na(coef.id4[,1])
	} else {
		group1.bad <- is.na(coef.id1[,1])
		group2.bad <- is.na(coef.id2[,1])
	}

	
	coef.id1 <- subset(coef.id1, !group1.bad)
	coef.id3 <- subset(coef.id3, !group1.bad)
	coef.id2 <- subset(coef.id2, !group2.bad)
	coef.id4 <- subset(coef.id4, !group2.bad)
	
	sdev.id1 <- subset(sdev.id1, !group1.bad)
	sdev.id3 <- subset(sdev.id3, !group1.bad)
	sdev.id2 <- subset(sdev.id2, !group2.bad)
	sdev.id4 <- subset(sdev.id4, !group2.bad)
	
	sigma.id1 <- subset(sigma.id1, !group1.bad)
	sigma.id3 <- subset(sigma.id3, !group1.bad)
	sigma.id2 <- subset(sigma.id2, !group2.bad)
	sigma.id4 <- subset(sigma.id4, !group2.bad)
	
	id.nums.g1 <- id.nums.g1[!group1.bad]
	id.nums.g2 <- id.nums.g2[!group2.bad]
	N.g1 <- N.g1 - sum(group1.bad)
	N.g2 <- N.g2 - sum(group2.bad)
  
  curve1.0 <-   matrix(NA, ncol = N.time, nrow = N.g1)
	mini3.ran  <- mini1.ran  <- rep(NA, N.g1)
	peak3.ran  <- peak1.ran  <- rep(NA, N.g1)
	slope3.ran <- slope1.ran <- rep(NA, N.g1)
	cross3.ran <- cross1.ran <- rep(NA, N.g1)
	
  curve2.0 <-   matrix(NA, ncol = N.time, nrow = N.g2)
	mini4.ran  <- mini2.ran  <- rep(NA, N.g2)
	peak4.ran  <- peak2.ran  <- rep(NA, N.g2)
	slope4.ran <- slope2.ran <- rep(NA, N.g2)
	cross4.ran <- cross2.ran <- rep(NA, N.g2)
  
	curve1.mat <- matrix(NA,ncol=N.time,nrow=N.iter)
	curve2.mat <- matrix(NA,ncol=N.time,nrow=N.iter)
	curve3.mat <- matrix(NA,ncol=N.time,nrow=N.iter)
  
	if(paired && N.g1 != N.g2) stop("Can't do paired test for different group sizes")
	
	#Target Curve
	curve.f <- function(mini,peak,slope,cross,t)
		mini + (peak - mini) / (1 + exp(4 * slope * (cross - t) / (peak-mini)))
		
	mini.1 <- peak.1 <- slope.1 <- cross.1 <-
		mini.2 <- peak.2 <- slope.2 <- cross.2 <- numeric(N.iter)
	
	##################
	##### 1 Core #####
	##################
	
	if(paired) {
		mini.cov.1 <- cov(coef.id1[,1], coef.id2[,1], use = "pairwise.complete.obs")
		peak.cov.1 <- cov(coef.id1[,2], coef.id2[,2], use = "pairwise.complete.obs")
		slope.cov.1 <- cov(coef.id1[,3], coef.id2[,3], use = "pairwise.complete.obs")
		cross.cov.1 <- cov(coef.id1[,4], coef.id2[,4], use = "pairwise.complete.obs")
		
		if(diffs) {
			mini.cov.2 <- cov(coef.id3[,1], coef.id4[,1], use = "pairwise.complete.obs")
			peak.cov.2 <- cov(coef.id3[,2], coef.id4[,2], use = "pairwise.complete.obs")
			slope.cov.2 <- cov(coef.id3[,3], coef.id4[,3], use = "pairwise.complete.obs")
			cross.cov.2 <- cov(coef.id3[,4], coef.id4[,4], use = "pairwise.complete.obs")
		}
	}
	
	if(cores == 1) {
		set.seed(seed)
		for(iter in 1:N.iter){
			if(paired) {
				for(i in 1:N.g1) {
					mini <- rmvnorm(1, mean = c(coef.id1[i,1], coef.id2[i,1]),
						sigma = matrix(c(sdev.id1[i,1] ^ 2, mini.cov.1, mini.cov.1, sdev.id2[i,1] ^ 2), nrow = 2))
					peak <- rmvnorm(1, mean = c(coef.id1[i,2], coef.id2[i,2]),
						sigma = matrix(c(sdev.id1[i,2] ^ 2, peak.cov.1, peak.cov.1, sdev.id2[i,2] ^ 2), nrow = 2))
					slope <- rmvnorm(1, mean = c(coef.id1[i,3], coef.id2[i,3]),
						sigma = matrix(c(sdev.id1[i,3] ^ 2, slope.cov.1, slope.cov.1, sdev.id2[i,3] ^ 2), nrow = 2))
					cross <- rmvnorm(1, mean = c(coef.id1[i,4], coef.id2[i,4]),
						sigma = matrix(c(sdev.id1[i,4] ^ 2, cross.cov.1, cross.cov.1, sdev.id2[i,4] ^ 2), nrow = 2))
						
					mini1.ran <- mini[1]; mini2.ran <- mini[2]
					peak1.ran <- peak[1]; peak2.ran <- peak[2]
					slope1.ran <- slope[1]; slope2.ran <- slope[2]
					cross1.ran <- cross[1]; cross2.ran <- cross[2]
				}
			} else {
				mini1.ran <-  rnorm(N.g1, coef.id1[,1], sdev.id1[,1])
				peak1.ran <-  rnorm(N.g1, coef.id1[,2], sdev.id1[,2])
				slope1.ran <- rnorm(N.g1, coef.id1[,3], sdev.id1[,3])
				cross1.ran <- rnorm(N.g1, coef.id1[,4], sdev.id1[,4])

				mini2.ran <-  rnorm(N.g2, coef.id2[,1], sdev.id2[,1])
				peak2.ran <-  rnorm(N.g2, coef.id2[,2], sdev.id2[,2])
				slope2.ran <- rnorm(N.g2, coef.id2[,3], sdev.id2[,3])
				cross2.ran <- rnorm(N.g2, coef.id2[,4], sdev.id2[,4])
			}
			
			mini.1[iter] <- mean(mini1.ran); mini.2[iter] <- mean(mini2.ran)
			peak.1[iter] <- mean(peak1.ran); peak.2[iter] <- mean(peak2.ran)
			slope.1[iter] <- mean(slope1.ran); slope.2[iter] <- mean(slope2.ran)
			cross.1[iter] <- mean(cross1.ran); cross.2[iter] <- mean(cross2.ran)
			
			if(diffs) {
				if(paired) {
					for(i in 1:N.g1) {
						mini <- rmvnorm(1, mean = c(coef.id3[i,1], coef.id4[i,1]),
							sigma = matrix(c(sdev.id3[i,1] ^ 2, mini.cov.1, mini.cov.1, sdev.id4[i,1] ^ 2), nrow = 2))
						peak <- rmvnorm(1, mean = c(coef.id3[i,2], coef.id4[i,2]),
							sigma = matrix(c(sdev.id3[i,2] ^ 2, peak.cov.1, peak.cov.1, sdev.id4[i,2] ^ 2), nrow = 2))
						slope <- rmvnorm(1, mean = c(coef.id3[i,3], coef.id4[i,3]),
							sigma = matrix(c(sdev.id3[i,3] ^ 2, slope.cov.1, slope.cov.1, sdev.id4[i,3] ^ 2), nrow = 2))
						cross <- rmvnorm(1, mean = c(coef.id3[i,4], coef.id4[i,4]),
							sigma = matrix(c(sdev.id3[i,4] ^ 2, cross.cov.1, cross.cov.1, sdev.id4[i,4] ^ 2), nrow = 2))
						
						mini3.ran <- mini[1]; mini4.ran <- mini[2]
						peak3.ran <- peak[1]; peak4.ran <- peak[2]
						slope3.ran <- slope[1]; slope4.ran <- slope[2]
						cross3.ran <- cross[1]; cross4.ran <- cross[2]
					}
				} else {
					mini3.ran <-  rnorm(N.g1, coef.id3[,1], sdev.id3[,1]) 
					peak3.ran <-  rnorm(N.g1, coef.id3[,2], sdev.id3[,2])
					slope3.ran <- rnorm(N.g1, coef.id3[,3], sdev.id3[,3])
					cross3.ran <- rnorm(N.g1, coef.id3[,4], sdev.id3[,4])

					mini4.ran <-  rnorm(N.g2, coef.id4[,1], sdev.id4[,1])
					peak4.ran <-  rnorm(N.g2, coef.id4[,2], sdev.id4[,2])
					slope4.ran <- rnorm(N.g2, coef.id4[,3], sdev.id4[,3])
					cross4.ran <- rnorm(N.g2, coef.id4[,4], sdev.id4[,4])
				}
			}
			
			if(diffs) {
				for(id in 1:N.g1){ #Get fixation level for each g1 subject
					curve1.0[id,] <- curve.f(mini1.ran[id], peak1.ran[id], slope1.ran[id],
																		cross1.ran[id], time.all) -
													 curve.f(mini3.ran[id], peak3.ran[id], slope3.ran[id],
																		cross3.ran[id], time.all)
				}
				for(id in 1:N.g2){ #Get fixation level for each g2 subject
					curve2.0[id,] <- curve.f(mini2.ran[id], peak2.ran[id], slope2.ran[id],
																		cross2.ran[id], time.all) - 
													 curve.f(mini4.ran[id], peak4.ran[id], slope4.ran[id],
																		cross4.ran[id], time.all)
				}
			} else {
				for(id in 1:N.g1){ #Get fixation level for each g1 subject
					curve1.0[id,] <- curve.f(mini1.ran[id], peak1.ran[id], slope1.ran[id],
																		cross1.ran[id], time.all)
				}
				for(id in 1:N.g2){ #Get fixation level for each g2 subject
					curve2.0[id,] <- curve.f(mini2.ran[id], peak2.ran[id], slope2.ran[id],
																		cross2.ran[id], time.all)
				}
			}

			curve1.mat[iter,] <- apply(curve1.0, 2, mean) #Mean fixations at each time point for CIs
			curve2.mat[iter,] <- apply(curve2.0, 2, mean) #Mean fixations at each time point for NHs
			
			if(paired) curve3.mat[iter,] <- apply(curve2.0 - curve1.0, 2, mean)
		}
	} else {
		####################
		##### 2+ Cores #####
		####################
	
		cl <- makePSOCKcluster(cores)
		registerDoParallel(cl)
		
		for.out <- foreach(iter = 1:N.iter, .combine = rbind, .options.RNG = seed) %dorng% {
			if(paired) {
				for(i in 1:N.g1) {
					mini <- rmvnorm(1, mean = c(coef.id1[i,1], coef.id2[i,1]),
						sigma = matrix(c(sdev.id1[i,1] ^ 2, mini.cov.1, mini.cov.1, sdev.id2[i,1] ^ 2), nrow = 2))
					peak <- rmvnorm(1, mean = c(coef.id1[i,2], coef.id2[i,2]),
						sigma = matrix(c(sdev.id1[i,2] ^ 2, peak.cov.1, peak.cov.1, sdev.id2[i,2] ^ 2), nrow = 2))
					slope <- rmvnorm(1, mean = c(coef.id1[i,3], coef.id2[i,3]),
						sigma = matrix(c(sdev.id1[i,3] ^ 2, slope.cov.1, slope.cov.1, sdev.id2[i,3] ^ 2), nrow = 2))
					cross <- rmvnorm(1, mean = c(coef.id1[i,4], coef.id2[i,4]),
						sigma = matrix(c(sdev.id1[i,4] ^ 2, cross.cov.1, cross.cov.1, sdev.id2[i,4] ^ 2), nrow = 2))
						
					mini1.ran <- mini[1]; mini2.ran <- mini[2]
					peak1.ran <- peak[1]; peak2.ran <- peak[2]
					slope1.ran <- slope[1]; slope2.ran <- slope[2]
					cross1.ran <- cross[1]; cross2.ran <- cross[2]
				}
			} else {
				mini1.ran <-  rnorm(N.g1, coef.id1[,1], sdev.id1[,1])
				peak1.ran <-  rnorm(N.g1, coef.id1[,2], sdev.id1[,2])
				slope1.ran <- rnorm(N.g1, coef.id1[,3], sdev.id1[,3])
				cross1.ran <- rnorm(N.g1, coef.id1[,4], sdev.id1[,4])

				mini2.ran <-  rnorm(N.g2, coef.id2[,1], sdev.id2[,1])
				peak2.ran <-  rnorm(N.g2, coef.id2[,2], sdev.id2[,2])
				slope2.ran <- rnorm(N.g2, coef.id2[,3], sdev.id2[,3])
				cross2.ran <- rnorm(N.g2, coef.id2[,4], sdev.id2[,4])
			}
			
			mini.temp.1 <- mean(mini1.ran); mini.temp.2 <- mean(mini2.ran)
			peak.temp.1 <- mean(peak1.ran); peak.temp.2 <- mean(peak2.ran)
			slope.temp.1 <- mean(slope1.ran); slope.temp.2 <- mean(slope2.ran)
			cross.temp.1 <- mean(cross1.ran); cross.temp.2 <- mean(cross2.ran)
			
			if(diffs) {
				if(paired) {
					for(i in 1:N.g1) {
						mini <- rmvnorm(1, mean = c(coef.id3[i,1], coef.id4[i,1]),
							sigma = matrix(c(sdev.id3[i,1] ^ 2, mini.cov.1, mini.cov.1, sdev.id4[i,1] ^ 2), nrow = 2))
						peak <- rmvnorm(1, mean = c(coef.id3[i,2], coef.id4[i,2]),
							sigma = matrix(c(sdev.id3[i,2] ^ 2, peak.cov.1, peak.cov.1, sdev.id4[i,2] ^ 2), nrow = 2))
						slope <- rmvnorm(1, mean = c(coef.id3[i,3], coef.id4[i,3]),
							sigma = matrix(c(sdev.id3[i,3] ^ 2, slope.cov.1, slope.cov.1, sdev.id4[i,3] ^ 2), nrow = 2))
						cross <- rmvnorm(1, mean = c(coef.id3[i,4], coef.id4[i,4]),
							sigma = matrix(c(sdev.id3[i,4] ^ 2, cross.cov.1, cross.cov.1, sdev.id4[i,4] ^ 2), nrow = 2))
						
						mini3.ran <- mini[1]; mini4.ran <- mini[2]
						peak3.ran <- peak[1]; peak4.ran <- peak[2]
						slope3.ran <- slope[1]; slope4.ran <- slope[2]
						cross3.ran <- cross[1]; cross4.ran <- cross[2]
					}
				} else {
					mini3.ran <-  rnorm(N.g1, coef.id3[,1], sdev.id3[,1]) 
					peak3.ran <-  rnorm(N.g1, coef.id3[,2], sdev.id3[,2])
					slope3.ran <- rnorm(N.g1, coef.id3[,3], sdev.id3[,3])
					cross3.ran <- rnorm(N.g1, coef.id3[,4], sdev.id3[,4])

					mini4.ran <-  rnorm(N.g2, coef.id4[,1], sdev.id4[,1])
					peak4.ran <-  rnorm(N.g2, coef.id4[,2], sdev.id4[,2])
					slope4.ran <- rnorm(N.g2, coef.id4[,3], sdev.id4[,3])
					cross4.ran <- rnorm(N.g2, coef.id4[,4], sdev.id4[,4])
				}
			}
			
			if(diffs) {
				for(id in 1:N.g1){ #Get fixation level for each g1 subject
					curve1.0[id,] <- curve.f(mini1.ran[id], peak1.ran[id], slope1.ran[id],
																		cross1.ran[id], time.all) -
													 curve.f(mini3.ran[id], peak3.ran[id], slope3.ran[id],
																		cross3.ran[id], time.all)
				}
				for(id in 1:N.g2){ #Get fixation level for each g2 subject
					curve2.0[id,] <- curve.f(mini2.ran[id], peak2.ran[id], slope2.ran[id],
																		cross2.ran[id], time.all) - 
													 curve.f(mini4.ran[id], peak4.ran[id], slope4.ran[id],
																		cross4.ran[id], time.all)
				}
			} else {
				for(id in 1:N.g1){ #Get fixation level for each g1 subject
					curve1.0[id,] <- curve.f(mini1.ran[id], peak1.ran[id], slope1.ran[id],
																		cross1.ran[id], time.all)
				}
				for(id in 1:N.g2){ #Get fixation level for each g2 subject
					curve2.0[id,] <- curve.f(mini2.ran[id], peak2.ran[id], slope2.ran[id],
																		cross2.ran[id], time.all)
				}
			}

			curve1 <- apply(curve1.0, 2, mean) #Mean fixations at each time point for CIs
			curve2 <- apply(curve2.0, 2, mean) #Mean fixations at each time point for NHs
			curve3 <- curve2 - curve1
			
			c(curve1, curve2, curve3, mini.temp.1, peak.temp.1, slope.temp.1, cross.temp.1,
				mini.temp.2, peak.temp.2, slope.temp.2, cross.temp.2)
		}
		
		curve1.mat <- for.out[,1:N.time]
		curve2.mat <- for.out[,(N.time + 1):(2 * N.time)]
		curve3.mat <- for.out[,(2 * N.time + 1):(3 * N.time)]
		mini.1 <- for.out[, 3 * N.time + 1]
		peak.1 <- for.out[, 3 * N.time + 2]
		slope.1 <- for.out[, 3 * N.time + 3]
		cross.1 <- for.out[, 3 * N.time + 4]
		mini.2 <- for.out[, 3 * N.time + 5]
		peak.2 <- for.out[, 3 * N.time + 6]
		slope.2 <- for.out[, 3 * N.time + 7]
		cross.2 <- for.out[, 3 * N.time + 8]
		
		stopCluster(cl)
	}
	
	curve.mean1 <- apply(curve1.mat, 2, mean)
	curve.mean2 <- apply(curve2.mat, 2, mean)
	
	curve.g1    <- curve.mean1
	curve.g2    <- curve.mean2
	
	curve.sd1   <- apply(curve1.mat, 2, sd)
	curve.sd2   <- apply(curve2.mat, 2, sd)

	if(paired) {
		diff.mean <- apply(curve3.mat, 2, mean)
		curve.sd <- apply(curve3.mat, 2, sd)
		
		t.val <- diff.mean / curve.sd
		p.values <- 2 * (1 - pt(abs(t.val), N.g1 - 1))
	} else {
		t.num <- (curve.mean1 - curve.mean2)
		t.den <- sqrt((N.g1 * (N.g1 - 1) * curve.sd1 ^ 2 + N.g2 * (N.g2 - 1) * curve.sd2 ^ 2) / (N.g1 + N.g2 - 2) * (1 / N.g1 + 1 / N.g2))
		t.val <- t.num / t.den
		p.values <- 2 * (1 - pt(abs(t.val), (N.g1 + N.g2 - 2)))
	}

	# compute t-test at each time point
	ticks <- seq(0, max(time.all), round(max(time.all) / 10))
	plot(NULL, ,xlim = c(0, max(time.all)), ylim = c(0,1), ylab='Proportion of Fixations',
		xlab='Time', axes=FALSE, main='Logistic Curve')
	axis(1, at = ticks)
	axis(2)
	box()
	legend('topleft', lty = 1:2, legend = groups)
	
	#Entries in tsmultcomp:
	#1 : Estimate of rho
	#2 : Overall Type I Error
	#3 : Total tests we will perform
	rho.est <- ar(t.val, FALSE, order.max = 1)$ar
	if(p.adj == "oleson") {
		if(paired) {
			alphastar <- tsmultcomp(rho.est, alpha, N.tests, df = N.g1 - 1)
		} else {
			alphastar <- tsmultcomp(rho.est, alpha, N.tests, df = N.g1 + N.g2 - 2)
		}
		sig <- p.values <= alphastar
	} else if(p.adj == "fdr") {
		sig <- p.adjust(p.values, "fdr") <= alpha
	} else if(p.adj == "none") {
		sig <- p.values <= alpha
	}

	#Make significant area yellow
	buck <- bucket(sig, time.all, ylim = c(0, .9))
	
	#Plot overall estimate of curves
	lines(time.all, curve.g1, lty = 1, lwd = 2)
	lines(time.all, curve.g2, lty = 2, lwd = 2)
	
	#Plot Confidence Interval for Group 1 curve
	lines(time.all, curve.g1 - curve.sd1 * qt(alpha / 2, N.g1 - 1), lty = 1, lwd = 1,
		col = "gray44")
	lines(time.all, curve.g1 + curve.sd1 * qt(alpha / 2, N.g1 - 1), lty = 1, lwd = 1,
		col = "gray44")
		
	#Plot Confidence Interval for Group 2 curve
	lines(time.all, curve.g2 - curve.sd2 * qt(alpha / 2, N.g2 - 1), lty = 2, lwd = 1,
		col = "gray44")
	lines(time.all, curve.g2 + curve.sd2 * qt(alpha / 2, N.g2 - 1), lty = 2, lwd = 1,
		col = "gray44")
		
	# Record confidence intervals
	curve.ci1 <- curve.ci2 <- matrix(NA, nrow = length(time.all), ncol = 4)
	curve.ci1[,1] <- curve.ci2[,1] <- time.all
	curve.ci1[,2] <- curve.g1 - curve.sd1 * qt(1 - alpha / 2, N.g1 - 1)
	curve.ci1[,3] <- curve.g1
	curve.ci1[,4] <- curve.g1 + curve.sd1 * qt(1 - alpha / 2, N.g1 - 1)
	curve.ci2[,2] <- curve.g2 - curve.sd2 * qt(1 - alpha / 2, N.g2 - 1)
	curve.ci2[,3] <- curve.g2
	curve.ci2[,4] <- curve.g2 + curve.sd2 * qt(1 - alpha / 2, N.g2 - 1)
	colnames(curve.ci1) <- colnames(curve.ci2) <- c("Time", "Lower CI", "Estimate", "Upper CI")
		
	if(!is.null(time.test)) {
		time.test <- which(time.all %in% time.test)
		
		cat("######################\n")
		cat("## Individual Tests ##\n")
		cat("######################\n")
		
		for(i in 1:length(time.test)) {
			time <- time.test[i]
			mean.1 <- curve.g1[time]
			mean.2 <- curve.g2[time]
			sd.1 <- curve.sd1[time]
			sd.2 <- curve.sd2[time]
			
			time.mean <- mean.1 - mean.2
			time.se <- sqrt((N.g1 * (N.g1 - 1) * sd.1 ^ 2 + N.g2 * (N.g2 - 1) * sd.2 ^ 2) / (N.g1 + N.g2 - 2) * (1 / N.g1 + 1 / N.g2))
			time.df <- N.g1 + N.g2 - 2
			time.t <- time.mean / time.se
			time.p <- pt(abs(time.t), time.df, lower.tail = FALSE) * 2
			pooled.sd <- sqrt((N.g1 * (N.g1 - 1) * sd.1 ^ 2 + N.g2 * (N.g2 - 1) * sd.2 ^ 2) / (N.g1 + N.g2 - 2))
			time.d <- time.mean / pooled.sd
			
			cat(paste0("Test # = ", i, " --- Time = ", time.all[time], "\n"))
			cat(paste0("Mean Diff = ", round(time.mean, 4), " --- SE = ", round(time.se, 4), "\n"))
			if(time.p < .0001) {
				cat(paste0("t = ", round(time.t, 2), " --- DF = ", round(time.df, 1), " --- p < 0.0001 \n"))
			} else {
				cat(paste0("t = ", round(time.t, 2), " --- DF = ", round(time.df, 1), " --- p = ", round(time.p, 4), "\n"))
			}
			cat(paste0("Pooled SD = ", round(pooled.sd, 4), " --- Cohen's d = ", round(time.d, 1), "\n\n"))
		}
	}
	
	if(test.params) {
		if(paired) {
			cat("######################\n")
			cat("## Parameter Tests  ##\n")
			cat("##  Paired t-test   ##\n")
			cat("######################\n")
			
			df <- N.g1 - 1
			
			mini <- mini.1 - mini.2
			mini.mean <- mean(mini)
			mini.se <- sd(mini)
			mini.t <- mini.mean / mini.se
			mini.p <- pt(abs(mini.t), df, lower.tail = FALSE) * 2
			
			peak <- peak.1 - peak.2
			peak.mean <- mean(peak)
			peak.se <- sd(peak)
			peak.t <- peak.mean / peak.se
			peak.p <- pt(abs(peak.t), df, lower.tail = FALSE) * 2
			
			slope <- slope.1 - slope.2
			slope.mean <- mean(slope)
			slope.se <- sd(slope)
			slope.t <- slope.mean / slope.se
			slope.p <- pt(abs(slope.t), df, lower.tail = FALSE) * 2
			
			cross <- cross.1 - cross.2
			cross.mean <- mean(cross)
			cross.se <- sd(cross)
			cross.t <- cross.mean / cross.se
			cross.p <- pt(abs(cross.t), df, lower.tail = FALSE) * 2
			
			cat(paste0("Mini -- Diff: ", round(mini.mean, 4), ", t: ", round(mini.t, 3), ", SE: ", round(mini.se, 3), ", df: ", df, ", p: ", round(mini.p, 4), "\n"))
			cat(paste0("Peak -- Diff: ", round(peak.mean, 4), ", t: ", round(peak.t, 3), ", SE: ", round(peak.se, 3), ", df: ", df,  ", p: ", round(peak.p, 4), "\n"))
			cat(paste0("Slope -- Diff: ", round(slope.mean, 4), ", t: ", round(slope.t, 3), ", SE: ", round(slope.se, 3), ", df: ", df,  ", p: ", round(slope.p, 4), "\n"))
			cat(paste0("Cross -- Diff: ", round(cross.mean, 4), ", t: ", round(cross.t, 3), ", SE: ", round(cross.se, 3), ", df: ", df,  ", p: ", round(cross.p, 4), "\n\n"))
		} else {
			cat("######################\n")
			cat("## Parameter Tests  ##\n")
			cat("## 2 Sample t-test  ##\n")
			cat("######################\n")
			
			df <- N.g1 + N.g2 - 2
		
			mini.mean <- mean(mini.1) - mean(mini.2)
			mini.se1 <- sd(mini.1)
			mini.se2 <- sd(mini.2)
			mini.se <- sqrt((N.g1 * (N.g1 - 1) * mini.se1 ^ 2 + N.g2 * (N.g2 - 1) * mini.se2 ^ 2) / (N.g1 + N.g2 - 2) * (1 / N.g1 + 1 / N.g2))
			mini.t <- mini.mean / mini.se
			mini.p <- pt(abs(mini.t), df, lower.tail = FALSE) * 2
			
			peak.mean <- mean(peak.1) - mean(peak.2)
			peak.se1 <- sd(peak.1)
			peak.se2 <- sd(peak.2)
			peak.se <- sqrt((N.g1 * (N.g1 - 1) * peak.se1 ^ 2 + N.g2 * (N.g2 - 1) * peak.se2 ^ 2) / (N.g1 + N.g2 - 2) * (1 / N.g1 + 1 / N.g2))
			peak.t <- peak.mean / peak.se
			peak.p <- pt(abs(peak.t), df, lower.tail = FALSE) * 2
			
			slope.mean <- mean(slope.1) - mean(slope.2)
			slope.se1 <- sd(slope.1)
			slope.se2 <- sd(slope.2)
			slope.se <- sqrt((N.g1 * (N.g1 - 1) * slope.se1 ^ 2 + N.g2 * (N.g2 - 1) * slope.se2 ^ 2) / (N.g1 + N.g2 - 2) * (1 / N.g1 + 1 / N.g2))
			slope.t <- slope.mean / slope.se
			slope.p <- pt(abs(slope.t), df, lower.tail = FALSE) * 2
			
			cross.mean <- mean(cross.1) - mean(cross.2)
			cross.se1 <- sd(cross.1)
			cross.se2 <- sd(cross.2)
			cross.se <- sqrt((N.g1 * (N.g1 - 1) * cross.se1 ^ 2 + N.g2 * (N.g2 - 1) * cross.se2 ^ 2) / (N.g1 + N.g2 - 2) * (1 / N.g1 + 1 / N.g2))
			cross.t <- cross.mean / cross.se
			cross.p <- pt(abs(cross.t), df, lower.tail = FALSE) * 2
			
			cat(paste0("Mini -- Diff: ", round(mini.mean, 4), ", t: ", round(mini.t, 3), ", SE: ", round(mini.se, 3), ", df: ", df, ", p: ", round(mini.p, 4), "\n"))
			cat(paste0("Peak -- Diff: ", round(peak.mean, 4), ", t: ", round(peak.t, 3), ", SE: ", round(peak.se, 3), ", df: ", df,  ", p: ", round(peak.p, 4), "\n"))
			cat(paste0("Slope -- Diff: ", round(slope.mean, 4), ", t: ", round(slope.t, 3), ", SE: ", round(slope.se, 3), ", df: ", df,  ", p: ", round(slope.p, 4), "\n"))
			cat(paste0("Cross -- Diff: ", round(cross.mean, 4), ", t: ", round(cross.t, 3), ", SE: ", round(cross.se, 3), ", df: ", df,  ", p: ", round(cross.p, 4), "\n\n"))
		}
	}
	
	params <- round(c(alpha = alpha, alpha.adj = alphastar, rho.est = rho.est), 4)
	print(list(alpha = params, significant = buck))
	invisible(list(alpha = params, significant = buck, time.all = time.all,sig = sig,
		curve.ci1 = curve.ci1, curve.ci2 = curve.ci2, curve.g1 = curve.g1, curve.g2 = curve.g2,
		curve.sd1 = curve.sd1, curve.sd2 = curve.sd2, N.g1 = N.g1, N.g2 = N.g2,
    curve1.mat = curve1.mat, curve2.mat = curve2.mat, groups = groups, seed = seed))
}
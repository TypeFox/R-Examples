doubleGauss.boot <- function(part1.list, seed = new.seed(), alpha = 0.05, paired = FALSE, N.iter = 1000, cores = 1, p.adj = "oleson", test.spots = NULL, time.test = NULL, test.params = FALSE) {
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
	
	group1.bad <- is.na(coef.id1[,1]) | is.na(coef.id3[,1])
	group2.bad <- is.na(coef.id2[,1]) | is.na(coef.id4[,1])
	
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

	curve1.0 <-  matrix(NA, ncol = N.time, nrow = N.g1)
	mu1.ran <-   rep(NA, N.g1)
	ht1.ran <-   rep(NA, N.g1)
	s11.ran <-   rep(NA, N.g1)
	s21.ran <-   rep(NA, N.g1)
	b11.ran <-   rep(NA, N.g1)
	b21.ran <-   rep(NA, N.g1)

	curve2.0 <-  matrix(NA, ncol = N.time, nrow = N.g2)
	mu2.ran <-   rep(NA, N.g2)
	ht2.ran <-   rep(NA, N.g2)
	s12.ran <-   rep(NA, N.g2)
	s22.ran <-   rep(NA, N.g2)
	b12.ran <-   rep(NA, N.g2)
	b22.ran <-   rep(NA, N.g2)

	curve1.mat <- matrix(NA, ncol = N.time, nrow = N.iter)
	curve2.mat <- matrix(NA, ncol = N.time, nrow = N.iter)
	curve3.mat <- matrix(NA, ncol = N.time, nrow = N.iter)
  
  #Target Curve
	curve.f <- function(mu, ht, sig1, sig2, base1, base2, x){
	  whichgauss <- x < mu
	  y1 <- exp(-1 * (x - mu) ^ 2 / (2 * sig1 ^ 2)) * (ht - base1) + base1
	  y2 <- exp(-1 * (x - mu) ^ 2 / (2 * sig2 ^ 2)) * (ht - base2) + base2
	
	  y <- whichgauss * y1 + (1 - whichgauss) * y2
		y
	}
	
	mu.1 <- ht.1 <- s1.1 <- s2.1 <- b1.1 <- b2.1 <-
		mu.2 <- ht.2 <- s1.2 <- s2.2 <- b1.2 <- b2.2 <- numeric(N.iter)
	
	##################
	##### 1 Core #####
	##################
	
	if(paired) {
		mu.cov.1 <- cov(coef.id1[,1], coef.id2[,1], use = "pairwise.complete.obs")
		ht.cov.1 <- cov(coef.id1[,2], coef.id2[,2], use = "pairwise.complete.obs")
		s1.cov.1 <- cov(coef.id1[,3], coef.id2[,3], use = "pairwise.complete.obs")
		s2.cov.1 <- cov(coef.id1[,4], coef.id2[,4], use = "pairwise.complete.obs")
		b1.cov.1 <- cov(coef.id1[,5], coef.id2[,5], use = "pairwise.complete.obs")
		b2.cov.1 <- cov(coef.id1[,6], coef.id2[,6], use = "pairwise.complete.obs")
		
		if(diffs) {
			mu.cov.2 <- cov(coef.id3[,1], coef.id4[,1], use = "pairwise.complete.obs")
			ht.cov.2 <- cov(coef.id3[,2], coef.id4[,2], use = "pairwise.complete.obs")
			s1.cov.2 <- cov(coef.id3[,3], coef.id4[,3], use = "pairwise.complete.obs")
			s2.cov.2 <- cov(coef.id3[,4], coef.id4[,4], use = "pairwise.complete.obs")
			b1.cov.2 <- cov(coef.id3[,5], coef.id4[,5], use = "pairwise.complete.obs")
			b2.cov.2 <- cov(coef.id3[,6], coef.id4[,6], use = "pairwise.complete.obs")
		}
	}

	if(cores == 1) {
		set.seed(seed)
		for(iter in 1:N.iter){
			if(paired) {
				for(i in 1:N.g1) {
					mu <- rmvnorm(1, mean = c(coef.id1[i,1], coef.id2[i,1]),
						sigma = matrix(c(sdev.id1[i,1] ^ 2, mu.cov.1, mu.cov.1, sdev.id2[i,1] ^ 2), nrow = 2))
					ht <- rmvnorm(1, mean = c(coef.id1[i,2], coef.id2[i,2]),
						sigma = matrix(c(sdev.id1[i,2] ^ 2, ht.cov.1, ht.cov.1, sdev.id2[i,2] ^ 2), nrow = 2))
					s1 <- rmvnorm(1, mean = c(coef.id1[i,3], coef.id2[i,3]),
						sigma = matrix(c(sdev.id1[i,3] ^ 2, s1.cov.1, s1.cov.1, sdev.id2[i,3] ^ 2), nrow = 2))
					s2 <- rmvnorm(1, mean = c(coef.id1[i,4], coef.id2[i,4]),
						sigma = matrix(c(sdev.id1[i,4] ^ 2, s2.cov.1, s2.cov.1, sdev.id2[i,4] ^ 2), nrow = 2))
					b1 <- rmvnorm(1, mean = c(coef.id1[i,5], coef.id2[i,5]),
						sigma = matrix(c(sdev.id1[i,5] ^ 2, b1.cov.1, b1.cov.1, sdev.id2[i,5] ^ 2), nrow = 2))
					b2 <- rmvnorm(1, mean = c(coef.id1[i,6], coef.id2[i,6]),
						sigma = matrix(c(sdev.id1[i,6] ^ 2, b2.cov.1, b2.cov.1, sdev.id2[i,6] ^ 2), nrow = 2))
						
					mu1.ran <- mu[1]; mu2.ran <- mu[2]
					ht1.ran <- ht[1]; ht2.ran <- ht[2]
					s11.ran <- s1[1]; s12.ran <- s1[2]
					s21.ran <- s2[1]; s22.ran <- s2[2]
					b11.ran <- b1[1]; b12.ran <- b1[2]
					b21.ran <- b2[1]; b22.ran <- b2[2]
				}
			} else {
				mu1.ran <- rnorm(N.g1, coef.id1[,1], sdev.id1[,1])
				ht1.ran <- rnorm(N.g1, coef.id1[,2], sdev.id1[,2])
				s11.ran <- rnorm(N.g1, coef.id1[,3], sdev.id1[,3])
				s21.ran <- rnorm(N.g1, coef.id1[,4], sdev.id1[,4])
				b11.ran <- rnorm(N.g1, coef.id1[,5], sdev.id1[,5])
				b21.ran <- rnorm(N.g1, coef.id1[,6], sdev.id1[,6])

				mu2.ran <- rnorm(N.g2, coef.id2[,1], sdev.id2[,1])
				ht2.ran <- rnorm(N.g2, coef.id2[,2], sdev.id2[,2])
				s12.ran <- rnorm(N.g2, coef.id2[,3], sdev.id2[,3])
				s22.ran <- rnorm(N.g2, coef.id2[,4], sdev.id2[,4])
				b12.ran <- rnorm(N.g2, coef.id2[,5], sdev.id2[,5])
				b22.ran <- rnorm(N.g2, coef.id2[,6], sdev.id2[,6])
			}
			
			mu.1[iter] <- mean(mu1.ran); mu.2[iter] <- mean(mu2.ran)
			ht.1[iter] <- mean(ht1.ran); ht.2[iter] <- mean(ht1.ran)
			s1.1[iter] <- mean(s11.ran); s1.2[iter] <- mean(s12.ran)
			s2.1[iter] <- mean(s21.ran); s2.2[iter] <- mean(s22.ran)
			b1.1[iter] <- mean(b11.ran); b1.2[iter] <- mean(b12.ran)
			b2.1[iter] <- mean(b21.ran); b2.2[iter] <- mean(b22.ran)
			
			if(diffs) {
				if(paired) {
					for(i in 1:N.g1) {
						mu <- rmvnorm(1, mean = c(coef.id3[i,1], coef.id4[i,1]),
							sigma = matrix(c(sdev.id1[i,1] ^ 2, mu.cov.2, mu.cov.2, sdev.id2[i,1] ^ 2), nrow = 2))
						ht <- rmvnorm(1, mean = c(coef.id3[i,2], coef.id4[i,2]),
							sigma = matrix(c(sdev.id1[i,2] ^ 2, ht.cov.2, ht.cov.2, sdev.id2[i,2] ^ 2), nrow = 2))
						s1 <- rmvnorm(1, mean = c(coef.id3[i,3], coef.id4[i,3]),
							sigma = matrix(c(sdev.id1[i,3] ^ 2, s1.cov.2, s1.cov.2, sdev.id2[i,3] ^ 2), nrow = 2))
						s2 <- rmvnorm(1, mean = c(coef.id3[i,4], coef.id4[i,4]),
							sigma = matrix(c(sdev.id1[i,4] ^ 2, s2.cov.2, s2.cov.2, sdev.id2[i,4] ^ 2), nrow = 2))
						b1 <- rmvnorm(1, mean = c(coef.id3[i,5], coef.id4[i,5]),
							sigma = matrix(c(sdev.id1[i,5] ^ 2, b1.cov.2, b1.cov.2, sdev.id2[i,5] ^ 2), nrow = 2))
						b2 <- rmvnorm(1, mean = c(coef.id3[i,6], coef.id4[i,6]),
							sigma = matrix(c(sdev.id1[i,6] ^ 2, b2.cov.2, b2.cov.2, sdev.id2[i,6] ^ 2), nrow = 2))
							
						mu3.ran <- mu[1]; mu4.ran <- mu[2]
						ht3.ran <- ht[1]; ht4.ran <- ht[2]
						s13.ran <- s1[1]; s14.ran <- s1[2]
						s23.ran <- s2[1]; s24.ran <- s2[2]
						b13.ran <- b1[1]; b14.ran <- b1[2]
						b23.ran <- b2[1]; b24.ran <- b2[2]
					}
				} else {
					mu3.ran <- rnorm(N.g1, coef.id3[,1], sdev.id3[,1])
					ht3.ran <- rnorm(N.g1, coef.id3[,2], sdev.id3[,2])
					s13.ran <- rnorm(N.g1, coef.id3[,3], sdev.id3[,3])
					s23.ran <- rnorm(N.g1, coef.id3[,4], sdev.id3[,4])
					b13.ran <- rnorm(N.g1, coef.id3[,5], sdev.id3[,5])
					b23.ran <- rnorm(N.g1, coef.id3[,6], sdev.id3[,6])

					mu4.ran <- rnorm(N.g2, coef.id4[,1], sdev.id4[,1])
					ht4.ran <- rnorm(N.g2, coef.id4[,2], sdev.id4[,2])
					s14.ran <- rnorm(N.g2, coef.id4[,3], sdev.id4[,3])
					s24.ran <- rnorm(N.g2, coef.id4[,4], sdev.id4[,4])
					b14.ran <- rnorm(N.g2, coef.id4[,5], sdev.id4[,5])
					b24.ran <- rnorm(N.g2, coef.id4[,6], sdev.id4[,6])
				}
			}

			if(diffs) {
				for(id in 1:N.g1) { #Get fixation level for each group 1 subject
					curve1.0[id,] <- curve.f(mu1.ran[id], ht1.ran[id], s11.ran[id], s21.ran[id],
																		b11.ran[id], b21.ran[id], time.all) -
													 curve.f(mu3.ran[id], ht3.ran[id], s13.ran[id], s23.ran[id],
																		b13.ran[id], b23.ran[id], time.all)
				}
				for(id in 1:N.g2) { #Get fixation level for each group 2 subject
					curve2.0[id,] <- curve.f(mu2.ran[id], ht2.ran[id], s12.ran[id], s22.ran[id],
																		b12.ran[id], b22.ran[id], time.all) -
													 curve.f(mu4.ran[id], ht4.ran[id], s14.ran[id], s24.ran[id],
																		b14.ran[id], b24.ran[id], time.all)
				}
			} else {
				for(id in 1:N.g1) { #Get fixation level for each group 1 subject
					curve1.0[id,] <- curve.f(mu1.ran[id], ht1.ran[id], s11.ran[id], s21.ran[id],
																		b11.ran[id], b21.ran[id], time.all)
				}
				for(id in 1:N.g2) { #Get fixation level for each group 2 subject
					curve2.0[id,] <- curve.f(mu2.ran[id], ht2.ran[id], s12.ran[id], s22.ran[id],
																		b12.ran[id], b22.ran[id], time.all)
				}
			}

			curve1.mat[iter,] <- apply(curve1.0, 2, mean) #Mean fixations at each time point for group 1
			curve2.mat[iter,] <- apply(curve2.0, 2, mean) #Mean fixations at each time point for group 2
			
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
					mu <- rmvnorm(1, mean = c(coef.id1[i,1], coef.id2[i,1]),
						sigma = matrix(c(sdev.id1[i,1] ^ 2, mu.cov.1, mu.cov.1, sdev.id2[i,1] ^ 2), nrow = 2))
					ht <- rmvnorm(1, mean = c(coef.id1[i,2], coef.id2[i,2]),
						sigma = matrix(c(sdev.id1[i,2] ^ 2, ht.cov.1, ht.cov.1, sdev.id2[i,2] ^ 2), nrow = 2))
					s1 <- rmvnorm(1, mean = c(coef.id1[i,3], coef.id2[i,3]),
						sigma = matrix(c(sdev.id1[i,3] ^ 2, s1.cov.1, s1.cov.1, sdev.id2[i,3] ^ 2), nrow = 2))
					s2 <- rmvnorm(1, mean = c(coef.id1[i,4], coef.id2[i,4]),
						sigma = matrix(c(sdev.id1[i,4] ^ 2, s2.cov.1, s2.cov.1, sdev.id2[i,4] ^ 2), nrow = 2))
					b1 <- rmvnorm(1, mean = c(coef.id1[i,5], coef.id2[i,5]),
						sigma = matrix(c(sdev.id1[i,5] ^ 2, b1.cov.1, b1.cov.1, sdev.id2[i,5] ^ 2), nrow = 2))
					b2 <- rmvnorm(1, mean = c(coef.id1[i,6], coef.id2[i,6]),
						sigma = matrix(c(sdev.id1[i,6] ^ 2, b2.cov.1, b2.cov.1, sdev.id2[i,6] ^ 2), nrow = 2))
						
					mu1.ran <- mu[1]; mu2.ran <- mu[2]
					ht1.ran <- ht[1]; ht2.ran <- ht[2]
					s11.ran <- s1[1]; s12.ran <- s1[2]
					s21.ran <- s2[1]; s22.ran <- s2[2]
					b11.ran <- b1[1]; b12.ran <- b1[2]
					b21.ran <- b2[1]; b22.ran <- b2[2]
				}
			} else {
				mu1.ran <- rnorm(N.g1, coef.id1[,1], sdev.id1[,1])
				ht1.ran <- rnorm(N.g1, coef.id1[,2], sdev.id1[,2])
				s11.ran <- rnorm(N.g1, coef.id1[,3], sdev.id1[,3])
				s21.ran <- rnorm(N.g1, coef.id1[,4], sdev.id1[,4])
				b11.ran <- rnorm(N.g1, coef.id1[,5], sdev.id1[,5])
				b21.ran <- rnorm(N.g1, coef.id1[,6], sdev.id1[,6])

				mu2.ran <- rnorm(N.g2, coef.id2[,1], sdev.id2[,1])
				ht2.ran <- rnorm(N.g2, coef.id2[,2], sdev.id2[,2])
				s12.ran <- rnorm(N.g2, coef.id2[,3], sdev.id2[,3])
				s22.ran <- rnorm(N.g2, coef.id2[,4], sdev.id2[,4])
				b12.ran <- rnorm(N.g2, coef.id2[,5], sdev.id2[,5])
				b22.ran <- rnorm(N.g2, coef.id2[,6], sdev.id2[,6])
			}
			
			mu.temp.1 <- mean(mu1.ran); mu.temp.2 <- mean(mu2.ran)
			ht.temp.1 <- mean(ht1.ran); ht.temp.2 <- mean(ht1.ran)
			s1.temp.1 <- mean(s11.ran); s1.temp.2 <- mean(s12.ran)
			s2.temp.1 <- mean(s21.ran); s2.temp.2 <- mean(s22.ran)
			b1.temp.1 <- mean(b11.ran); b1.temp.2 <- mean(b12.ran)
			b2.temp.1 <- mean(b21.ran); b2.temp.2 <- mean(b22.ran)
			
			if(diffs) {
				if(paired) {
					for(i in 1:N.g1) {
						mu <- rmvnorm(1, mean = c(coef.id3[i,1], coef.id4[i,1]),
							sigma = matrix(c(sdev.id1[i,1] ^ 2, mu.cov.2, mu.cov.2, sdev.id2[i,1] ^ 2), nrow = 2))
						ht <- rmvnorm(1, mean = c(coef.id3[i,2], coef.id4[i,2]),
							sigma = matrix(c(sdev.id1[i,2] ^ 2, ht.cov.2, ht.cov.2, sdev.id2[i,2] ^ 2), nrow = 2))
						s1 <- rmvnorm(1, mean = c(coef.id3[i,3], coef.id4[i,3]),
							sigma = matrix(c(sdev.id1[i,3] ^ 2, s1.cov.2, s1.cov.2, sdev.id2[i,3] ^ 2), nrow = 2))
						s2 <- rmvnorm(1, mean = c(coef.id3[i,4], coef.id4[i,4]),
							sigma = matrix(c(sdev.id1[i,4] ^ 2, s2.cov.2, s2.cov.2, sdev.id2[i,4] ^ 2), nrow = 2))
						b1 <- rmvnorm(1, mean = c(coef.id3[i,5], coef.id4[i,5]),
							sigma = matrix(c(sdev.id1[i,5] ^ 2, b1.cov.2, b1.cov.2, sdev.id2[i,5] ^ 2), nrow = 2))
						b2 <- rmvnorm(1, mean = c(coef.id3[i,6], coef.id4[i,6]),
							sigma = matrix(c(sdev.id1[i,6] ^ 2, b2.cov.2, b2.cov.2, sdev.id2[i,6] ^ 2), nrow = 2))
							
						mu3.ran <- mu[1]; mu4.ran <- mu[2]
						ht3.ran <- ht[1]; ht4.ran <- ht[2]
						s13.ran <- s1[1]; s14.ran <- s1[2]
						s23.ran <- s2[1]; s24.ran <- s2[2]
						b13.ran <- b1[1]; b14.ran <- b1[2]
						b23.ran <- b2[1]; b24.ran <- b2[2]
					}
				} else {
					mu3.ran <- rnorm(N.g1, coef.id3[,1], sdev.id3[,1])
					ht3.ran <- rnorm(N.g1, coef.id3[,2], sdev.id3[,2])
					s13.ran <- rnorm(N.g1, coef.id3[,3], sdev.id3[,3])
					s23.ran <- rnorm(N.g1, coef.id3[,4], sdev.id3[,4])
					b13.ran <- rnorm(N.g1, coef.id3[,5], sdev.id3[,5])
					b23.ran <- rnorm(N.g1, coef.id3[,6], sdev.id3[,6])

					mu4.ran <- rnorm(N.g2, coef.id4[,1], sdev.id4[,1])
					ht4.ran <- rnorm(N.g2, coef.id4[,2], sdev.id4[,2])
					s14.ran <- rnorm(N.g2, coef.id4[,3], sdev.id4[,3])
					s24.ran <- rnorm(N.g2, coef.id4[,4], sdev.id4[,4])
					b14.ran <- rnorm(N.g2, coef.id4[,5], sdev.id4[,5])
					b24.ran <- rnorm(N.g2, coef.id4[,6], sdev.id4[,6])
				}
			}
			
			if(diffs) {
				for(id in 1:N.g1) { #Get fixation level for each group 1 subject
					curve1.0[id,] <- curve.f(mu1.ran[id], ht1.ran[id], s11.ran[id], s21.ran[id],
																		b11.ran[id], b21.ran[id], time.all) -
													 curve.f(mu3.ran[id], ht3.ran[id], s13.ran[id], s23.ran[id],
																		b13.ran[id], b23.ran[id], time.all)
				}
				for(id in 1:N.g2) { #Get fixation level for each group 2 subject
					curve2.0[id,] <- curve.f(mu2.ran[id], ht2.ran[id], s12.ran[id], s22.ran[id],
																		b12.ran[id], b22.ran[id], time.all) -
													 curve.f(mu4.ran[id], ht4.ran[id], s14.ran[id], s24.ran[id],
																		b14.ran[id], b24.ran[id], time.all)
				}
			} else {
				for(id in 1:N.g1) { #Get fixation level for each group 1 subject
					curve1.0[id,] <- curve.f(mu1.ran[id], ht1.ran[id], s11.ran[id], s21.ran[id],
																		b11.ran[id], b21.ran[id], time.all)
				}
				for(id in 1:N.g2) { #Get fixation level for each group 2 subject
					curve2.0[id,] <- curve.f(mu2.ran[id], ht2.ran[id], s12.ran[id], s22.ran[id],
																		b12.ran[id], b22.ran[id], time.all)
				}
			}
			
			curve1 <- apply(curve1.0, 2, mean) #Mean fixations at each time point for CIs
			curve2 <- apply(curve2.0, 2, mean) #Mean fixations at each time point for NHs
			curve3 <- curve2 - curve1
			
			c(curve1, curve2, curve3, mu.temp.1, ht.temp.1, s1.temp.1, s2.temp.1, b1.temp.1, b2.temp.1,
				mu.temp.2, ht.temp.2, s1.temp.2, s2.temp.2, b1.temp.2, b2.temp.2)
		}
		
		curve1.mat <- for.out[,1:N.time]
		curve2.mat <- for.out[,(N.time + 1):(2 * N.time)]
		curve3.mat <- for.out[,(2 * N.time + 1):(3 * N.time)]
		mu.1 <- for.out[, 3 * N.time + 1]
		ht.1 <- for.out[, 3 * N.time + 2]
		s1.1 <- for.out[, 3 * N.time + 3]
		s2.1 <- for.out[, 3 * N.time + 4]
		b1.1 <- for.out[, 3 * N.time + 5]
		b2.1 <- for.out[, 3 * N.time + 6]
		mu.2 <- for.out[, 3 * N.time + 7]
		ht.2 <- for.out[, 3 * N.time + 8]
		s1.2 <- for.out[, 3 * N.time + 9]
		s2.2 <- for.out[, 3 * N.time + 10]
		b1.2 <- for.out[, 3 * N.time + 11]
		b2.2 <- for.out[, 3 * N.time + 12]
		
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

	par(mfrow = c(1,1))
	# compute t-test at each time point
	ticks <- seq(0, max(time.all), round(max(time.all) / 10))
	plot(NULL, ,xlim = c(0, max(time.all)), ylim = c(0,1), ylab = 'Proportion of Fixations',
		xlab = 'Time', axes = FALSE, main = 'Double-Gauss Curve')
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
			
			mu <- mu.1 - mu.2
			mu.mean <- mean(mu)
			mu.se <- sd(mu)
			mu.t <- mu.mean / mu.se
			mu.p <- pt(abs(mu.t), df, lower.tail = FALSE) * 2
			
			ht <- ht.1 - ht.2
			ht.mean <- mean(ht)
			ht.se <- sd(ht)
			ht.t <- ht.mean / ht.se
			ht.p <- pt(abs(ht.t), df, lower.tail = FALSE) * 2
			
			s1 <- s1.1 - s1.2
			s1.mean <- mean(s1)
			s1.se <- sd(s1)
			s1.t <- s1.mean / s1.se
			s1.p <- pt(abs(s1.t), df, lower.tail = FALSE) * 2
			
			s2 <- s2.1 - s2.2
			s2.mean <- mean(s2)
			s2.se <- sd(s2)
			s2.t <- s2.mean / s2.se
			s2.p <- pt(abs(s2.t), df, lower.tail = FALSE) * 2
			
			b1 <- b1.1 - b1.2
			b1.mean <- mean(b1)
			b1.se <- sd(b1)
			b1.t <- b1.mean / b1.se
			b1.p <- pt(abs(b1.t), df, lower.tail = FALSE) * 2
			
			b2 <- b2.1 - b2.2
			b2.mean <- mean(b2)
			b2.se <- sd(b2)
			b2.t <- b2.mean / b2.se
			b2.p <- pt(abs(b2.t), df, lower.tail = FALSE) * 2
			
			cat(paste0("Mu -- Diff: ", round(mu.mean, 4), ", t: ", round(mu.t, 3), ", SE: ", round(mu.se, 3), ", df: ", df, ", p: ", round(mu.p, 4), "\n"))
			cat(paste0("Height -- Diff: ", round(ht.mean, 4), ", t: ", round(ht.t, 3), ", SE: ", round(ht.se, 3), ", df: ", df,  ", p: ", round(ht.p, 4), "\n"))
			cat(paste0("SD 1 -- Diff: ", round(s1.mean, 4), ", t: ", round(s1.t, 3), ", SE: ", round(s1.se, 3), ", df: ", df,  ", p: ", round(s1.p, 4), "\n"))
			cat(paste0("SD 2 -- Diff: ", round(s2.mean, 4), ", t: ", round(s2.t, 3), ", SE: ", round(s2.se, 3), ", df: ", df,  ", p: ", round(s2.p, 4), "\n"))
			cat(paste0("Base 1 -- Diff: ", round(b1.mean, 4), ", t: ", round(b1.t, 3), ", SE: ", round(b1.se, 3), ", df: ", df,  ", p: ", round(b1.p, 4), "\n"))
			cat(paste0("Base 2 -- Diff: ", round(b2.mean, 4), ", t: ", round(b2.t, 3), ", SE: ", round(b2.se, 3), ", df: ", df,  ", p: ", round(b2.p, 4), "\n\n"))
		} else {
			cat("######################\n")
			cat("## Parameter Tests  ##\n")
			cat("## 2 Sample t-test  ##\n")
			cat("######################\n")
			
			df <- N.g1 + N.g2 - 2
		
			mu.mean <- mean(mu.1) - mean(mu.2)
			mu.se1 <- sd(mu.1)
			mu.se2 <- sd(mu.2)
			mu.se <- sqrt((N.g1 * (N.g1 - 1) * mu.se1 ^ 2 + N.g2 * (N.g2 - 1) * mu.se2 ^ 2) / (N.g1 + N.g2 - 2) * (1 / N.g1 + 1 / N.g2))
			mu.t <- mu.mean / mu.se
			mu.p <- pt(abs(mu.t), df, lower.tail = FALSE) * 2
			
			ht.mean <- mean(ht.1) - mean(ht.2)
			ht.se1 <- sd(ht.1)
			ht.se2 <- sd(ht.2)
			ht.se <- sqrt((N.g1 * (N.g1 - 1) * ht.se1 ^ 2 + N.g2 * (N.g2 - 1) * ht.se2 ^ 2) / (N.g1 + N.g2 - 2) * (1 / N.g1 + 1 / N.g2))
			ht.t <- ht.mean / ht.se
			ht.p <- pt(abs(ht.t), df, lower.tail = FALSE) * 2
			
			s1.mean <- mean(s1.1) - mean(s1.2)
			s1.se1 <- sd(s1.1)
			s1.se2 <- sd(s1.2)
			s1.se <- sqrt((N.g1 * (N.g1 - 1) * s1.se1 ^ 2 + N.g2 * (N.g2 - 1) * s1.se2 ^ 2) / (N.g1 + N.g2 - 2) * (1 / N.g1 + 1 / N.g2))
			s1.t <- s1.mean / s1.se
			s1.p <- pt(abs(s1.t), df, lower.tail = FALSE) * 2
			
			s2.mean <- mean(s2.1) - mean(s2.2)
			s2.se1 <- sd(s2.1)
			s2.se2 <- sd(s2.2)
			s2.se <- sqrt((N.g1 * (N.g1 - 1) * s2.se1 ^ 2 + N.g2 * (N.g2 - 1) * s2.se2 ^ 2) / (N.g1 + N.g2 - 2) * (1 / N.g1 + 1 / N.g2))
			s2.t <- s2.mean / s2.se
			s2.p <- pt(abs(s2.t), df, lower.tail = FALSE) * 2
			
			b1.mean <- mean(b1.1) - mean(b1.2)
			b1.se1 <- sd(b1.1)
			b1.se2 <- sd(b1.2)
			b1.se <- sqrt((N.g1 * (N.g1 - 1) * b1.se1 ^ 2 + N.g2 * (N.g2 - 1) * b1.se2 ^ 2) / (N.g1 + N.g2 - 2) * (1 / N.g1 + 1 / N.g2))
			b1.t <- b1.mean / b1.se
			b1.p <- pt(abs(b1.t), df, lower.tail = FALSE) * 2
			
			b2.mean <- mean(b2.1) - mean(b2.2)
			b2.se1 <- sd(b2.1)
			b2.se2 <- sd(b2.2)
			b2.se <- sqrt((N.g1 * (N.g1 - 1) * b2.se1 ^ 2 + N.g2 * (N.g2 - 1) * b2.se2 ^ 2) / (N.g1 + N.g2 - 2) * (1 / N.g1 + 1 / N.g2))
			b2.t <- b2.mean / b2.se
			b2.p <- pt(abs(b2.t), df, lower.tail = FALSE) * 2
			
			cat(paste0("Mu -- Diff: ", round(mu.mean, 4), ", t: ", round(mu.t, 3), ", SE: ", round(mu.se, 3), ", df: ", df, ", p: ", round(mu.p, 4), "\n"))
			cat(paste0("Height -- Diff: ", round(ht.mean, 4), ", t: ", round(ht.t, 3), ", SE: ", round(ht.se, 3), ", df: ", df,  ", p: ", round(ht.p, 4), "\n"))
			cat(paste0("SD 1 -- Diff: ", round(s1.mean, 4), ", t: ", round(s1.t, 3), ", SE: ", round(s1.se, 3), ", df: ", df,  ", p: ", round(s1.p, 4), "\n"))
			cat(paste0("SD 2 -- Diff: ", round(s2.mean, 4), ", t: ", round(s2.t, 3), ", SE: ", round(s2.se, 3), ", df: ", df,  ", p: ", round(s2.p, 4), "\n"))
			cat(paste0("Base 1 -- Diff: ", round(b1.mean, 4), ", t: ", round(b1.t, 3), ", SE: ", round(b1.se, 3), ", df: ", df,  ", p: ", round(b1.p, 4), "\n"))
			cat(paste0("Base 2 -- Diff: ", round(b2.mean, 4), ", t: ", round(b2.t, 3), ", SE: ", round(b2.se, 3), ", df: ", df,  ", p: ", round(b2.p, 4), "\n\n"))
		}
	}
	
	params <- round(c(alpha = alpha, alpha.adj = alphastar, rho.est = rho.est), 4)
	print(list(alpha = params, significant = buck))
	invisible(list(alpha = params, significant = buck, time.all = time.all,sig = sig,
		curve.ci1 = curve.ci1, curve.ci2 = curve.ci2, curve.g1 = curve.g1, curve.g2 = curve.g2,
		curve.sd1 = curve.sd1, curve.sd2 = curve.sd2, N.g1 = N.g1, N.g2 = N.g2,
    curve1.mat = curve1.mat, curve2.mat = curve2.mat, groups = groups, seed = seed))
}
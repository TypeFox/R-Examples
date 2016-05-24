doubleGauss.fit <- function(data, col, concave = TRUE, diffs = FALSE, rho.0 = 0.9, cor = TRUE, cores = 1) {
	#Check that all required columns are included
	if(!("Subject"  %in% names(data)))       stop("data needs to include a 'Subject' column")
	if(!("Time"     %in% names(data)))       stop("data needs to include a 'Time' column")
	if(!("Group"    %in% names(data)))       stop("data needs to include a 'Group' column")
	if(diffs && !("Curve" %in% names(data))) stop("data needs to include a 'Curve' column if diffs=TRUE")
	if(ncol(data) < col)                     stop("specified col is too large for given data")
	if(!is.numeric(data[,col]))              stop("specified column is not numeric")
	if(cor && (rho.0 < 0 || rho.0 > 1)) 		 stop("rho.0 should be in the interval [0,1]")
	
	if(length(concave) == 1) concave <- rep(concave, 2)
	
	#Set overall time
	time.all <- sort(unique(data$Time))
	
	#Set groups
	groups <- unique(data$Group)
	if(length(groups) != 2) stop(paste("Expecting 2 unique groups. Actual number:", length(groups)))
	
	#Subjects per group
	#Group 1
	id.nums.g1 <- unique(data$Subject[data$Group == groups[1]])
	#Group 2
	id.nums.g2 <- unique(data$Subject[data$Group == groups[2]])
	
	# Check obs per subject
	if(!diffs) {
		counts <- aggregate(data[,col], list(data$Subject, data$Time, data$Group), length)
		if(any(counts$x > 1)) {
			warning("Some subjects have multiple observations at time slots, will average these")
			agg <- aggregate(data[,col], list(data$Subject, data$Time, data$Group), mean)
			names(agg) <- c("Subject", "Time", "Group", names(data)[col])
			data <- agg
		}
	} else {
		counts <- aggregate(data[,col], list(data$Subject, data$Time, data$Group, data$Curve), length)
		if(any(counts$x > 1)) {
			warning("Some subjects have multiple observations per curve at time slots, will average these")
			agg <- aggregate(data[,col], list(data$Subject, data$Time, data$Group, data$Curve), mean)
			names(agg) <- c("Subject", "Time", "Group", "Curve", names(data)[col])
			data <- agg
		}
	}

	# Size of groups
	N.g1 <- length(id.nums.g1)
	N.g2 <- length(id.nums.g2)
	
	#   pick starting values
	N.time <- 1 + max(data$Time) / 4 # 751 for 0-3000 ms
	N.sub1 <- N.g1
	N.sub2 <- N.g2

	# Save parameter estimates
	coef.id1 <-  matrix(NA, ncol = 6, nrow = N.sub1)
	sdev.id1 <-  matrix(NA, ncol = 6, nrow = N.sub1)
	sigma.id1 <- matrix(NA, ncol = 1, nrow = N.sub1)

	coef.id2 <-  matrix(NA, ncol = 6, nrow = N.sub2)
	sdev.id2 <-  matrix(NA, ncol = 6, nrow = N.sub2)
	sigma.id2 <- matrix(NA, ncol = 1, nrow = N.sub2)
	
	coef.id3 <-  matrix(NA, ncol = 6, nrow = N.sub1)
	sdev.id3 <-  matrix(NA, ncol = 6, nrow = N.sub1)
	sigma.id3 <- matrix(NA, ncol = 1, nrow = N.sub1)
	
	coef.id4 <-  matrix(NA, ncol = 6, nrow = N.sub2)
	sdev.id4 <-  matrix(NA, ncol = 6, nrow = N.sub2)
	sigma.id4 <- matrix(NA, ncol = 1, nrow = N.sub2)
	
	dgauss <- function(time, mu, ht, sig1, sig2, base1, base2) {
		(time < mu) * (exp(-1 * (time - mu) ^ 2 / (2 * sig1 ^ 2)) * (ht - base1) +
		base1) + (mu <= time) * (exp(-1 * (time - mu) ^ 2 / (2 * sig2 ^ 2)) * (ht - base2) +
		base2)
	}
	
	R2.g1.1 <- R2.g1.2 <- numeric(N.g1)
	R2.g2.1 <- R2.g2.2 <- numeric(N.g2)
	
	cor.1 <- cor.3 <- rep(cor, N.sub1)
	cor.2 <- cor.4 <- rep(cor, N.sub2)

	##### 1 Core #####
	if(cores == 1) {
		for(id in 1:N.g1){
			if(diffs) {
				y1id <- subset(data, data$Subject == id.nums.g1[id] & data$Group == groups[1] & data$Curve == 1)
			} else {
				y1id <- subset(data, data$Subject == id.nums.g1[id] & data$Group == groups[1])
			}
			
			y.fix <- y1id[,col]
			
			fit.curve <- est.dgauss.curve(time.all, y.fix, rho.0, concave[1], cor = cor)
			
			if(is.null(fit.curve$fit)) {
				coef.id1[id,] <- rep(NA, 6)
				sdev.id1[id,] <- rep(NA, 6)
				sigma.id1[id,] <- NA
				cor.1[id] <- NA
				R2.g1.1[id] <- NA
			} else {
				cor.1[id] <- fit.curve$cor
				fit.curve <- fit.curve$fit
			
				coef.id1[id,] <- coef(fit.curve)
				sdev.id1[id,] <- sqrt(diag(fit.curve$varBeta))
				sigma.id1[id,] <- fit.curve$sigma
						
				SSY <- sum((y.fix - mean(y.fix)) ^ 2)
				y.fit <- dgauss(time.all, coef(fit.curve)[1], coef(fit.curve)[2], coef(fit.curve)[3],
					coef(fit.curve)[4], coef(fit.curve)[5], coef(fit.curve)[6])
				y.err <- y.fit - y.fix
				SSE <- sum(y.err ^ 2)
				R2.g1.1[id] <- 1 - SSE / SSY
				
				print(paste0("Group = ", groups[1], ", ID = ", id, ", Subject = ", id.nums.g1[id],
											", Curve = 1, R2 = ", round(R2.g1.1[id], 3), ", AR1 = ", as.logical(cor.1[id])))
			}
			
			if(diffs) {
				y1id <- subset(data, data$Subject == id.nums.g1[id] & data$Group == groups[1] & data$Curve == 2)
				
				y.fix <- y1id[,col]
				
				fit.curve <- est.dgauss.curve(time.all, y.fix, rho.0, concave[1], cor = cor)
					
				if(is.null(fit.curve$fit)) {
					coef.id3[id,] <- rep(NA, 6)
					sdev.id3[id,] <- rep(NA, 6)
					sigma.id3[id,] <- NA
					cor.3[id] <- NA
					R2.g1.2[id] <- NA
				} else {
					cor.3[id] <- fit.curve$cor
					fit.curve <- fit.curve$fit
				
					coef.id3[id,] <- coef(fit.curve)
					sdev.id3[id,] <- sqrt(diag(fit.curve$varBeta))
					sigma.id3[id,] <- fit.curve$sigma
							
					SSY <- sum((y.fix - mean(y.fix)) ^ 2)
					y.fit <- dgauss(time.all, coef(fit.curve)[1], coef(fit.curve)[2], coef(fit.curve)[3],
						coef(fit.curve)[4], coef(fit.curve)[5], coef(fit.curve)[6])
					y.err <- y.fit - y.fix
					SSE <- sum(y.err ^ 2)
					R2.g1.2[id] <- 1 - SSE / SSY
					
					print(paste0("Group = ", groups[1], ", ID = ", id, ", Subject = ", id.nums.g1[id],
												", Curve = 2, R2 = ", round(R2.g1.2[id], 3), ", AR1 = ", as.logical(cor.3[id])))
				}
			}
		}

		for(id in 1:N.g2){
			if(diffs) {
				y2id <- subset(data, data$Subject == id.nums.g2[id] & data$Group == groups[2] & data$Curve == 1)
			} else {
				y2id <- subset(data, data$Subject == id.nums.g2[id] & data$Group == groups[2])
			}
			
			y.fix <- y2id[,col]
			
			fit.curve <- est.dgauss.curve(time.all, y.fix, rho.0, concave[2], cor = cor)
				
			if(is.null(fit.curve$fit)) {
				coef.id2[id,] <- rep(NA, 6)
				sdev.id2[id,] <- rep(NA, 6)
				sigma.id2[id,] <- NA
				cor.2[id] <- NA
				R2.g2.1[id] <- NA
			} else {
				cor.2[id] <- fit.curve$cor
				fit.curve <- fit.curve$fit
			
				coef.id2[id,] <- coef(fit.curve)
				sdev.id2[id,] <- sqrt(diag(fit.curve$varBeta))
				sigma.id2[id,] <- fit.curve$sigma
						
				SSY <- sum((y.fix - mean(y.fix)) ^ 2)
				y.fit <- dgauss(time.all, coef(fit.curve)[1], coef(fit.curve)[2], coef(fit.curve)[3],
					coef(fit.curve)[4], coef(fit.curve)[5], coef(fit.curve)[6])
				y.err <- y.fit - y.fix
				SSE <- sum(y.err ^ 2)
				R2.g2.1[id] <- 1 - SSE / SSY
				
				print(paste0("Group = ", groups[2], ", ID = ", id, ", Subject = ", id.nums.g2[id],
											", Curve = 1, R2 = ", round(R2.g2.1[id], 3), ", AR1 = ", as.logical(cor.2[id])))
			}

			if(diffs) {
				y2id <- subset(data, data$Subject == id.nums.g2[id] & data$Group == groups[2] & data$Curve == 2)
				
				y.fix <- y2id[,col]
			
				fit.curve <- est.dgauss.curve(time.all, y.fix, rho.0, concave[2], cor = cor)
					
				if(is.null(fit.curve$fit)) {
					coef.id4[id,] <- rep(NA, 6)
					sdev.id4[id,] <- rep(NA, 6)
					sigma.id4[id,] <- NA
					cor.4[id] <- NA
					R2.g2.2[id] <- NA
				} else {
					cor.4[id] <- fit.curve$cor
					fit.curve <- fit.curve$fit
				
					coef.id4[id,] <- coef(fit.curve)
					sdev.id4[id,] <- sqrt(diag(fit.curve$varBeta))
					sigma.id4[id,] <- fit.curve$sigma
							
					SSY <- sum((y.fix - mean(y.fix)) ^ 2)
					y.fit <- dgauss(time.all, coef(fit.curve)[1], coef(fit.curve)[2], coef(fit.curve)[3],
						coef(fit.curve)[4], coef(fit.curve)[5], coef(fit.curve)[6])
					y.err <- y.fit - y.fix
					SSE <- sum(y.err ^ 2)
					R2.g2.2[id] <- 1 - SSE / SSY
					
					print(paste0("Group = ", groups[2], ", ID = ", id, ", Subject = ", id.nums.g2[id],
												", Curve = 1, R2 = ", round(R2.g2.2[id], 3), ", AR1 = ", as.logical(cor.4[id])))
				}
			}
		}
	##### 2+ Cores #####
	} else {
		cl <- makePSOCKcluster(cores)
		registerDoParallel(cl)
		
		for.out <- foreach(id = 1:N.g1, .combine = rbind) %dopar% {
			if(diffs) {
				y1id <- subset(data, data$Subject == id.nums.g1[id] & data$Group == groups[1] & data$Curve == 1)
			} else {
				y1id <- subset(data, data$Subject == id.nums.g1[id] & data$Group == groups[1])
			}
			
			y.fix <- y1id[,col]
			
			fit.curve <- est.dgauss.curve(time.all, y.fix, rho.0, concave[1], cor = cor)
			
			if(is.null(fit.curve$fit)) {
				coef.1 <- rep(NA, 6)
				sdev.1 <- rep(NA, 6)
				sigma.1 <- NA
				cor.temp.1 <- NA
				R2.1 <- NA
			} else {
				cor.temp.1 <- fit.curve$cor
				fit.curve <- fit.curve$fit
			
				coef.1 <- coef(fit.curve)
				sdev.1 <- sqrt(diag(fit.curve$varBeta))
				sigma.1 <- fit.curve$sigma
						
				SSY <- sum((y.fix - mean(y.fix)) ^ 2)
				y.fit <- dgauss(time.all, coef(fit.curve)[1], coef(fit.curve)[2], coef(fit.curve)[3],
					coef(fit.curve)[4], coef(fit.curve)[5], coef(fit.curve)[6])
				y.err <- y.fit - y.fix
				SSE <- sum(y.err ^ 2)
				R2.1 <- 1 - SSE / SSY
			}
					
			if(diffs) {
				y1id <- subset(data, data$Subject == id.nums.g1[id] & data$Group == groups[1] & data$Curve == 2)
				
				y.fix <- y1id[,col]
				
				fit.curve <- est.dgauss.curve(time.all, y.fix, rho.0, concave[1], cor = cor)
				
				if(is.null(fit.curve$fit)) {
					coef.2 <- rep(NA, 6)
					sdev.2 <- rep(NA, 6)
					sigma.2 <- NA
					cor.temp.2 <- NA
					R2.2 <- NA
				} else {
					cor.temp.2 <- fit.curve$cor
					fit.curve <- fit.curve$fit
				
					coef.2 <- coef(fit.curve)
					sdev.2 <- sqrt(diag(fit.curve$varBeta))
					sigma.2 <- fit.curve$sigma
							
					SSY <- sum((y.fix - mean(y.fix)) ^ 2)
					y.fit <- dgauss(time.all, coef(fit.curve)[1], coef(fit.curve)[2], coef(fit.curve)[3],
						coef(fit.curve)[4], coef(fit.curve)[5], coef(fit.curve)[6])
					y.err <- y.fit - y.fix
					SSE <- sum(y.err ^ 2)
					R2.2 <- 1 - SSE / SSY
				}
			}
			if(diffs) {
				out <- c(coef.1, sdev.1, sigma.1, R2.1, cor.temp.1, coef.2, sdev.2, sigma.2, R2.2, cor.temp.2)
			} else {
				out <- c(coef.1, sdev.1, sigma.1, R2.1, cor.temp.1)
			}
			out
			
		}
		
		coef.id1 <- for.out[,1:6]
		sdev.id1 <- for.out[,7:12]
		sigma.id1[,1] <- as.numeric(for.out[,13])
		R2.g1.1 <- as.numeric(for.out[,14])
		cor.1 <- as.numeric(for.out[,15])
		if(diffs) {
			coef.id3 <- for.out[,16:21]
			sdev.id3 <- for.out[,22:27]
			sigma.id3[,1] <- as.numeric(for.out[,28])
			R2.g1.2 <- as.numeric(for.out[,29])
			cor.3 <- as.numeric(for.out[,30])
		}

		for.out <- foreach(id = 1:N.g2, .combine = rbind) %dopar% {
			if(diffs) {
				y2id <- subset(data, data$Subject == id.nums.g2[id] & data$Group == groups[2] & data$Curve == 1)
			} else {
				y2id <- subset(data, data$Subject == id.nums.g2[id] & data$Group == groups[2])
			}
			
			y.fix <- y2id[,col]
			
			fit.curve <- est.dgauss.curve(time.all, y.fix, rho.0, concave[2], cor = cor)
			
			if(is.null(fit.curve$fit)) {
				coef.1 <- rep(NA, 6)
				sdev.1 <- rep(NA, 6)
				sigma.1 <- NA
				cor.temp.1 <- NA
				R2.1 <- NA
			} else {
				cor.temp.1 <- fit.curve$cor
				fit.curve <- fit.curve$fit
			
				coef.1 <- coef(fit.curve)
				sdev.1 <- sqrt(diag(fit.curve$varBeta))
				sigma.1 <- fit.curve$sigma
						
				SSY <- sum((y.fix - mean(y.fix)) ^ 2)
				y.fit <- dgauss(time.all, coef(fit.curve)[1], coef(fit.curve)[2], coef(fit.curve)[3],
					coef(fit.curve)[4], coef(fit.curve)[5], coef(fit.curve)[6])
				y.err <- y.fit - y.fix
				SSE <- sum(y.err ^ 2)
				R2.1 <- 1 - SSE / SSY
			}

			if(diffs) {
				y2id <- subset(data, data$Subject == id.nums.g2[id] & data$Group == groups[2] & data$Curve == 2)
				
				y.fix <- y2id[,col]
			
				fit.curve <- est.dgauss.curve(time.all, y.fix, rho.0, concave[2], cor = cor)
				
				if(is.null(fit.curve$fit)) {
					coef.2 <- rep(NA, 6)
					sdev.2 <- rep(NA, 6)
					sigma.2 <- NA
					cor.temp.2 <- NA
					R2.2 <- NA
				} else {
					cor.temp.2 <- fit.curve$cor
					fit.curve <- fit.curve$fit
				
					coef.2 <- coef(fit.curve)
					sdev.2 <- sqrt(diag(fit.curve$varBeta))
					sigma.2 <- fit.curve$sigma
							
					SSY <- sum((y.fix - mean(y.fix)) ^ 2)
					y.fit <- dgauss(time.all, coef(fit.curve)[1], coef(fit.curve)[2], coef(fit.curve)[3],
						coef(fit.curve)[4], coef(fit.curve)[5], coef(fit.curve)[6])
					y.err <- y.fit - y.fix
					SSE <- sum(y.err ^ 2)
					R2.2 <- 1 - SSE / SSY
				}
			}
			
			if(diffs) {
				out <- c(coef.1, sdev.1, sigma.1, R2.1, cor.temp.1, coef.2, sdev.2, sigma.2, R2.2, cor.temp.2)
			} else {
				out <- c(coef.1, sdev.1, sigma.1, R2.1, cor.temp.1)
			}
			out
		}
		
		coef.id2 <- for.out[,1:6]
		sdev.id2 <- for.out[,7:12]
		sigma.id2[,1] <- as.numeric(for.out[,13])
		R2.g2.1 <- as.numeric(for.out[,14])
		cor.2 <- as.numeric(for.out[,15])
		if(diffs) {
			coef.id4 <- for.out[,16:21]
			sdev.id4 <- for.out[,22:27]
			sigma.id4[,1] <- as.numeric(for.out[,28])
			R2.g2.2 <- as.numeric(for.out[,29])
			cor.4 <- as.numeric(for.out[,30])
		}
		
		stopCluster(cl)
	}
	
	nofit <- sum(is.na(cor.1)) + sum(is.na(cor.2)) + sum(is.na(cor.3)) + sum(is.na(cor.4))
	if(nofit > 0) warning(paste("There were", nofit, "eyetracks that could not be fit"))
	if(cor == TRUE) {
		nocor <- sum(!cor.1, na.rm = TRUE) + sum(!cor.2, na.rm = TRUE) + sum(!cor.3, na.rm = TRUE) + sum(!cor.4, na.rm = TRUE)
		if(nocor > 0) warning(paste("There were", nocor, "eyetracks that were fit without the AR1 assumption"))
	}
	
	ar1.good <- sum(c(R2.g1.1, R2.g1.2, R2.g2.1, R2.g2.2) >= .95 &
		c(cor.1, cor.3, cor.2, cor.4), na.rm = TRUE)
	ar1.ok <- sum(c(R2.g1.1, R2.g1.2, R2.g2.1, R2.g2.2) < .95 &
		c(R2.g1.1, R2.g1.2, R2.g2.1, R2.g2.2) >= .8 &
		c(cor.1, cor.3, cor.2, cor.4), na.rm = TRUE)
	ar1.bad <- sum(c(R2.g1.1, R2.g1.2, R2.g2.1, R2.g2.2) < .8 &
		c(R2.g1.1, R2.g1.2, R2.g2.1, R2.g2.2) > 0 &
		c(cor.1, cor.3, cor.2, cor.4), na.rm = TRUE)
		
	nonar1.good <- sum(c(R2.g1.1, R2.g1.2, R2.g2.1, R2.g2.2) >= .95 &
		!c(cor.1, cor.3, cor.2, cor.4), na.rm = TRUE)
	nonar1.ok <- sum(c(R2.g1.1, R2.g1.2, R2.g2.1, R2.g2.2) < .95 &
		c(R2.g1.1, R2.g1.2, R2.g2.1, R2.g2.2) >= .8 &
		!c(cor.1, cor.3, cor.2, cor.4), na.rm = TRUE)
	nonar1.bad <- sum(c(R2.g1.1, R2.g1.2, R2.g2.1, R2.g2.2) < .8 &
		c(R2.g1.1, R2.g1.2, R2.g2.1, R2.g2.2) > 0 &
		!c(cor.1, cor.3, cor.2, cor.4), na.rm = TRUE)
		
	if(diffs) {
		nofit <- sum(is.na(c(cor.1, cor.3, cor.2, cor.4)))
	} else {
		nofit <- sum(is.na(c(cor.1, cor.2)))
	}

	cat("########################################\n")
	cat("############### FITS ###################\n")
	cat("########################################\n")
	cat(paste("AR1,       R2>=0.95   --", ar1.good, "\n"))
	cat(paste("AR1,     0.95>R2>=0.8 --", ar1.ok, "\n"))
	cat(paste("AR1,       0.8>R2     --", ar1.bad, "\n"))
	cat(paste("Non-AR1,   R2>=0.95   --", nonar1.good, "\n"))
	cat(paste("Non-AR1, 0.95>R2>=0.8 --", nonar1.ok, "\n"))
	cat(paste("Non-AR1,   0.8>R2     --", nonar1.bad, "\n"))
	cat(paste("No Fit                --", nofit, "\n"))
	cat("########################################\n\n")
	
	cat("Next Steps: Check goodness of fits (ests.plot, subs.plot),\n\trefit bad fits (doubleGauss.refit),\n\tbootstrap and t-test (doubleGauss.boot)\n\n")
  
	names(cor.1) <- 1:length(cor.1)
	names(cor.2) <- 1:length(cor.2)
	names(R2.g1.1) <- 1:length(R2.g1.1)
	names(R2.g2.1) <- 1:length(R2.g2.1)
	if(diffs) {
		names(cor.3) <- 1:length(cor.3)
		names(cor.4) <- 1:length(cor.4)
		names(R2.g1.2) <- 1:length(R2.g1.2)
		names(R2.g2.2) <- 1:length(R2.g2.2)
	}
	
  list(data = data, col = col, rho.0 = rho.0, N.time = N.time, N.sub1 = N.sub1, N.sub2 = N.sub2,
    coef.id1 = coef.id1, coef.id2 = coef.id2, sdev.id1 = sdev.id1, sdev.id2 = sdev.id2, sigma.id1 = sigma.id1, sigma.id2 = sigma.id2,
		coef.id3 = coef.id3, coef.id4 = coef.id4, sdev.id3 = sdev.id3, sdev.id4 = sdev.id4, sigma.id3 = sigma.id3, sigma.id4 = sigma.id4,
		id.nums.g1 = id.nums.g1, id.nums.g2 = id.nums.g2,
    groups = groups, time.all = time.all, N.g1 = N.g1, N.g2 = N.g2, concave = concave, model = "dgauss",
		R2.g1.1 = R2.g1.1, R2.g2.1 = R2.g2.1, R2.g1.2 = R2.g1.2, R2.g2.2 = R2.g2.2, diffs = diffs, cor = cor,
		cor.1 = cor.1, cor.2 = cor.2, cor.3 = cor.3, cor.4 = cor.4)
}
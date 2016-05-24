cinid.table <- function(HCW, mu4, sd4, threshold = 0.95, file = NULL, w = c(1, 1, 1, 1)) {
  
	HCW <- na.omit(as.vector(unlist(HCW)))
  
	## check validity
	if(any(c(HCW, mu4, sd4) < 0))
		stop("'HCW', 'mu4' and 'sd4' must be higher than zero")
	if((threshold < 0.25) | (threshold > 1))
		return("'threshold' must be higher than 0.25 and lower than 1")
  if(length(w) != 4)
    return("'w' must be a four-length vector")
  
	n <- length(HCW)
	instar_determ <- rep(NA, n)
	instar_stoch <- rep(NA, n)
  
	mu1 <- 0.348 * mu4 + 17.24
	mu2 <- 0.498 * mu4 + 11.65
	mu3 <- 0.693 * mu4 + 32.34
  
	sd1 <- 0.200 * sd4 + 13.63
	sd2 <- 0.421 * sd4 + 7.80
	sd3 <- 0.699 * sd4 + 3.90
  
	X1 <- (HCW - mu1) / sd1
	X2 <- (HCW - mu2) / sd2
	X3 <- (HCW - mu3) / sd3
	X4 <- (HCW - mu4) / sd4
  
	p1 <- 2 * (1 - pnorm(abs(X1)))
	p2 <- 2 * (1 - pnorm(abs(X2)))
	p3 <- 2 * (1 - pnorm(abs(X3)))
	p4 <- 2 * (1 - pnorm(abs(X4)))
	p <- cbind(p1, p2, p3, p4)
  
	d1 <- dnorm(abs(X1)) * w[1]
	d2 <- dnorm(abs(X2)) * w[2]
	d3 <- dnorm(abs(X3)) * w[3]
	d4 <- dnorm(abs(X4)) * w[4]
	d <- cbind(d1, d2, d3, d4)
	d_tot <- apply(d, 1, sum)
	d_freq <- d / d_tot
  
	for(i in 1:n)	{	
		if((d_tot[i] == 0) | (length(which(d_freq[i, ] == max(d_freq[i, ]))) > 1)) {
			instar_determ[i] <- "Indet"
			instar_stoch[i] <- "Indet"
		} else {
			if(max(d_freq[i, ]) < as.numeric(threshold))
				instar_determ[i] <- "Indet"
			else
				instar_determ[i] <- as.numeric(which(d[i, ] == max(d[i, ])))
			K <- as.null()
			K <- rmultinom(1, 1, d_freq[i, ])
			instar_stoch[i] <- which(K == 1)
		}
	}
  
	## output for each individuals
	table_indiv <- as.null()
	table_indiv <- as.data.frame(cbind(HCW, instar_determ, instar_stoch, round(p, 4), round(d_freq, 4)))
	colnames(table_indiv) <- c("HCW", "instar_determ", "instar_stoch", "p1", "p2", "p3", "p4", "rd1", "rd2", "rd3", "rd4")
  
	## output for the population
	number_instar_determ <- table(instar_determ)
	freq_instar_determ <- round(number_instar_determ / sum(number_instar_determ), 4)
	number_instar_stoch <- table(instar_stoch)
	freq_instar_stoch <- round(number_instar_stoch / sum(number_instar_stoch), 4)
	table_pop <- matrix(c(mu1, sd1, mu2, sd2, mu3, sd3, mu4, sd4, NA, NA), 2)
	colnames(table_pop) <- c("Instar1", "Instar2", "Instar3", "Instar4", "InstarIndet")
	table_pop <- rbind(table_pop, c(max(0, as.numeric(number_instar_determ[which(names(number_instar_determ) == 1)])), max(0, as.numeric(number_instar_determ[which(names(number_instar_determ) == 2)])), max(0, as.numeric(number_instar_determ[which(names(number_instar_determ) == 3)])), max(0, as.numeric(number_instar_determ[which(names(number_instar_determ) == 4)])), max(0, as.numeric(number_instar_determ[which(names(number_instar_determ) == "Indet")]))))
	table_pop <- rbind(table_pop, c(max(0, as.numeric(freq_instar_determ[which(names(freq_instar_determ) == 1)])), max(0, as.numeric(freq_instar_determ[which(names(freq_instar_determ) == 2)])), max(0, as.numeric(freq_instar_determ[which(names(freq_instar_determ) == 3)])), max(0, as.numeric(freq_instar_determ[which(names(freq_instar_determ) == 4)])), max(0, as.numeric(freq_instar_determ[which(names(freq_instar_determ) == "Indet")]))))
	table_pop <- rbind(table_pop, c(max(0, as.numeric(number_instar_stoch[which(names(number_instar_stoch) == 1)])), max(0, as.numeric(number_instar_stoch[which(names(number_instar_stoch) == 2)])), max(0, as.numeric(number_instar_stoch[which(names(number_instar_stoch) == 3)])), max(0, as.numeric(number_instar_stoch[which(names(number_instar_stoch) == 4)])), max(0, as.numeric(number_instar_stoch[which(names(number_instar_stoch) == "Indet")]))))
	table_pop <- rbind(table_pop, c(max(0, as.numeric(freq_instar_stoch[which(names(freq_instar_stoch) == 1)])), max(0, as.numeric(freq_instar_stoch[which(names(freq_instar_stoch) == 2)])), max(0, as.numeric(freq_instar_stoch[which(names(freq_instar_stoch) == 3)])), max(0, as.numeric(freq_instar_stoch[which(names(freq_instar_stoch) == 4)])), max(0, as.numeric(freq_instar_stoch[which(names(freq_instar_stoch) == "Indet")]))))
	rownames(table_pop) <- c("mu", "sd", "N_determ", "F_determ", "N_stoch", "F_stoch")
  
	if(is.null(file) == F) {
		write.table(table_indiv, paste(file, "_indiv.txt", sep = ""), quote = FALSE, row.names = FALSE)
		write.table(table_pop, paste(file, "_pop.txt", sep = ""), quote = FALSE, row.names = FALSE)
	}
  
	l <- new.env()
	l$indiv <- table_indiv
	l$pop <- table_pop
	l <- as.list(l)
	return(l)
}

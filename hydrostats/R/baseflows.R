baseflows <- function(flow.ts, a = 0.975, n.reflected = 30, ts = "mean") {
	
	full.flow.ts <- flow.ts[c("Date", "Q")]
	record.year <- strftime(full.flow.ts[["Date"]], format = "%Y")
	n.years <- length(levels(factor(record.year)))
	
	red.flow.ts <- full.flow.ts[complete.cases(full.flow.ts[["Date"]], full.flow.ts[["Q"]]), ]
	
	Date <- red.flow.ts[[1]]
	Q <- red.flow.ts[[2]]
	
	seq.length <- length(Q)
	
	Q <- c(Q[n.reflected:1], Q, Q[seq.length:(seq.length-n.reflected)])
	
	lh3 <- function(Q, a) {
		
		qb1 <- lh(Q, a)
		qb2 <- lh(rev(qb1), a)
		qb3 <- lh(rev(qb2), a)
		
		return(qb3)
		
	}
	
	
	lh <- function(Q, a) {
		
		qf <- rep(Q[1], length(Q))
		qb <- rep(Q[1], length(Q))
		
		qf[1] <- (Q[2] - Q[1]) * ((1 + a)/2)
		
		for (i in 2:length(Q)) {
			
			qf[i] <- a * qf[i - 1] + ((Q[i] - Q[i - 1]) * ((1 + a)/2))
			
		}
		for (i in 1:length(Q)) {
			
			if (qf[i] < 0) {
				qb[i] <- Q[i]
			} else if (Q[i] > qf[i]) {
				qb[i] <- Q[i] - qf[i]
			} else qb[i] <- 0
			
		}
		return(qb)
		
	}
	
	
	
	bf <- lh3(Q, a)
	bf <- bf[-c(1:n.reflected, (length(bf)-n.reflected):length(bf))]
	Q <- Q[-c(1:n.reflected, (length(bf)-n.reflected):length(bf))]
	bfi <- ifelse(Q == 0, 0, bf/Q)
	
	
	out <- data.frame(Date, bf, bfi)
	out <- merge(full.flow.ts, out, by = "Date", all.x = T)
	names(out) <- c("Date", "Q", "bf", "bfi")
	
	if (ts == "daily") {
		return(out)
	} else {
		
		if (ts == "annual") {
			a.obs <- aggregate(full.flow.ts$Q, by = list(year = strftime(full.flow.ts[["Date"]], format = "%Y")), function(x) {
				sum(!is.na(x))
			})
			names(a.obs) <- c("year", "no.obs")
			a.bf <- aggregate(out[2:4], by = list(year = strftime(out[["Date"]], format = "%Y")), mean, na.rm = T)
			out <- merge(a.obs, a.bf, by = "year", all.x = T)
			
			return(out)
			
		}
		
		
		if (ts == "mean") {
			all.days <- nrow(full.flow.ts)
			obs <- nrow(red.flow.ts)
			prop.obs = obs/all.days
			
			return(data.frame(n.years = n.years, prop.obs = prop.obs, MDF = mean(out[, 2], na.rm = T), Q50 = median(out[, 2], na.rm = T), mean.bf = mean(out[, 3], na.rm = T), mean.bfi = mean(out[, 
																																																																																														 3], na.rm = T)/mean(out[, 2], na.rm = T)))
			
		}
	}
} 

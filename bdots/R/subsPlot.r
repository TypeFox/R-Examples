subs.plot <- function(part1.list, legend.spot = "topright", ylim = NULL) {
	par(mfrow = c(1,1))
	
	data       <- part1.list$data
  col        <- part1.list$col
  N.time     <- part1.list$N.time
  N.sub1     <- part1.list$N.sub1
  N.sub2     <- part1.list$N.sub2
  coef.id1   <- part1.list$coef.id1
  coef.id2   <- part1.list$coef.id2
	coef.id3   <- part1.list$coef.id3
  coef.id4   <- part1.list$coef.id4
  id.nums.g1 <- part1.list$id.nums.g1
  id.nums.g2 <- part1.list$id.nums.g2
  groups     <- part1.list$groups
  time.all   <- part1.list$time.all
  N.g1       <- part1.list$N.g1
  N.g2       <- part1.list$N.g2
	model      <- part1.list$model
	diffs      <- part1.list$diffs
	
	if(is.null(ylim)) {
		ylim <- c(min(data[,col]), max(data[,col]))
	}
	xticks <- seq(0, max(time.all), round(max(time.all) / 10))
	yticks <- seq(ylim[1], ylim[2], round(diff(ylim) / 10, 2))
	gray.subs <- rgb(0, 0, 0, alpha = 20, maxColorValue = 255)
	
	dgauss <- function(time, mu, ht, sig1, sig2, base1, base2) {
		(time < mu) * (exp(-1 * (time - mu) ^ 2 / (2 * sig1 ^ 2)) * (ht - base1) +
		base1) + (mu <= time) * (exp(-1 * (time - mu) ^ 2 / (2 * sig2 ^ 2)) * (ht - base2) +
		base2)
	}
	
	logist <- function(time, mini, peak, slope, cross) {
		mini + (peak-mini) / (1 + exp(4 * slope * (cross - (time)) / (peak - mini)))
	}
	
	for(g in 1:2) {
		if(g == 1) {
			N <- N.g1
			id.nums <- id.nums.g1
			coef1 <- coef.id1
			coef2 <- coef.id3
		} else {
			N <- N.g2
			id.nums <- id.nums.g2
			coef1 <- coef.id2
			coef2 <- coef.id4
		}
		
		i <- 0
		while(i < N) {
			i <- i + 1
			if(i < 1) i <- 1
			if(diffs) {
				y1id <- subset(data, data$Subject == id.nums[i] & data$Group == groups[g] & data$Curve == 1)
			} else {
				y1id <- subset(data, data$Subject == id.nums[i] & data$Group == groups[g])
			}
			y.fix <- y1id[,col]
			
			if(model == "logistic") {
				mini  <- coef1[i,1]
				peak  <- coef1[i,2]
				slope <- coef1[i,3]
				cross <- coef1[i,4]
				y.fix.fitted <- logist(time.all, mini, peak, slope, cross)
			} else if (model == "dgauss") {
				mu    <- coef1[i,1]
				ht    <- coef1[i,2]
				sig1  <- coef1[i,3]
				sig2  <- coef1[i,4]
				base1 <- coef1[i,5]
				base2 <- coef1[i,6]
				y.fix.fitted <- dgauss(time.all, mu, ht, sig1, sig2, base1, base2)
			}

			plot(time.all, y.fix, type = "l", lwd = 2, xlab = "Time", ylab = "Proportion of Fixations",
				ylim = ylim, axes = FALSE)
			axis(1, at = xticks)
			axis(2, at = yticks)
			box()
			abline(v = xticks, h = yticks, col = gray.subs)
			if(diffs) {
				title(paste0("Group = ", groups[g], ", Subject = ", id.nums[i], ", Dataset # =", i, ", Curve = 1"))
			} else {
				title(paste0("Group = ", groups[g], ", Subject = ", id.nums[i], ", Dataset # =", i))
			}
			if(model == "logistic") {
				title(sub = paste0("Mini = ", round(mini, 2), ", Peak = ", round(peak, 2), ", Slope = ", round(slope, 4), ", Cross = ", round(cross)))
			} else if(model == "dgauss") {
				title(sub = paste0("Mu = ", round(mu), ", Height = ", round(ht, 2), ", Sigma1 = ", round(sig1), ", Sigma2 = ", round(sig2),
														", Base1 = ", round(base1, 3), ", Base2 = ", round(base2, 3)))
			}
			lines(time.all, y.fix.fitted, lty = 2, lwd = 2, col = "red")
			legend(legend.spot, lty = 1:2, col = c("black", "red"),
				legend = c("Observed", "Fitted"))
			
			read <- readline("Press 'Enter' to go to next plot  ")
			if(read == "back") {
				i <- i - 2
				next
			}
			
			if(diffs) {
				y1id <- subset(data, data$Subject == id.nums[i] & data$Group == groups[g] & data$Curve == 2)
				y.fix <- y1id[,col]
				
				if(model == "logistic") {
					mini  <- coef2[i,1]
					peak  <- coef2[i,2]
					slope <- coef2[i,3]
					cross <- coef2[i,4]
					y.fix.fitted <- logist(time.all, mini, peak, slope, cross)
				} else if (model == "dgauss") {
					mu    <- coef2[i,1]
					ht    <- coef2[i,2]
					sig1  <- coef2[i,3]
					sig2  <- coef2[i,4]
					base1 <- coef2[i,5]
					base2 <- coef2[i,6]
					y.fix.fitted <- dgauss(time.all, mu, ht, sig1, sig2, base1, base2)
				}

				plot(time.all, y.fix, type = "l", lwd = 2, xlab = "Time", ylab = "Proportion of Fixations",
					ylim = ylim, axes = FALSE)
				axis(1, at = xticks)
				axis(2, at = yticks)
				box()
				abline(v = xticks, h = yticks, col = gray.subs)
				title(paste0("Group = ", groups[g], ", Subject = ", id.nums[i], ", Dataset # =", i, ", Curve = 2"))
				if(model == "logistic") {
					title(sub = paste0("Mini = ", round(mini, 2), ", Peak = ", round(peak, 2), ", Slope = ", round(slope, 4), ", Cross = ", round(cross)))
				} else if(model == "dgauss") {
					title(sub = paste0("Mu = ", round(mu), ", Height = ", round(ht, 2), ", Sigma1 = ", round(sig1), ", Sigma2 = ", round(sig2),
														", Base1 = ", round(base1, 3), ", Base2 = ", round(base2, 3)))
				}
				lines(time.all, y.fix.fitted, lty = 2, lwd = 2, col = "red")
				legend(legend.spot, lty = 1:2, col = c("black", "red"),
					legend = c("Observed", "Fitted"))
					
				read <- readline("Press 'Enter' to go to next plot  ")
				if(read == "back") {
					i <- i - 2
					next
				}
			}
		}
	}
}

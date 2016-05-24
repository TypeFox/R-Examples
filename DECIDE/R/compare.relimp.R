`compare.relimp` <-
function(dataset1, dataset2) {

	ci.normal <- function(f, se.f) {
    		c(f - (1.96 * se.f), f + (1.96 * se.f))           
    		}

# relative importance of secondary effects - average

	ri1 <- relative.importance(dataset1)
	ri2 <- relative.importance(dataset2)

	ci.diff.ri.avg <- matrix(NA, nrow=ri1$no_classes, ncol=ri1$no_classes)
	ci.diff.ri.avg <- as.data.frame(ci.diff.ri.avg)

	ci.diff.ri.avg.lower <- matrix(NA, nrow=ri1$no_classes, ncol=ri1$no_classes)
	ci.diff.ri.avg.upper <- matrix(NA, nrow=ri1$no_classes, ncol=ri1$no_classes)


	for (i in 1:dim(ri1$rel_imp_sec_avg)[1]) {
		for (j in 1:dim(ri1$rel_imp_sec_avg)[2]) {
			if (i < j) {
				ci.diff.ri.avg.lower[i, j] <- ci.normal(ri1$rel_imp_sec_avg[i, j] - ri2$rel_imp_sec_avg[i, j],
								sqrt(((ri1$se.ri.avg[i, j])^2) + ((ri2$se.ri.avg[i, j])^2)))[1]
				ci.diff.ri.avg.upper[i, j] <- ci.normal(ri1$rel_imp_sec_avg[i, j] - ri2$rel_imp_sec_avg[i, j],
								sqrt(((ri1$se.ri.avg[i, j])^2) + ((ri2$se.ri.avg[i, j])^2)))[2]
	    			ci.diff.ri.avg[i, j] <- paste("(", format(ci.diff.ri.avg.lower[i, j], digits=5), ", ",
							 format(ci.diff.ri.avg.upper[i, j], digits=5), ")", sep="")
            			}
        		}
    		}
		
# relative importance - 1			

	ci.diff.ri.1 <- matrix(NA, nrow=ri1$no_classes, ncol=ri1$no_classes)
	ci.diff.ri.1 <- as.data.frame(ci.diff.ri.1)

	ci.diff.ri.1.lower <- matrix(NA, nrow=ri1$no_classes, ncol=ri1$no_classes)
	ci.diff.ri.1.upper <- matrix(NA, nrow=ri1$no_classes, ncol=ri1$no_classes)


	for (i in 1:dim(ri1$rel_imp_sec1)[1]) {
		for (j in 1:dim(ri1$rel_imp_sec1)[2]) {
			if (i < j) {
				ci.diff.ri.1.lower[i, j] <- ci.normal(ri1$rel_imp_sec1[i, j] - ri2$rel_imp_sec1[i, j],
								sqrt(((ri1$se.ri.1[i, j])^2) + ((ri2$se.ri.1[i, j])^2)))[1]
				ci.diff.ri.1.upper[i, j] <- ci.normal(ri1$rel_imp_sec1[i, j] - ri2$rel_imp_sec1[i, j],
								sqrt(((ri1$se.ri.1[i, j])^2) + ((ri2$se.ri.1[i, j])^2)))[2]
	    			ci.diff.ri.1[i, j] <- paste("(", format(ci.diff.ri.1.lower[i, j], digits=5), ", ",
							 format(ci.diff.ri.1.upper[i, j], digits=5), ")", sep="")
            			}
        		}
    		}
				
# relative importance - 2			

	ci.diff.ri.2 <- matrix(NA, nrow=ri1$no_classes, ncol=ri1$no_classes)
	ci.diff.ri.2 <- as.data.frame(ci.diff.ri.1)

	ci.diff.ri.2.lower <- matrix(NA, nrow=ri1$no_classes, ncol=ri1$no_classes)
	ci.diff.ri.2.upper <- matrix(NA, nrow=ri1$no_classes, ncol=ri1$no_classes)


	for (i in 1:dim(ri1$rel_imp_sec2)[1]) {
		for (j in 1:dim(ri1$rel_imp_sec2)[2]) {
			if (i < j) {
				ci.diff.ri.2.lower[i, j] <- ci.normal(ri1$rel_imp_sec2[i, j] - ri2$rel_imp_sec2[i, j],
								sqrt(((ri1$se.ri.2[i, j])^2) + ((ri2$se.ri.2[i, j])^2)))[1]
				ci.diff.ri.2.upper[i, j] <- ci.normal(ri1$rel_imp_sec2[i, j] - ri2$rel_imp_sec2[i, j],
								sqrt(((ri1$se.ri.2[i, j])^2) + ((ri2$se.ri.2[i, j])^2)))[2]
	    			ci.diff.ri.2[i, j] <- paste("(", format(ci.diff.ri.2.lower[i, j], digits=5), ", ",
							 format(ci.diff.ri.2.upper[i, j], digits=5), ")", sep="")
            			}
        		}
    		}
				


# differences in log odds

	ci.diff.lo.lower <- matrix(NA, nrow=ri1$no_classes, ncol=ri1$no_classes)
	ci.diff.lo.upper <- matrix(NA, nrow=ri1$no_classes, ncol=ri1$no_classes)

	ci.diff.lo <- matrix(NA, nrow=ri1$no_classes, ncol=ri1$no_classes)	
	ci.diff.lo <- as.data.frame(ci.diff.lo)

	for (i in 1:dim(ri1$log_odds)[1]) {
		for (j in 1:dim(ri1$log_odds)[2]) {
			ci.diff.lo.lower[i, j] <- ci.normal(ri1$log_odds[i, j] - ri2$log_odds[i, j],
								sqrt(((ri1$se_logodds[i, j])^2) + ((ri2$se_logodds[i, j])^2)))[1]
			ci.diff.lo.upper[i, j] <- ci.normal(ri1$log_odds[i, j] - ri2$log_odds[i, j],
								sqrt(((ri1$se_logodds[i, j])^2) + ((ri2$se_logodds[i, j])^2)))[2]
				
	    		ci.diff.lo[i, j] <- paste("(", format(ci.diff.lo.lower[i, j], digits=5), ", ",
							 format(ci.diff.lo.upper[i, j], digits=5), ")", sep="")
        		}
    		}
				

# confidence intervals for differences in log odds ratios

	ci.diff.lor.lower <- matrix(NA, nrow=ri1$no_classes, ncol=ri1$no_classes)
	ci.diff.lor.upper <- matrix(NA, nrow=ri1$no_classes, ncol=ri1$no_classes)

	ci.diff.lor <- matrix(NA, nrow=ri1$no_classes, ncol=ri1$no_classes)	
	ci.diff.lor <- as.data.frame(ci.diff.lor)

	for (i in 1:dim(ri1$log_oddsratios)[1]) {
		for (j in 1:dim(ri1$log_oddsratios)[2]) {
			if (i < j) {
				ci.diff.lor.lower[i, j] <- ci.normal(ri1$log_oddsratios[i, j] - ri2$log_oddsratios[i, j],
								sqrt(((ri1$se_logoddsratios[i, j])^2) + ((ri2$se_logoddsratios[i, j])^2)))[1]
				ci.diff.lor.upper[i, j] <- ci.normal(ri1$log_oddsratios[i, j] - ri2$log_oddsratios[i, j],
								sqrt(((ri1$se_logoddsratios[i, j])^2) + ((ri2$se_logoddsratios[i, j])^2)))[2]
				
	    			ci.diff.lor[i, j] <- paste("(", format(ci.diff.lor.lower[i, j], digits=5), ", ", format(ci.diff.lor.upper[i, j], digits=5), ")", sep="")
        			}
    			}
		}		



# significance testing


# > chisq.3df.a
# [1] 7.814728


# 95% quantile of chi-square distribution with df=no_classes

	chisq.ndf.a <- qchisq(p=0.95, df=ri1$no_classes, lower.tail = TRUE)	


# differences in log odds

	lo.diff <- matrix(NA, nrow=1, ncol=ri1$no_classes)

	for (i in 1:ri1$no_classes) {
		lo.diff[i] <- ri1$log_odds[i, i] - ri2$log_odds[i, i]
		}

	lo.diff.cov.mat <- matrix(rep(0, (ri1$no_classes)^2), nrow=ri1$no_classes, ncol=ri1$no_classes)

	for (i in 1:ri1$no_classes) {
		lo.diff.cov.mat[i, i] <- (ri1$se_logodds[i, i]^2) + (ri2$se_logodds[i, i]^2)
		}

	lo.diff.cov.mat.inv <- solve(lo.diff.cov.mat)

   # test statistic for differences in log odds

	test.diff.lo <- as.numeric(lo.diff %*% lo.diff.cov.mat.inv %*% t(lo.diff))

   # p-value

	test.diff.lo.pvalue <- pchisq(test.diff.lo, df=ri1$no_classes, lower.tail=FALSE)



# differences in log odds ratios

	no_oddsratios <- function(no_classes) {
    		N <- no_classes
    		no_oddsratios <- 0
		for (j in 1:(N-1)) {
        		no_oddsratios <- no_oddsratios + j
        		}
    		no_oddsratios
    		}

	lor.diff.mat <- matrix(NA, nrow=ri1$no_classes, ncol=ri1$no_classes)

	for (i in 1:dim(ri1$log_oddsratios)[1]) {
		for (j in 1:dim(ri1$log_oddsratios)[1]) {
			if (i < j) {
				lor.diff.mat[i, j] <- ri1$log_oddsratios[i, j] - ri2$log_oddsratios[i, j]
				}
			}
		}

	lor.diff <- matrix(NA, nrow=1, ncol=no_oddsratios(ri1$no_classes))

	counter1 <- 1
	for (i in 1:dim(ri1$log_oddsratios)[1]) {
		for (j in 1:dim(ri1$log_oddsratios)[1]) {
			if (!is.na(lor.diff.mat[i, j])) {
				lor.diff[counter1] <- lor.diff.mat[i, j]
				counter1 <- counter1 + 1
				}
			}
		}

# variance of log odds ratios of dataset 1

	lor.var.1 <- matrix(NA, nrow=1, ncol=no_oddsratios(ri1$no_classes))

	counter2 <- 1
	for (i in 1:dim(ri1$se_logoddsratios)[1]) {
		for (j in 1:dim(ri1$se_logoddsratios)[1]) {
			if (!is.na(ri1$se_logoddsratios[i, j])) {
				lor.var.1[1, counter2] <- (ri1$se_logoddsratios[i, j])^2
				counter2 <- counter2 + 1
				}
			}
		}

# variance of log odds ratios of dataset 2

	lor.var.2 <- matrix(NA, nrow=1, ncol=no_oddsratios(ri2$no_classes))

	temp1 <- matrix(NA, nrow=2, ncol=no_oddsratios(ri1$no_classes))

	counter3 <- 1
	for (i in 1:dim(ri2$se_logoddsratios)[1]) {
		for (j in 1:dim(ri2$se_logoddsratios)[1]) {
			if (!is.na(ri2$se_logoddsratios[i, j])) {
				lor.var.2[1, counter3] <- (ri2$se_logoddsratios[i, j])^2
				temp1[1, counter3] <- i
				temp1[2, counter3] <- j
				counter3 <- counter3 + 1
				}
			}
		}
					
	lor.diff.cov.mat <- matrix(rep(0, (no_oddsratios(ri1$no_classes))^2), nrow=no_oddsratios(ri1$no_classes), 
										ncol=no_oddsratios(ri1$no_classes))

	for (i in 1:no_oddsratios(ri1$no_classes)) {
		lor.diff.cov.mat[i, i] <- lor.var.1[i] + lor.var.2[i]
		}

	temp2 <- as.numeric(NA)

	for (i in 1:no_oddsratios(ri1$no_classes)) {
		for (j in 1:no_oddsratios(ri1$no_classes)) {
			if (i != j) {	
				if (temp1[1, i]==temp1[1, j]) {
					temp2 <- temp1[1, i]
					lor.diff.cov.mat[i, j] <- (ri1$se_logodds[temp2, temp2])^2 + (ri2$se_logodds[temp2, temp2])^2
					}
				if (temp1[2, i]==temp1[2, j]) {
					temp2 <- temp1[2, i]
					lor.diff.cov.mat[i, j] <- (ri1$se_logodds[temp2, temp2])^2 + (ri2$se_logodds[temp2, temp2])^2
					}
				if (temp1[2, i]==temp1[1, j]) {
					temp2 <- temp1[2, i]
					lor.diff.cov.mat[i, j] <- (ri1$se_logodds[temp2, temp2])^2 + (ri2$se_logodds[temp2, temp2])^2
					}
				if (temp1[1, i]==temp1[2, j]) {
					temp2 <- temp1[1, i]
					lor.diff.cov.mat[i, j] <- (ri1$se_logodds[temp2, temp2])^2 + (ri2$se_logodds[temp2, temp2])^2
					}
				}
			}
		}

	lor.diff.cov.mat.inv <- solve(lor.diff.cov.mat)

   # test statistic for differences in log odds ratios

	test.diff.lor <- as.numeric(lor.diff %*% lor.diff.cov.mat.inv %*% t(lor.diff))

   # p-value

	test.diff.lor.pvalue <- pchisq(test.diff.lor, df=no_oddsratios(ri1$no_classes), lower.tail=FALSE)


return(list("ci.diff.lo"=ci.diff.lo, "test.diff.lo"=test.diff.lo, "test.diff.lo.pvalue"=test.diff.lo.pvalue, "ci.diff.lor"=ci.diff.lor, "test.diff.lor"=test.diff.lor, "test.diff.lor.pvalue"=test.diff.lor.pvalue, "ci.diff.ri.1"=ci.diff.ri.1, "ci.diff.ri.2"=ci.diff.ri.2, "ci.diff.ri.avg"=ci.diff.ri.avg))
}


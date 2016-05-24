`relative.importance` <-
function(dataset) {

	dataset <- prepare.data(dataset)

	# how many classes are there?
	no_classes <- length(levels(dataset$class))

	data_class <- create.classdata(dataset)

	# sample size
	n <- rep(as.numeric(NA), each=no_classes)

	# n[i]: sample size for class i

	for (i in 1:no_classes) {
    		n[i] <- dim(data_class[[i]])[1]
    		}

	n.total <- dim(dataset)[1]

# initialise data frame with parameters for each class

	parameters <- matrix(data = NA, nrow = no_classes, ncol = 4)
	parameters <- as.data.frame(parameters)
	names(parameters) <- c("alpha", "beta", "mu", "sigma")

# run logistic regression for each class

	fit <- list()

	for (i in 1:no_classes) {
    		fit[[i]] <- glm(transition ~ performance, data=data_class[[i]], family="binomial")
    		}

	for (i in 1:no_classes) {
   		 parameters$alpha[i] <- as.numeric(fit[[i]]$coefficients[1])
  		 parameters$beta[i] <- as.numeric(fit[[i]]$coefficients[2])
    		 parameters$mu[i] <- mean(data_class[[i]]$performance)
   		 parameters$sigma[i] <- sd(data_class[[i]]$performance)
    		 }

####	relative.importance$parameters <- parameters


# percentage of each class that made the transition

	perc.class <- rep(as.numeric(NA), no_classes)

	for (i in 1:no_classes) {
    		temp1 <- data_class[[i]]
    		temp2 <- dim(temp1[(temp1$transition==1),])[1]
    		perc.class[i] <- temp2 / (dim(temp1)[1])
    		}

# total percentage that made the transition

	perc.total <- (dim(dataset[(dataset$transition==1),])[1]) / (dim(dataset)[1])


# 50% point of transition (-a/b) for each class

	fifty.point <- rep(as.numeric(NA), no_classes)

	for (i in 1:no_classes) {
    		fifty.point[i] <- - parameters$alpha[i] / parameters$beta[i]
   		}


# gamma (g=a+bm) for each class

	gamma <- rep(as.numeric(NA), no_classes)

	for (i in 1:no_classes) {
    		gamma[i] <- parameters$alpha[i] + (parameters$beta[i] * parameters$mu[i])
    		}

	k <- 0.61

# function calculating log odds using formula 3

	log.odds <- function(g, s, b) as.numeric(g/sqrt(1+(k*s*b)^2))      


# log odds for each class
	log_odds <- rep(as.numeric(NA), no_classes)
	odds <- rep(as.numeric(NA), no_classes)

	for (i in 1:no_classes) {
    		log_odds[i] <- log.odds(gamma[i], parameters$sigma[i], parameters$beta[i])
    		}

	for (i in 1:no_classes) {
		odds[i] <- exp(log_odds[i])
    		}

# how many odds ratios can be formed given the number of classes (sum of i from 1 to N-1)

no_oddsratios <- function(no_classes) {
    N <- no_classes
    no_oddsratios <- 0

    for (j in 1:(N-1)) {
        no_oddsratios <- no_oddsratios + j
        }
    no_oddsratios
    }

# table of odds ratios

	oddsratios <- matrix(NA, ncol=no_classes, nrow=no_classes)

	for (i in 1:dim(oddsratios)[1]) {
    		for (j in 1:dim(oddsratios)[2]) {
        		if (i < j) {
            			oddsratios[i, j] <- odds[i] / odds[j]
            			}
        		}
    		}

##### save.image("C:/Users/sony/Desktop/primary&secondary/ycs9.RData")

### str()
### dev.new()

# logistic function
	logistic <- function(x) exp(x)/(1+exp(x))       

# logit function
	logit <- function(x) log(x/(1-x))


# transition probabilities

	trans.prob <- rep(as.numeric(NA), no_classes)

	for (i in 1:no_classes) {
		trans.prob[i] <- logistic(log_odds[i])
		}


# transition probabilities using eq. 2

	tr.pr <- function(g, s, b) {
		pnorm(as.numeric(k*g / sqrt(1 + (k*s*b)^2)))
		}

#	trans.prob1 <- rep(as.numeric(NA), no_classes)

#	for (i in 1:no_classes) {
#		trans.prob1[i] <- tr.pr(gamma[i], parameters$sigma[i], parameters$beta[i])
#		}
      

# number of counterfactual combinations for logodds: 2*choose(no_classes, 2)

	logodds <- matrix(NA, nrow = no_classes, ncol = no_classes)

	for (i in 1:no_classes) {
		for (j in 1:no_classes) {
			logodds[i, j] <- log.odds(parameters$alpha[j] + (parameters$beta[j] * parameters$mu[i]), parameters$sigma[i], parameters$beta[j])
			}
		}

# function log_odds not really needed anymore; log_odds are the diagonal elements of logodds


# transition probabilities (using eq. 2)

	transprob <- matrix(NA, nrow = no_classes, ncol = no_classes)

	for (i in 1:no_classes) {
		for (j in 1:no_classes) {
			transprob[i, j] <- tr.pr(parameters$alpha[j] + (parameters$beta[j] * parameters$mu[i]), parameters$sigma[i], parameters$beta[j])
			}
		}

# trans.prob1 not needed; they are the diagonal elements of transprob



# relative importance of secondary effects 1

# number of estimates for relative importance = no_oddsratios

	rel.imp.sec1 <- matrix(NA, ncol=no_classes, nrow=no_classes)

	for (i in 1:dim(rel.imp.sec1)[1]) {
    		for (j in 1:dim(rel.imp.sec1)[2]) {
        		if (i < j) {
            			rel.imp.sec1[i, j] <- (logodds[i, i] - logodds[i, j]) / (logodds[i, i] - logodds[j, j])
            			}
        		}
    		}


# relative importance of secondary effects 2

	rel.imp.sec2 <- matrix(NA, ncol=no_classes, nrow=no_classes)

	for (i in 1:dim(rel.imp.sec2)[1]) {
	    for (j in 1:dim(rel.imp.sec2)[2]) {
	        if (i < j) {
        	    rel.imp.sec2[i, j] <- (logodds[j, i] - logodds[j, j]) / (logodds[i, i] - logodds[j, j])
         		   }
       			}
    		}

# relative importance of secondary effects - average

	rel.imp.sec.avg <- matrix(NA, ncol=no_classes, nrow=no_classes)

	for (i in 1:dim(rel.imp.sec.avg)[1]) {
		for (j in 1:dim(rel.imp.sec.avg)[2]) {
        		if (i < j) {
            			rel.imp.sec.avg[i, j] <- (rel.imp.sec1[i, j] + rel.imp.sec2[i, j]) / 2
            			}
        		}
    		}


# relative importance of primary effects 1

	rel.imp.prim1 <- matrix(NA, ncol=no_classes, nrow=no_classes)

	for (i in 1:dim(rel.imp.prim1)[1]) {
    		for (j in 1:dim(rel.imp.prim1)[2]) {
        		if (i < j) {
            			rel.imp.prim1[i, j] <- 1 - rel.imp.sec1[i, j]
            			}
        		}
    		}


# relative importance of primary effects 2

	rel.imp.prim2 <- matrix(NA, ncol=no_classes, nrow=no_classes)

	for (i in 1:dim(rel.imp.prim2)[1]) {
    		for (j in 1:dim(rel.imp.prim2)[2]) {
        		if (i < j) {
            			rel.imp.prim2[i, j] <- 1 - rel.imp.sec2[i, j]
            			}
        		}
    		}


# relative importance of primary effects - average

	rel.imp.prim.avg <- matrix(NA, ncol=no_classes, nrow=no_classes)

	for (i in 1:dim(rel.imp.prim.avg)[1]) {
    		for (j in 1:dim(rel.imp.prim.avg)[2]) {
        		if (i < j) {
            			rel.imp.prim.avg[i, j] <- 1 - rel.imp.sec.avg[i, j]
            			}
        		}
    		}



# inverse of information matrix for each class

	I <- list()

	for (i in 1:no_classes) {
    		I[[i]] <- vcov(fit[[i]])
    		}



# function calculating variance of log odds

	var.logodds <- function(a, b, m, s, i11, i12, i22, n) {       
    		var1 <- i11/(1 + (k*s*b)^2)
    		var2 <- 2*i12*(m - (2*a*b*(k*s)^2))/((1+(k*s*b)^2)^2)
    		var3 <- i22*((m - (2*a*b*(k*s)^2))^2)/((1+(k*s*b)^2)^3)
    		var4 <- ((s*b)^2) / (n*(1 + (k*s*b)^2))
    		var5 <- (2*((k*b*s)^4)*((a + b*m)^2) / (n*((1+(k*s*b)^2)^3)))
    		as.numeric(var1+var2+var3+var4+var5)
    		}


# variance of log odds

	var.lo <- matrix(NA, ncol=no_classes, nrow=no_classes)

	for (i in 1:dim(var.lo)[1]) {
    		for (j in 1:dim(var.lo)[2]) {
			var.lo[i, j] <- var.logodds(parameters$alpha[j], parameters$beta[j], parameters$mu[i], parameters$sigma[i], I[[j]][1,1], I[[j]][1,2], I[[j]][2,2], n[i])
			}
		}



# standard error of log odds of transition by the delta method

	se.lo <- matrix(NA, ncol=no_classes, nrow=no_classes)

	for (i in 1:dim(var.lo)[1]) {
    		for (j in 1:dim(var.lo)[2]) {
			se.lo[i, j] <- sqrt(var.lo[i, j])
			}
		}

	

# function calculating covariance of log odds ii, ij

	cov.logodds.ii.ij <- function(a.i, a.j, b.i, b.j, m.i, s.i, n.i) {        
    		as.numeric((s.i^2*b.i*b.j/(n.i*sqrt(1+(k*s.i*b.i)^2)*sqrt(1+(k*s.i*b.j)^2)) + (s.i*k)^4*(b.i*b.j)^2*(a.i+b.i*m.i)*(a.j + 
			b.j*m.i)/(2*n.i*((1+(k*s.i*b.i)^2)^(3/2))*((1+(k*s.i*b.j)^2)^(3/2)))))
		}

	cov.lo.ii.ij <- matrix(NA, nrow=no_classes, ncol=no_classes)

	for (i in 1:dim(cov.lo.ii.ij)[1]) {
		for (j in 1:dim(cov.lo.ii.ij)[2]) {
			if (i != j) {
				cov.lo.ii.ij[i, j] <- cov.logodds.ii.ij(parameters$alpha[i], parameters$alpha[j], parameters$beta[i],
				parameters$beta[j], parameters$mu[i], parameters$sigma[i], n[i])
				}
			}
		}


# function calculating covariance of log odds ii, ji

	cov.logodds.ii.ji <- function(a.i, b.i, m.i, s.i, i11, i12, i22, m.j, s.j) {     
    		var1 <- i11/sqrt(1 + ((k*s.i*b.i)^2)) + ((m.i-(a.i*b.i*(k*s.i)^2))*i12/((1 + (k*s.i*b.i)^2)^(3/2)))
    		var2 <- i12/sqrt(1 + ((k*s.i*b.i)^2)) + ((m.i-(a.i*b.i*(k*s.i)^2))*i22/((1 + (k*s.i*b.i)^2)^(3/2)))
    		as.numeric((var1/sqrt(1 + (k*s.j*b.i)^2) ) + (var2 * ((m.j-a.i*b.i*(k*s.j)^2)/((1 + (k*s.j*b.i)^2))^(3/2))))
    		}

	cov.lo.ii.ji <- matrix(NA, nrow=no_classes, ncol=no_classes)

	for (i in 1:dim(cov.lo.ii.ji)[1]) {
		for (j in 1:dim(cov.lo.ii.ji)[2]) {
			if (i != j) {
				cov.lo.ii.ji[i, j] <- cov.logodds.ii.ji(parameters$alpha[i], parameters$beta[i], parameters$mu[i],
				 	parameters$sigma[i], I[[i]][1, 1], I[[i]][1, 2], I[[i]][2, 2],
				 	parameters$mu[j], parameters$sigma[j])
					}
			}
		}




# variance of relative importance 1

	var.relimp.1 <- function(fii, fij, fjj, var.fii, var.fij, var.fjj, cov.fii.fij, cov.fij.fjj) {
		as.numeric(((fij-fjj)^2*var.fii/((fii-fjj)^4)) + ((fii-fij)^2*var.fjj/((fii-fjj)^4)) + (var.fij/((fii-fjj)^2))
			- (2*(fij-fjj)*cov.fii.fij/((fii-fjj)^3)) - (2*(fii-fij)*cov.fij.fjj/((fii-fjj)^3)))
		}


	var.ri.1 <- matrix(NA, ncol=no_classes, nrow=no_classes)

	for (i in 1:dim(var.ri.1)[1]) {
    		for (j in 1:dim(var.ri.1)[2]) {
        		if (i < j) {
            			var.ri.1[i, j] <- var.relimp.1(logodds[i, i], logodds[i, j], logodds[j, j], var.lo[i, i], var.lo[i, j],
							var.lo[j, j], cov.lo.ii.ij[i, j], cov.lo.ii.ij[j, i])
            			}
        		}
    		}


# variance of relative importance 2

	var.relimp.2 <- function(fii, fji, fjj, var.fii, var.fji, var.fjj, cov.fii.fji, cov.fji.fjj) {
		as.numeric(((fji-fjj)^2*var.fii/((fii-fjj)^4)) + ((fii-fji)^2*var.fjj/((fii-fjj)^4)) + (var.fji/((fii-fjj)^2))
			- (2*(fji-fjj)*cov.fii.fji/((fii-fjj)^3)) - (2*(fii-fji)*cov.fji.fjj/((fii-fjj)^3)))
		}

	var.ri.2 <- matrix(NA, ncol=no_classes, nrow=no_classes)

	for (i in 1:dim(var.ri.2)[1]) {
    		for (j in 1:dim(var.ri.2)[2]) {
        		if (i < j) {
            			var.ri.2[i, j] <- var.relimp.2(logodds[i, i], logodds[j, i], logodds[j, j], var.lo[i, i], var.lo[j, i],
							var.lo[j, j], cov.lo.ii.ji[i, j], cov.lo.ii.ji[j, i])
            			}
        		}
    		}

# standard error of relative importance 1

	se.ri.1 <- matrix(NA, ncol=no_classes, nrow=no_classes)

	for (i in 1:dim(se.ri.1)[1]) {
    		for (j in 1:dim(se.ri.1)[2]) {
        		if (i < j) {
            			se.ri.1[i, j] <- sqrt(var.ri.1[i, j])
            			}
        		}
    		}

# standard error of relative importance 2

	se.ri.2 <- matrix(NA, ncol=no_classes, nrow=no_classes)

	for (i in 1:dim(se.ri.2)[1]) {
    		for (j in 1:dim(se.ri.2)[2]) {
        		if (i < j) {
            			se.ri.2[i, j] <- sqrt(var.ri.2[i, j])
            			}
        		}
    		}

# > qnorm(1-0.025)
# [1] 1.959964


# function calculating 95% confidence intervals

	ci.normal <- function(f, se.f) {
    		c(f - (1.96 * se.f), f + (1.96 * se.f))           
    		}

# CI for log odds

	ci.lo <- matrix(NA, nrow = no_classes, ncol = 2)

	for (i in 1:no_classes) {
		ci.lo[i, 1] <- ci.normal(logodds[i, i], se.lo[i, i])[1]
		ci.lo[i, 2] <- ci.normal(logodds[i, i], se.lo[i, i])[2]	
		}

# log odds ratios

	log.oddsratios <- matrix(NA, ncol=no_classes, nrow=no_classes)

	for (i in 1:dim(log.oddsratios)[1]) {
    		for (j in 1:dim(log.oddsratios)[2]) {
        		if (i < j) {
            			log.oddsratios[i, j] <- logodds[i, i] - logodds[j, j]
            			}
        		}
    		}


# standard errors for log odds ratios

	se.log.oddsratios <- matrix(NA, ncol=no_classes, nrow=no_classes)

	for (i in 1:dim(se.log.oddsratios)[1]) {
    		for (j in 1:dim(se.log.oddsratios)[2]) {
        		if (i < j) {
            			se.log.oddsratios[i, j] <- sqrt(var.lo[i, i] + var.lo[j, j])
            			}
        		}
    		}


# confidence intervals for log odds ratios

	ci.log.oddsratios <- matrix(NA, ncol=no_classes, nrow=no_classes)
	ci.log.oddsratios <- as.data.frame(ci.log.oddsratios)

	ci.log.oddsratios.lower <- matrix(NA, ncol=no_classes, nrow=no_classes)
	ci.log.oddsratios.upper <- matrix(NA, ncol=no_classes, nrow=no_classes)

	for (i in 1:dim(ci.log.oddsratios)[1]) {
    		for (j in 1:dim(ci.log.oddsratios)[2]) {
        		if (i < j) {
            			ci.log.oddsratios.lower[i, j] <- ci.normal(log.oddsratios[i, j], se.log.oddsratios[i, j])[1]
            			ci.log.oddsratios.upper[i, j] <- ci.normal(log.oddsratios[i, j], se.log.oddsratios[i, j])[2]
            			ci.log.oddsratios[i, j] <- paste("(", format(ci.log.oddsratios.lower[i, j], digits=5), ", ", 
				format(ci.log.oddsratios.upper[i, j], digits=5), ")", sep="")
            			}
        		}
    		}


# CI for relative importance of secondary effects 1

	ci.ri.1 <- matrix(NA, nrow=no_classes, ncol=no_classes)
	ci.ri.1 <- as.data.frame(ci.ri.1)

	ci.ri.1.lower <- matrix(NA, nrow=no_classes, ncol=no_classes)
	ci.ri.1.upper <- matrix(NA, nrow=no_classes, ncol=no_classes)

	for (i in 1:dim(ci.ri.1)[1]) {
    		for (j in 1:dim(ci.ri.1)[2]) {
			if (i < j) {
            			ci.ri.1.lower[i, j] <- ci.normal(rel.imp.sec1[i, j], se.ri.1[i, j])[1]
	    			ci.ri.1.upper[i, j] <- ci.normal(rel.imp.sec1[i, j], se.ri.1[i, j])[2]
	    			ci.ri.1[i, j] <- paste("(", format(ci.ri.1.lower[i, j], digits=5), ", ", format(ci.ri.1.upper[i, j], digits=5), ")", sep="")
            			}
        		}
    		}

# CI for relative importance of secondary effects 2

	ci.ri.2 <- matrix(NA, nrow=no_classes, ncol=no_classes)
	ci.ri.2 <- as.data.frame(ci.ri.2)

	ci.ri.2.lower <- matrix(NA, nrow=no_classes, ncol=no_classes)
	ci.ri.2.upper <- matrix(NA, nrow=no_classes, ncol=no_classes)

	for (i in 1:dim(ci.ri.2)[1]) {
    		for (j in 1:dim(ci.ri.2)[2]) {
			if (i < j) {
            			ci.ri.2.lower[i, j] <- ci.normal(rel.imp.sec2[i, j], se.ri.2[i, j])[1]
	    			ci.ri.2.upper[i, j] <- ci.normal(rel.imp.sec2[i, j], se.ri.2[i, j])[2]
	    			ci.ri.2[i, j] <- paste("(", format(ci.ri.2.lower[i, j], digits=5), ", ", format(ci.ri.2.upper[i, j], digits=5), ")", sep="")
            			}
        		}
    		}



# function that calculates the variance of the average of the two estimates of relative importance

	var.relimp.avg <- function(fji, fij, fii, fjj, var.fji, var.fij, var.fii, var.fjj, cov.fji.fjj, cov.fij.fii, cov.fji.fii, cov.fij.fjj) {
    			as.numeric(((var.fji + var.fij)/(4*((fii-fjj)^2))) + ((fji-fij)^2*(var.fii+var.fjj)/(4*((fii-fjj)^4)))
    				+ ((fji-fij)*(cov.fji.fjj+cov.fij.fii-cov.fji.fii-cov.fij.fjj)/(2*((fii-fjj)^3))))
			}


	var.ri.avg <- matrix(NA, ncol=no_classes, nrow=no_classes)

	for (i in 1:dim(var.ri.avg)[1]) {
    		for (j in 1:dim(var.ri.avg)[2]) {
        		if (i < j) {
            			var.ri.avg[i, j] <- var.relimp.avg(logodds[j, i], logodds[i, j], logodds[i, i], logodds[j, j], var.lo[j, i],
						var.lo[i, j], var.lo[i, i], var.lo[j, j], cov.lo.ii.ij[j, i], cov.lo.ii.ij[i, j],
						cov.lo.ii.ji[i, j], cov.lo.ii.ji[j, i])
            			}
        		}
    		}



# standard error of relative importance - average

	se.ri.avg <- matrix(NA, ncol=no_classes, nrow=no_classes)

	for (i in 1:dim(se.ri.avg)[1]) {
    		for (j in 1:dim(se.ri.avg)[2]) {
        		if (i < j) {
            			se.ri.avg[i, j] <- sqrt(var.ri.avg[i, j])
            			}
        		}
    		}


# CI for relative importance of secondary effects - average

	ci.ri.avg <- matrix(NA, nrow=no_classes, ncol=no_classes)
	ci.ri.avg <- as.data.frame(ci.ri.avg)

	ci.ri.avg.lower <- matrix(NA, nrow=no_classes, ncol=no_classes)
	ci.ri.avg.upper <- matrix(NA, nrow=no_classes, ncol=no_classes)

	for (i in 1:dim(ci.ri.avg)[1]) {
    		for (j in 1:dim(ci.ri.avg)[2]) {
			if (i < j) {
            			ci.ri.avg.lower[i, j] <- ci.normal(rel.imp.sec.avg[i, j], se.ri.avg[i, j])[1]
	    			ci.ri.avg.upper[i, j] <- ci.normal(rel.imp.sec.avg[i, j], se.ri.avg[i, j])[2]
	    			ci.ri.avg[i, j] <- paste("(", format(ci.ri.avg.lower[i, j], digits=5), ", ", format(ci.ri.avg.upper[i, j], digits=5), ")", sep="")
            			}
        		}
    		}

return(list("sample_size"=n.total, "no_classes"=no_classes, "class_size"=n, "percentage_overall"=perc.total,
	"percentage_class"=perc.class, "fifty_point"=fifty.point, "parameters"=parameters, "transition_prob"=transprob, 
	"log_odds"=logodds, "se_logodds"=se.lo, "ci_logodds"=ci.lo, "odds"=odds, "log_oddsratios"=log.oddsratios, 
	"se_logoddsratios"=se.log.oddsratios, "ci_logoddsratios"=ci.log.oddsratios, "oddsratios"=oddsratios, 
	"rel_imp_prim1"=rel.imp.prim1, "rel_imp_prim2"=rel.imp.prim2, "rel_imp_prim_avg"=rel.imp.prim.avg,
	"rel_imp_sec1"=rel.imp.sec1, "rel_imp_sec2"=rel.imp.sec2, "rel_imp_sec_avg"=rel.imp.sec.avg,
	"se.ri.1"=se.ri.1, "ci.ri.1"=ci.ri.1, "se.ri.2"=se.ri.2, "ci.ri.2"=ci.ri.2, "se.ri.avg"=se.ri.avg, "ci.ri.avg"=ci.ri.avg))

}


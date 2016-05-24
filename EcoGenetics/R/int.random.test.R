#' random test
#' 
#' @author Leandro Roser \email{leandroroser@@ege.fcen.uba.ar}
#' @keywords internal

int.random.test <- function(repsim, obs, nsim, 
												test = c("permutation", "bootstrap"), 
												alternative = c("auto","two.sided", 
																				"greater", "less"),
												adjust = "fdr") {

test <- match.arg(test)
alternative <- match.arg(alternative)


clase <- class(repsim)
if(clase == "matrix" | clase == "data.frame") {
	multi <- TRUE 
	N <- nrow(obs)
} else if(class(repsim) == "vector" | class(obs) == "integer" | class(obs) == "numeric") {
	multi <- FALSE
	N <- length(obs)
} else {
	stop("obs is not of class matrix, data.frame, vector, numeric or integer")
}


if(!multi) {
#single test
if(test == "permutation") {

if(!is.na(obs)) {
expected <- median(repsim, na.rm = TRUE)
if(alternative == "auto") {
	alter <- expected - obs
	if(alter > 0) {
		alternative <- "less"
	} else if (alter <= 0) {
		alternative <- "greater"
	}
}

if(alternative == "greater") {
	howmany <- sum(repsim >= obs, na.rm = TRUE)
	p.val <- (howmany + 1) / (nsim + 1)
	
} else if(alternative == "less") {
	howmany <- sum(repsim <= obs, na.rm = TRUE)
	p.val <- (howmany + 1) / (nsim + 1)
	
} else if(alternative == "two.sided") {
	howmany <- sum(abs(repsim) >= abs(obs), na.rm = TRUE)
	p.val <- (howmany + 1) / (nsim + 1)
}
confint <- quantile(repsim, probs = c(0.05, 0.95), na.rm = TRUE)
} else {
	expected <- NA
	alter <- NA
	p.val <- NA
	confint <- NA
}
result <- list(obs = obs, 
							 exp = expected, 
							 alter = alternative, 
							 p.val = p.val, 
							 CI = confint)

} else if(test == "bootstrap") {
	if(!is.na(obs)){
	result <- quantile(repsim, probs = c(0.05, 0.95), na.rm = TRUE)
	} else {
		CI <- NA
	}
	result <- list(obs = obs, CI = result)
}

#multiple tests
#variables in columns, repetitions in rows
} else {
	
	p.val <- numeric()
	expected <- numeric()
	altern <- rep("", ncol(repsim))
	alternative.i <- alternative
	
	if(test == "permutation") {
		
		for (i in 1:ncol(repsim)) {
			
			if(!is.na(obs[i])) {
			alternative <- alternative.i
			
			repet <- repsim[, i]
			expected[i] <- median(repet, na.rm = TRUE)
			
			if(is.na(expected[i])) {
				alternative <- "NA"
			}
			
			if(alternative == "auto") {
				alter <- expected[i] - obs[i]
				if(alter > 0) {
					alternative <- "less"
				} else if (alter <= 0) {
					alternative <- "greater"
				}
			}
			
			if(!is.na(alternative)) {
				altern[i] <- alternative
			} else {
				altern[i] <- NA
			}
			
			if(alternative == "greater") {
				howmany <- sum(repet >= obs[i], na.rm = TRUE)
				p.val[i] <- (howmany + 1) / (nsim + 1)
				
			} else if(alternative == "less") {
				howmany <- sum(repet <= obs[i], na.rm = TRUE)
				p.val[i] <- (howmany + 1) / (nsim + 1)
				
			} else if(alternative == "two.sided") {
				howmany <- sum(abs(repet) >= abs(obs[i]), na.rm = TRUE)
				p.val[i] <- (howmany + 1) / (nsim + 1)
				
			} else if(alternative == "NA") {
				p.val[i] <- NA
			}
			
			}  else {
			altern[i] <- NA
			p.val[i] <- NA
			expected[i] <- NA
			}
		}
		
		p.val <- p.adjust(p.val, method = adjust)
		
		if(!is.null(dim(repsim))) {
			qq <- apply(repsim, 2,	quantile, probs = c(0.05, 0.95),  
									na.rm = TRUE)
		} else {
			qq <- quantile (repsim, probs = c(0.05, 0.95),  na.rm = TRUE)
		}
		
		
		tab <- data.frame(outer(obs, c(1:4)))
		tab[, 1] <- round(obs, 4)
		tab[, 2] <- round(expected, 4)
		tab[, 3] <- altern
		tab[, 4] <- round(p.val, 4)
		tab[, 5:6] <- round(t(qq), 4)
		
		colnames(tab) <- c("obs", "exp", "alter", "p.val", "lwr", "uppr")
		
		result <- tab
		
	} else if(test == "bootstrap") {
		
		
		if(!is.null(dim(repsim))) {
			qq <- apply(repsim, 2,	quantile, probs = c(0.05, 0.95),  
									na.rm = TRUE)
		} else {
			qq <- quantile (repsim, probs = c(0.05, 0.95),  na.rm = TRUE)
		}
		
		tab <- data.frame(outer(obs, c(1:3)))
		tab[, 1] <- round(obs, 4)
		tab[, 2:3] <- round(t(qq), 4)
		
	
		colnames(tab) <- c("obs", "lwr", "uppr")
		result <- tab
		
	}
}


result

}

	

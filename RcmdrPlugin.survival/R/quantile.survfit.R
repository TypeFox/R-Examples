# last modified 2015-02-03 by J. Fox

# quantile.survfit <-
# 	function(x, quantiles=c(.25, .5, .75), ...){
# 	quants <- function(surv, lower, upper, t){
# 		if (length(surv) == 1) return(NA)
# 		warn <- options(warn=-1)
# 		on.exit(options(warn))
# 		select <- sapply(quantiles, function(q) q >= surv)
# 		posn <- apply(select, 2, function(x) min(which(x)))
# 		q <- t[posn]
# 		select <- sapply(quantiles, function(q) q >= lower)
# 		posn <- apply(select, 2, function(x) min(which(x)))
# 		low <- t[posn]
# 		select <- sapply(quantiles, function(q) q >= upper)
# 		posn <- apply(select, 2, function(x) min(which(x)))
# 		up <- t[posn]
# 		rbind(low, q, up)
# 	}
# 	summary <- summary(x)
# 	conf.level <- x$conf.int
# 	strata <- summary$strata
# 	if (is.null(strata)) {
# 		table <- quants(summary$surv, summary$lower, summary$upper, summary$time)
# 		dimnames(table) <- list(Estimate=c(paste("lower", conf.level, "CL"), "quantile", paste("upper", conf.level, "CL")), 
# 			"Survival Probability"=as.character(round(quantiles, 3)))
# 	} else {
# 		levels <- levels(strata)
# 		table <- array(0, c(3, length(quantiles), length(levels)))
# 		dimnames(table) <- list(Estimate=c(paste("lower", conf.level, "CL"), "quantile", paste("upper", conf.level, "CL")), 
# 			"Survival Probability"=as.character(round(quantiles, 3)), Stratum=levels)
# 		for (s in levels){
# 			select <- strata == s
# 			table[, , s] <- quants(summary$surv[select], summary$lower[select], 
# 				summary$upper[select], summary$time[select])
# 		}
# 	}
# 	table
# }


# weighted variance
var.wt <- function(x, w, na.rm = FALSE) {
    if (na.rm) {
        w <- w[i <- !is.na(x)]
        x <- x[i]
    }
    sum.w <- sum(w)
    return((sum(w*x^2) * sum.w - sum(w*x)^2) / (sum.w^2 - sum(w^2)))
}


weighted.t.test <- function(x, w, mu, conf.level = 0.95, alternative="two.sided", na.rm=TRUE) {
	
	if(!missing(conf.level) &&
       (length(conf.level) != 1 || !is.finite(conf.level) ||
        conf.level < 0 || conf.level > 1))
        stop("'conf.level' must be a single number between 0 and 1")
	
	if (na.rm) { 
		w <- w[i <- !is.na(x)] 
		x <- x[i] 
	}
	
	# to achieve consistent behavior in loops, return NA-structure in case of complete missings
	if (sum(is.na(x)) == length(x)) return(list(estimate=NA, se=NA, conf.int=NA, statistic=NA, df=NA, p.value=NA))
	
	# if only one value is present: this is the best estimate, no significance test provided
	if (sum(!is.na(x)) == 1) {
		warning("Warning weighted.t.test: only one value provided; this value is returned without test of significance!", call.=FALSE)
		return(list(estimate=x[which(!is.na(x))], se=NA, conf.int=NA, statistic=NA, df=NA, p.value=NA))
	}
	
	x.w <- weighted.mean(x,w, na.rm=na.rm)
	var.w <- var.wt(x,w, na.rm=na.rm)
	df <- length(x)-1
	t.value <- sqrt(length(x))*((x.w-mu)/sqrt(var.w))
	se <- sqrt(var.w)/sqrt(length(x))
	
	if (alternative == "less") {
		pval <- pt(t.value, df)
		cint <- c(-Inf, x.w + se*qt(conf.level, df) )
    }
    else if (alternative == "greater") {
		pval <- pt(t.value, df, lower.tail = FALSE)
		cint <- c(x.w - se * qt(conf.level, df), Inf)
    }
    else {
		pval <- 2 * pt(-abs(t.value), df)
		alpha <- 1 - conf.level
        cint <- x.w + se*qt(1 - alpha/2, df)*c(-1,1)
    }
	
	names(t.value) <- "t"
	return(list(estimate=x.w, se=se, conf.int=cint, statistic=t.value, df=df, p.value=pval))
}

## estimation functions should have
## INPUTS
##   ~ dat = data as a data.frame with columns 
##          - id, clust, per, trt, mean.y, y, at.risk.time
##   ~ incl.period.effect = indicator of whether to include a period effect
##   ~ outcome.type = one of "gaussian", "binomial", "poisson"
##   ~ alpha = the type 1 error rate
## OUTPUTS
##   ~ a vector with three elements, in order:
##       [1] a point estimate for the treatment effect
##       [2] lower bound of (1-alpha) confidence interval
##       [3] lower bound of (1-alpha) confidence interval


fixed.effect <- function(dat, incl.period.effect, outcome.type, alpha) {
	if(outcome.type=="poisson"){
		offsets <- log(dat[,"at.risk.time"])
	} else {
		offsets <- rep(0, nrow(dat))
	}

	if(incl.period.effect==0){
		fit <- glm(y ~ trt + clust,
			   data=dat,
			   family=outcome.type,
			   offset=offsets)
	} else {
		fit <- glm(y ~ trt + per + clust - 1,
			   data=dat,
			   family=outcome.type,
			   offset=offsets)
	}

	## get CI
	Z <- qnorm(1-alpha/2)*c(-1,1)
	est <- summary(fit)$coef["trt",]
	ci <- est["Estimate"] + Z*est["Std. Error"]
	return(c(est["Estimate"], ci))

}

fixed.effect.cluster.level <- function(dat, incl.period.effect, outcome.type, alpha) {
	cols <- c("y", "at.risk.time")

	clust.dat <- aggregate(dat[,cols], FUN=sum,
			       list(clust=dat[,"clust"],
				    per=dat[,"per"],
				    trt=dat[,"trt"]))
	if(outcome.type=="poisson"){
		offsets <- log(clust.dat[,"at.risk.time"])
	} else {
		offsets <- rep(0, nrow(clust.dat))
	}

	## set the formula
	if(incl.period.effect==0 & outcome.type!="binomial")
		form <- formula(y ~ trt + clust)
	if(incl.period.effect==0 & outcome.type=="binomial"){
		successes <- clust.dat[,"y"]
		failures <- clust.dat[,"at.risk.time"]-clust.dat[,"y"]
		form <- formula(cbind(successes, failures) ~ trt + clust)
	}
	if(incl.period.effect!=0 & outcome.type!="binomial")
		form <- formula(y ~ trt + per + clust - 1)
	if(incl.period.effect!=0 & outcome.type=="binomial") {
		successes <- clust.dat[,"y"]
		failures <- clust.dat[,"at.risk.time"]-clust.dat[,"y"]
		form <- formula(cbind(successes, failures) ~ trt + per + clust - 1)
	}


	fit <- glm(form, data=clust.dat, family=outcome.type, offset=offsets)

	## get CI
	Z <- qnorm(1-alpha/2)*c(-1,1)
	est <- summary(fit)$coef["trt",]
	ci <- est["Estimate"] + Z*est["Std. Error"]
	return(c(est["Estimate"], ci))

}


random.effect <- function(dat, incl.period.effect, outcome.type, alpha) {
	if(outcome.type=="poisson"){
		offsets <- log(dat[,"at.risk.time"])
	} else {
		offsets <- rep(0, nrow(dat))
	}

	if(incl.period.effect==0){
		fit <- glmer(y ~ trt + (1|clust),
			     data=dat,
			     family=outcome.type,
			     offset=offsets)
	} else {
		fit <- glmer(y ~ trt + per + (1|clust) - 1,
			     data=dat,
			     family=outcome.type,
			     offset=offsets)
	}

	n.clust <- length(unique(dat$clust))
	df <- n.clust - 2 ## based on k-2 in Donner & Klar p.118
	t <- qt(1 - alpha/2, df=df) * c(-1, 1)
	est <- coef(summary(fit))["trt", ]
	ci <- est["Estimate"] + t * est["Std. Error"]
	return(c(est["Estimate"], ci))
}

## based on Turner (1997) method C2
weighted.crossover.cluster.level <- function(dat, incl.period.effect, outcome.type, alpha) {
        if(outcome.type %in% c("poisson", "binomial"))
                stop("Currently, only Gaussian models are supported with this method.")
        if(unique(dat[,"per"])!=2)
                stop("This method only works for crosover studies with 2 periods.")

        ## takes the mean outcome by cluster-period
        clust.dat <- aggregate(dat[,"y"], FUN=mean,
                               list(clust=dat[,"clust"],
                                    per=as.numeric(dat[,"per"]),
                                    trt=dat[,"trt"]))
        ## makes a wide file with separate colums for treatment/control means
        clust.means <- reshape(clust.dat, idvar="clust", 
                               v.names=c("x", "per"), 
                               timevar="trt", direction="wide")
        clust.diffs <- clust.means[,"x.1"]-clust.means[,"x.0"]
        ws <- clust.means[,"per.1"]-clust.means[,"per.0"]
        
        clust.sizes <- table(dat$clust)
        wts <- clust.sizes/2 ## section 3.2.2 in Turner (1997)
        
        ## fit linear model
        fit <- lm(clust.diffs ~ ws, weights=wts)

        
}
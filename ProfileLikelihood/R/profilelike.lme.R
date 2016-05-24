profilelike.lme <-
function(formula, data, subject, random, correlation=NULL, profile.theta, method="ML", lo.theta, hi.theta, length=300, round=2, subset=NULL, weights=NULL, ...){
	if(!is.null(subset)){
		stop("Warning message: 'subset' should not be provided")
	}
	if(!is.null(weights)){
		stop("Warning message: 'weights' should not be provided")
	}
m <- model.frame(formula, data)
X <- model.matrix(formula, m)
y <- model.response(m)
theta.off <- data[,names(data)==profile.theta]
id <- data[,names(data)==subject]

	if(!is.numeric(theta.off)){
		stop("Warning message: 'profile.theta' must be a numeric variable")
	}
	if( ( length(theta.off)!= length(y) | length(theta.off)!= length(X[,1]) | length(y)!= length(X[,1]) ) ){
		cat("Warning message: remove missing data \n")
		}
		if( ( is.null(lo.theta) | is.null(hi.theta) )){
			cat("Warning message: provide lo.theta and hi.theta \n")
			fit <- lm(y ~ -1 + X + theta.off, na.action=na.fail)
			mle <- summary(fit)$coefficient["theta.off",1]
			se <- summary(fit)$coefficient["theta.off",2]
			lo.theta <- round(mle - 4*se, round)
			hi.theta <- round(mle + 4*se, round)
			}

theta <- seq(from =lo.theta, to=hi.theta, length=length)
log.lik <- rep(NA, length)

for(i in 1:length){
pi <- theta[i]
y.off <- y - pi*theta.off
fit <- lme(y.off ~ -1 + X, random = random, correlation=correlation, method="ML", na.action=na.fail)
log.lik[i] <- logLik(fit)
}

theta <- theta[is.na(log.lik)!=1]
log.lik <- log.lik[is.na(log.lik)!=1]
profile.lik <- exp(log.lik)

mm <- max(log.lik, na.rm=TRUE)
log.norm.lik <- log.lik - mm
profile.lik.norm <- exp(log.norm.lik)

return(list(theta=theta, profile.lik=profile.lik, profile.lik.norm=profile.lik.norm))
}


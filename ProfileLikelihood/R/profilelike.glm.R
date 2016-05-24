profilelike.glm <-
function(formula, data, profile.theta, family=gaussian, offset.glm=NULL, lo.theta=NULL, hi.theta=NULL, length=300, round=2, subset=NULL, weights=NULL, offset=NULL, ...){
	if(!is.null(subset)){
		stop("Warning message: 'subset' should not be provided")
	}
	if(!is.null(weights)){
		stop("Warning message: 'weights' should not be provided")
	}
	if(!is.null(offset)){
		stop("Warning message: do not use 'offset'; use 'offset.glm' instead of 'offset' ")
	}
m <- model.frame(formula, data)
X <- model.matrix(formula, m)
y <- model.response(m)
theta.off <- data[,names(data)==profile.theta]
	if(!is.numeric(theta.off)){
		stop("Warning message: 'profile.theta' must be a numeric variable")
	}
	if( ( length(theta.off)!= length(y) | length(theta.off)!= length(X[,1]) | length(y)!= length(X[,1]) ) ){
		cat("Warning message: remove missing data \n")
		}
		if( ( is.null(lo.theta) | is.null(hi.theta) )){
			cat("Warning message: provide lo.theta and hi.theta \n")
			fit <- glm(y ~ -1 + X + theta.off, family=family, na.action=na.fail)
			mle <- summary(fit)$coefficient["theta.off",1]
			se <- summary(fit)$coefficient["theta.off",2]
			lo.theta <- round(mle - 4*se, round)
			hi.theta <- round(mle + 4*se, round)
			}

theta <- seq(from =lo.theta, to=hi.theta, length=length)
log.lik <- rep(NA, length)

for(i in 1:length){
pi <- theta[i]
fit <- glm(y ~ -1 + X + offset(pi*theta.off), family=family, na.action=na.fail)
if(!is.null(offset.glm)){
glm.off <- data[,names(data)==offset.glm]
fit <- glm(y ~ -1 + X + offset(pi*theta.off) + offset(glm.off), family=family, na.action=na.fail)
}
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


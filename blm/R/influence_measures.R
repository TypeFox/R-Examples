influence.blm <- function(formula, data, weights=NULL,initial){
	
	fit <- do.call("glm", args=list(formula=formula, 
												data=data,
												start=initial,
												weights=weights,
												family=binomial(link=make.link("identity"))))

	lm.influence(fit)$coefficients
}

influence.lexpit <- function(formula, data, weights=NULL){
	
	fit <- do.call("glm", args=list(formula=formula, data=data, 
												weights=weights,family="binomial"))
	lm.influence(fit)$coefficients
}

vcov.influence.blm <- function(formula, data, initial){
	influence <- influence.blm(formula, data, weights=NULL, initial)
t(influence)%*%influence
}

vcov.influence.lexpit <- function(formula.linear, formula.expit, data, initial){

	influence.linear <- influence.blm(formula.linear, data, weights=NULL, initial)
	influence.expit <- influence.lexpit(formula.expit, data, initial)
	
	influence <- cbind(influence.linear[,-1],influence.expit)
	
t(influence)%*%influence
}


vcov.influence.blm.strata <- function(formula, data, weights, strata, initial){
	
	influence <- influence.blm(formula, data, weights, initial)
	means <- influence
	size <- table(strata)[strata]
	size <- ifelse(size==1,size, size/(size-1))

	for(i in 1:ncol(influence)){
		means[,i] <- (tapply(influence[,i], strata, mean)[strata])
		influence[,i] <- (influence[,i]-means[,i])*size
	}
	
t(influence)%*%influence
}


vcov.influence.lexpit.strata <- function(formula.linear, formula.expit, data, weights, strata, initial){

	influence.linear <- influence.blm(formula.linear, data, weights, initial)
	influence.expit <- influence.lexpit(formula.expit, data, weights)
	influence <- cbind(influence.linear[,-1],influence.expit)

	means <- influence
	size <- table(strata)[strata]
	size <- ifelse(size==1,size, size/(size-1))
	
	for(i in 1:ncol(influence)){
		means[,i] <- (tapply(influence[,i], strata, mean))[strata]
		influence[,i] <- (influence[,i]-means[,i])*size
	}
	
t(influence)%*%influence
}

vcov.blm.big <- function(formula, data, weights=NULL, initial){
	fit <- do.call("glm", args=list(formula=formula, data=data, start=initial,
												weights=weights,family=binomial(link=make.link("identity"))))
vcov(fit)
}

vcov.lexpit.big <- function(formula.linear, formula.expit, data, weights=NULL, initial){

	fit <- do.call("glm", args=list(formula=formula.linear, data=data, start = initial,
												weights=weights,family=binomial(link=make.link("identity"))))

   vcov.linear <- vcov(fit)
   
   	fit <- do.call("glm", args=list(formula=formula.expit, data=data, 
												weights=weights,family="binomial"))
	vcov.expit <- vcov(fit)

	p <- nrow(vcov.linear)	
	q <- nrow(vcov.expit)
	
	V <- matrix(0, p+q-1, p+q-1)
	
	V[(1:(p-1)),(1:(p-1))] <- vcov.linear[2:p,2:p] 
	V[p:(p+q-1),p:(p+q-1)] <- vcov.expit
	
V
}



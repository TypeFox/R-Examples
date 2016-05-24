speff <- function(formula, endpoint=c("quantitative", "dichotomous"), data, 
	postrandom=NULL, force.in=NULL, nvmax=9, method=c("exhaustive", "forward", 
	"backward"), optimal=c("cp", "bic", "rsq"), trt.id, conf.level=0.95, 
	missCtrl=NULL, missTreat=NULL, endCtrlPre=NULL, endTreatPre=NULL, 
	endCtrlPost=NULL, endTreatPost=NULL){
	require(leaps, quietly=TRUE)
	options(na.action=na.pass)
	if (missing(trt.id)) stop("Treatment indicator in 'trt.id' is missing.")
	endpoint <- match.arg(endpoint)
	method <- match.arg(method)
	optimal <- match.arg(optimal)
	mf <- match.call()
	mf$trt.id <- mf$endpoint <- mf$method <- mf$trt.id <- mf$conf.level <- 
	mf$optimal <- mf$postrandom <- mf$force.in <- mf$nvmax <- NULL
	mf[[1]] <- as.name("model.frame")
	mf <- eval(mf, parent.frame())
	X <- model.matrix(terms(formula), mf)
	Y <- Y0 <- model.response(mf)
	response <- as.character(formula[[2]])
	family <- ifelse(endpoint=="quantitative", "gaussian", "binomial")
	if (sum(c(is.null(endCtrlPre), is.null(endCtrlPost)))==1) endCtrlPost <- endCtrlPre
	if (sum(c(is.null(endTreatPre), is.null(endTreatPost)))==1) endTreatPost <- endTreatPre
	n <- table(ind <- ind0 <- as.factor(data[,trt.id]))
	levels(ind) <- 0:1
	naProb <- list(R <- ifelse(is.na(Y), 0, 1))
	naProb[[2]] <- R

	if (sum(is.na(Y))>0){
		data$R <- R
		miss <- list(missCtrl, missTreat)
		for (i in 1:2){
			if (is.null(miss[[i]])){ 
				optmod <- modSearch(as.formula(paste("R~",c(formula[[3]]))), X[ind==i-1,-1], 
				R[ind==i-1], endpoint, method, optimal, force.in, nvmax)
				naProb[[i]] <- predict(optmod$mod, data, type="response")
			} else {
				naProb[[i]] <- miss[[i]]
			}
		}
		Y0 <- na.omit(Y)
		naidx <- attributes(Y0)$na.action
		X <- X[-naidx,]; ind0 <- ind[-naidx]
	}
		
	optmod <- pVal <- pValBase <- rsq <- NULL
	epre <- list(endCtrlPre, endTreatPre)
	epost <- list(endCtrlPost, endTreatPost)
	for (i in 1:2){
		if (is.null(epre[[i]])){
			optmod[[i]] <- modSearch(formula, X[ind0==i-1,-1], Y0[ind0==i-1],
			endpoint, method, optimal, force.in, nvmax)
			namesBase <- setdiff(optmod[[i]]$names, postrandom)
			optmodBase <- glm(as.formula(paste(response,"~",paste(namesBase,collapse="+"))),
			family=family, data=na.omit(data[ind==i-1,]))
			rsq[i] <- optmod[[i]]$rsq
			pVal[[i]] <- predict(optmod[[i]]$mod, data, type="response")
			pValBase[[i]] <- predict(optmodBase, data, type="response")
		} else {
			pVal[[i]] <- epost[[i]]
			pValBase[[i]] <- epre[[i]]
		}
	}

	Y <- ifelse(is.na(Y), 0, Y)
	ind <- as.numeric(as.vector(ind))
	d <- n[2]/sum(n)
	if (endpoint=="quantitative"){
		mu0 <- sum(R*(1-ind)*Y/naProb[[1]] + (ind-d)*pValBase[[1]] - 
		(R-naProb[[1]])*(1-ind)*pVal[[1]]/naProb[[1]])/n[1]
		mu1 <- sum(R*ind*Y/naProb[[2]] - (ind-d)*pValBase[[2]] - 
		(R-naProb[[2]])*ind*pVal[[2]]/naProb[[2]])/n[2]
		beta <- mu1 - mu0
		vcov11 <- R*ind*(Y-mu1)/(d*naProb[[2]]) - (ind-d)*(pValBase[[2]]-mu1)/d -
		(R-naProb[[2]])*ind*(pVal[[2]]-mu1)/(d*naProb[[2]]) 
		vcov00 <- R*(1-ind)*(Y-mu0)/((1-d)*naProb[[1]]) - (ind-d)*(pValBase[[2]]-mu1)/(1-d) -
		(R-naProb[[1]])*(1-ind)*(pVal[[1]]-mu0)/((1-d)*naProb[[1]])
		vcov01 <- sum(vcov11*vcov00)/sum(n)^2
		vcov11 <- sum(vcov11^2)/sum(n)^2
		vcov00 <- sum(vcov00^2)/sum(n)^2
		vcov <- matrix(c(vcov00,vcov01,vcov01,vcov11),2,2)
		varbeta <- vcov11+vcov00-2*vcov01
	} else {
		expit <- function(x) exp(x)/(1+exp(x))
		logit <- function(p) log(p/(1-p))
		mu0 <- logit(sum((1-ind)*(Y-(R-naProb[[1]])*pVal[[1]])/naProb[[1]] +
		(ind-d)*pValBase[[1]])/n[1])
		mu1 <- logit(sum(ind*(Y-(R-naProb[[2]])*pVal[[2]])/naProb[[2]] -
		(ind-d)*pValBase[[2]])/n[2])
		beta <- mu1 - mu0
		m <- cbind((1-ind)*Y/naProb[[1]] + (ind-d)*pValBase[[1]] -
		(R-naProb[[1]])*(1-ind)*pVal[[1]]/naProb[[1]] - (1-d)*expit(mu0),
		ind*Y/naProb[[2]] - (ind-d)*pValBase[[2]] - (R-naProb[[2]])*ind*pVal[[2]]/
		naProb[[2]] - d*expit(mu1))
		B <- crossprod(m)/sum(n)
		A <- matrix(c((1-d)*exp(mu0)/(1+exp(mu0))^2, 0,
		rep(d*exp(mu1)/(1+exp(mu1))^2, 2)), 2, 2, byrow=TRUE)
		invA <- solve(A)
		v <- invA %*% B %*% t(invA)/sum(n)
		vcov <- matrix(c(v[1,1],sum(v[1,]),sum(v[1,]),sum(v)),2,2)
		varbeta <- v[2,2]
	}
		
	options(na.action=na.omit)
	unmod <- summary(glm(as.formula(paste(response,"~",trt.id)), family=family, data=data))
	unest <- unmod$coef[,1]
	coefs <- rbind(c(unest[1], sum(unest), unest[2]), c(mu0, mu1, beta))
	unv <- unmod$cov.scaled
	
	fits <- list(coef=coefs)
	class(fits) <- "speff"
	fits$cov <- list(naive=matrix(c(unv[1,1],sum(unv[1,]),sum(unv[1,]),sum(unv)),2,2), semi=vcov)
	fits$varbeta <- c(unv[2,2], varbeta)
	names(fits$varbeta) <- rownames(fits$coef) <- c("Naive", "Speff")
	colnames(fits$coef) <- c("Est Ctrl", "Est Treat", "Treat Effect")
	if (!is.null(optmod)) fits$formula <- list(control=formula(optmod[[1]]$mod), treatment=formula(optmod[[2]]$mod))
	fits$rsq <- rsq
	if (!is.null(rsq)) names(fits$rsq) <- c("Control", "Treatment")
	fits$endpoint <- endpoint
	fits$postrandom <- postrandom
	fits$predicted <- c(is.null(epre[[1]]), is.null(epre[[2]]))
	fits$conf.level <- conf.level
	fits$method <- method
	fits$n <- n
	fits
}

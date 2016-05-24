###############
#Find partial correlations of Xs
#(correlation of X with residuals from GLM/LMER fit without X)
###############

####
#Master function: chooses which subfunction (GLM/LMER) to use 
#arguments: 
#Y: response
#Z: treatment
#X: matrix of covariates
#resp.family: GLM family or "LMER" for LMER fit in response model
#trt.family: GLM family or "LMER" for LMER fit in treatment model
####

X.partials <- function(Y, Z, X, resp.family, trt.family) {
	if(class(resp.family) == "function"){
		fname = "X.partials.GLM"
	}else{
		fname <- "X.partials.GLM"
	}
	do.call(fname, list(Y, Z, X, resp.family, trt.family))
}

####
#Calculate partials for GLM
####

X.partials.GLM <- function(Y, Z, X, resp.family, trt.family) {
	nX <- dim(X)[2]
	if(is.null(nX))
		return(NULL)
	if(nX == 1) {
		XcorZ = cor(X, Z-mean(Z))
		fit.resp <- glm(Y~Z, resp.family)
		Yr <- Y-fit.resp$fitted.values
		XcorY <- cor(X, Yr)
	}else{
		XcorY <- XcorZ <- vector()
		for(i in 1:nX) {
			fit.resp <- glm(Y~X[,-i]+Z, resp.family)
			fit.trt <- glm(Z~X[,-i], trt.family)

			Yr <- Y-fit.resp$fitted.values
			Zr <- Z-fit.trt$fitted.values
		
			XcorY[i] <- cor(X[,i], Yr)
			XcorZ[i] <- cor(X[,i], Zr)
		}
	}
	return(cbind(XcorZ, XcorY))
}



# helpers.R

# compute interaction, squared, cubic, dummy variables, etc. for RSA
add.variables <- function(formula, df) {
	IV1 <- all.vars(formula)[2]
	IV2 <- all.vars(formula)[3]
	
	IV12 <- paste0(IV1, "2")
	IV22 <- paste0(IV2, "2")
	IV13 <- paste0(IV1, "3")
	IV23 <- paste0(IV2, "3")
	IV_IA <- paste0(IV1, "_", IV2)
	IV_IA2 <- paste0(IV1, "_", IV2, "2")
	IV_IA3 <- paste0(IV1, "2", "_", IV2)
	
	
	df[, IV12] <- df[, IV1]^2
	df[, IV22] <- df[, IV2]^2
	df[, IV_IA] <- df[, IV1]*df[, IV2]
	
	# three new variables for piecewise regression (test absolute difference score) - Edwards (2002) model
	df$W.JRE <- ifelse(df[, IV1] >= df[, IV2], 0, 1)
	df[, paste0("W.JRE_", IV1)] <- df$W.JRE*df[, IV1]
	df[, paste0("W.JRE_", IV2)] <- df$W.JRE*df[, IV2]
	
	# three new variables for piecewise regression (test absolute difference score) - new model Schoenbrodt 2012
	df$W <- ifelse(df[, IV1] >= df[, IV2], 1, -1)
	df$W[df[, IV1] == df[, IV2]] <- 0
	df[, paste0("W_", IV1)] <- df$W*df[, IV1]
	df[, paste0("W_", IV2)] <- df$W*df[, IV2]
	
	df$diff <- df[, IV2] - df[, IV1]
	df$SD <- df$diff^2
	df$absdiff <- abs(df$diff)
	
	# cubic terms
	df[, IV13] <- df[, IV1]^3
	df[, IV_IA2] <- df[, IV1]*df[, IV2]^2
	df[, IV_IA3] <- df[, IV1]^2*df[, IV2]
	df[, IV23] <- df[, IV2]^3
	
	return(df)
}


# helper function: takes a list of lavaan models (can include NULLs), and returns the usual anova object
anovaList <- function(modellist) {
	mods <- modellist[!sapply(modellist, function(x) is.null(x))]
	mods <- mods[!sapply(mods, function(x) !inspect(x, "converged"))]
	
	if (length(mods) == 0) {
		return(list(n.mods=0))
	}
	
    # put them in order (using df)
    DF <- sapply(mods, fitmeasures, "df")
    mods <- mods[order(DF, decreasing = FALSE)]
		
	pStr <- sapply(1:length(mods), function(x){ 
		if(x==1) {
			paste("mods[[",x,"]]",sep = "")
		} else {
			paste("force(mods[[",x,"]])",sep = "")
		}
	})
	#pStr2 <- paste0("lavTestLRT(", paste(pStr, collapse=", "), ", method='satorra.bentler.2010')")
	pStr2 <- paste0("lavTestLRT(", paste(pStr, collapse=", "), ", method='default')")
	
	a1 <- eval(parse(text = pStr2))
	
	if (length(mods) > 1) {
		rownames(a1) <- names(mods)
	}
	
	attr(a1, "n.mods") <- length(mods)
	return(list(ANOVA=a1, models=mods, n.mods=length(mods)))
}


## internal helper function: compare models
# mL = model list
# set = label that is attached to the results
cModels <- function(mL, set, free.max) {
	aL1 <- anovaList(mL)
	if (aL1$n.mods > 1) {
		N <- lavaan::nobs(aL1$models[[1]])
		a1 <- cbind(aL1$ANOVA[, c(1, 4:7)], plyr::ldply(aL1$models, function(X) {
			F <- fitmeasures(X)
			R <- inspect(X, "r2")
			names(R) <- "R2"
			n <- lavaan::nobs(X)
			k <- free.max - F["df"]		
			
			suppressWarnings({		
				R2.p <- ifelse(k==0,
					NA,
					pf(((n-k-1)*R)/(k*(1-R)), k, n-k-1, lower.tail=FALSE))
			})
			
			names(R2.p) <- "R2.p"
			
			# compute AICc
	      	AICc <- F["aic"] + 2*(k*(k+1))/(n-k-1)
			names(AICc) <- NULL
			
			return(c(AICc=AICc, F[c("cfi", "srmr")], R, R2.p))
		}))
		a1 <- a1[, !grepl(".id", colnames(a1))]
		a1$k <- free.max - a1$Df
		a1$R2.adj <- 1 - ((1-a1$R2))*((N-1)/(N-a1$k-1))
		a1$delta.R2 <- c(NA, a1$R2[1:(nrow(a1)-1)] - a1$R2[2:(nrow(a1))])			
		a1$model <- rownames(a1)
		a1$set <- set
		return(a1)
	}
}


# simple wrapper: formats a number in f.2 format
f2 <- function(x, digits=2, prepoint=0, skipZero=FALSE) {
	
	if (skipZero == TRUE) {zero <- "."} else {zero <- "0."}
	
	if (length(dim(x)) == 2) {
		apply(x, 2, function(x2) {gsub("0.", zero, sprintf(paste("%",prepoint,".",digits,"f",sep=""), x2) , fixed=TRUE)})
	} else {
		gsub("0.", zero, sprintf(paste("%",prepoint,".",digits,"f",sep=""), x) , fixed=TRUE)
	}
}

# converts p values in stars
p2star <- function(val) {
	
	res <- val
	
	for (i in 1:length(val)) {
		res[i] <- ""
		if (is.na(val[i])) next();
		if (val[i] <= 0.1) res[i] <- "\U2020"
		if (val[i] <= 0.05) res[i] <- "*"
		if (val[i] <= 0.01) res[i] <- "**"
		if (val[i] <= 0.001) res[i] <- "***"
	}
	
	return(res)
}

# nicely formats a p-value
p0 <- function(x) {
	if (is.na(x)) return("NA")
	if (x >= .001) return(paste0("p = ", f2(x, 3, skipZero=TRUE)))
	if (x <  .001) return("p <.001")	
}
p <- Vectorize(p0)

# returns number of maximum free parameters of a regression model
getFreeParameters <- function(model) {
	VARS <- nrow(inspect(model, "free")$beta)	# number of variables
	df.max <- (VARS*(VARS+1))/2		# maximum df
	df.pred <- ((VARS-1)*(VARS))/2 + 1 # df bound in the predictors (i.e., (co)variances of the predictors & variance of DV)
	free.max <- df.max - df.pred	# maximum of free parameters
	return(free.max)
}




# computes the coordinates of an arbitrary intersection of the surface,
# defined by a line on the X-Y plane (p0 = intercept, p1=slope)
getIntersect <- function(b0=0, x=0, y=0, x2=0, xy=0, y2=0, p0, p1, xlim=c(-2, 2), grid=21) {
	X <- seq(min(xlim), max(xlim), length.out=grid)
	Y <- p0 + p1*X
	n <- data.frame(X, Y)
	n2 <- add.variables(z~X+Y, n)
	n2$Z <- b0 + colSums(c(x, y, x2, y2, xy)*t(n2[, c(1:5)]))
	return(n2[, c("X", "Y", "Z")])
}


model <- function(x, model="full") x$models[[model]]

syntax <- function(x, model="full") cat(x$models[[model]]@Options$model)




# transforms p-values to colors
pRamp <- function(p, sig=.05, borderline=.10, bias=.8) {
	# calculate bias that the color transition is at the borderline value
	bias2 <- .33/(borderline/(1 - sig))
	cR1 <- colorRamp(c("red", "red", "orange"), bias=bias, space="Lab")
	cR2 <- colorRamp(c("orange", "green", "green"), bias=bias2, space="Lab")
	
	p2 <- rep("#FFFFFF", length(p))
	if (length(p[p < sig])>0) {
		p2[p < sig] <- rgb(cR1(p[p < sig]/sig), maxColorValue=255)
	}
	if (length(p[p >= sig])>0) {
		p2[p >= sig] <- rgb(cR2((p[p >= sig] - sig) / (1 - sig)), maxColorValue=255)
	}
	return(p2)
}


# helper function: find closest value in vector
f0 <- function (vec, target, unique = TRUE) {
    ret <- vec[sapply(target, function(x) which.min(abs(x - vec)))]
    if (unique) { ret <- unique(ret) }
    ret
}

# compute the predicted value from a single pair of predictors
predictRSA <- function(object, X, Y, model="full") {
	C <- coef(object$models[[model]])
	if (object$models[[model]]@Options$estimator != "DWLS") {
		b0 <- as.numeric(ifelse(is.na(C[paste0(object$DV, "~1")]), b0, C[paste0(object$DV, "~1")]))
		} else {
			# the threshold is the negative of the intercept ...
			b0 <- -as.numeric(ifelse(is.na(C[paste0(object$DV, "|t1")]), b0, C[paste0(object$DV, "|t1")]))
		}
	x <- as.numeric(ifelse(is.na(C["b1"]), 0, C["b1"]))
	y <- as.numeric(ifelse(is.na(C["b2"]), 0, C["b2"]))
	x2 <- as.numeric(ifelse(is.na(C["b3"]), 0, C["b3"]))
	y2 <- as.numeric(ifelse(is.na(C["b5"]), 0, C["b5"]))
	xy <- as.numeric(ifelse(is.na(C["b4"]), 0, C["b4"]))
	w <- as.numeric(ifelse(is.na(C["b6"]), 0, C["b6"]))
	wx <- as.numeric(ifelse(is.na(C["b7"]), 0, C["b7"]))
	wy <- as.numeric(ifelse(is.na(C["b8"]), 0, C["b8"]))
	
	# cubic parameters
	x3 <- as.numeric(ifelse(is.na(C["b9"]), 0, C["b9"]))
	xy2 <- as.numeric(ifelse(is.na(C["b10"]), 0, C["b10"]))
	x2y <- as.numeric(ifelse(is.na(C["b11"]), 0, C["b11"]))
	y3 <- as.numeric(ifelse(is.na(C["b12"]), 0, C["b12"]))
	
	
	C <- c(x, y, x2, y2, xy, w, wx, wy,x3, xy2, x2y, y3)
	
	# compute predicted value
	Z <- b0 + colSums(C*t(cbind(X, Y, X^2, Y^2, X*Y, 0, 0, 0, X^3, X*Y^2, X^2*Y, Y^3)))
	return(Z)
}


# fills up the long edges of a polygon with intermediate points
# If an edge is longer than minDist, new points re inserted.
# @param x Vector of x values
# @param y Vector of y values
interpolatePolygon <- function(x, y, minDist, plot=FALSE) {
	minDist <- minDist^2	# compare with squared x^2 + y^2 (faster)
	interp <- data.frame()
	pol <- data.frame(x, y)
	colnames(pol) <- c("x", "y")
	for (i in 1:(nrow(pol)-1)) {
		# get distance
		D <- (pol[i, 1] - pol[i+1, 1])^2 + (pol[i, 2] - pol[i+1, 2])^2
		if (D > minDist) {
			N <- ceiling(sqrt(D)/sqrt(minDist)) # number of interpolations
			APPROX <- data.frame(
				x = seq(pol[i, 1], pol[i+1, 1], length.out=N),
				y = seq(pol[i, 2], pol[i+1, 2], length.out=N)
			)
			interp <- rbind(interp, APPROX)
		} else if (D>0 & D <= minDist){
			interp <- rbind(interp, pol[i, ])
			if (i==1) colnames(interp) <- c("x", "y")
		}
	}
	interp <- rbind(interp, pol[nrow(pol), ])
	if (plot==TRUE) {
		plot(pol, col="red")
		points(interp[, 1], interp[, 2], col="green", pch=20)
		lines(interp[, 1], interp[, 2], col="darkgreen")
		text(pol, label=1:nrow(pol))
	}
	return(interp)
}


# agd.r
#
# Tools for course Analysis of Growth Data
### SvB 20dec2014

#'Convert standard deviation scores (SDS) to measurements 
#'
#'Converts standard deviation score (SDS) into measurements using 
#'an age- and sex-conditional external reference.
#'
#'Functions \code{z2y()} and \code{y2z()} are the inverse of each other.
#'
#'The argument \code{dist} determines the statistical distribution.  The
#'possibilities are as follows: \describe{ \item{list("\"NO\"")}{\code{ref}
#'should contain columns \code{mean} and \code{sd}, containing the mean and the
#'standard deviation in the external reference population.}
#'\item{list("\"LMS\"")}{\code{ref} should contain columns \code{L}, \code{S}
#'and \code{M} containing the LMS parameters.}
#'\item{list("\"BCCG\"")}{\code{ref} should contain columns \code{mu},
#'\code{sigma} and \code{nu} containing the Box-Cox Cole-Green parameters.}
#'\item{list("\"BCPE\"")}{\code{ref} should contain columns \code{mu},
#'\code{sigma}, \code{nu} and \code{tau} containing the Box-Cox Power
#'Exponential parameters.} \item{list("\"BCT\"")}{\code{ref} should contain
#'columns \code{mu}, \code{sigma}, \code{nu} and \code{tau} containing the
#'Box-Cox T distribution parameters.} }
#'
#'@aliases z2y
#'@param z A numerical vector containing standard deviation scores that are to
#'be converted. The length \code{length(z)} determines the size of the output
#'vector.
#'@param x A vector containing the values of the numerical covariate (typically
#'decimal age or height) at which conversion is desired.  Values are replicated
#'to match \code{length(y)}.
#'@param sex A character vector indicating whether the male (\code{"M"}) of
#'female (\code{"F"})reference should be used.  Values are replicated to match
#'\code{length(y)}.
#'@param sub A character vector indicating the level of the \code{sub} field of
#'the reference standard defined in \code{ref}
#'@param ref A data frame containing a factor \code{sex}, a numerical variable
#'\code{age} containing the tabulated decimal point ages, and two or more
#'numerical variables with reference values. See details.
#'@param dist A string identifying the type of distribution. Values values are:
#'\code{"NO"}, \code{"BCCG"}, \code{"LMS"}, \code{"BCPE"} and \code{"BCT"}.
#'The default is \code{"LMS"}.
#'@param dec A scalar value indicating the number of decimals used to round the
#'value.
#'@param sex.fallback The level of the \code{sex} field used when no match is
#'found.  The default is \code{"M"} for males. Specify \code{sex.fallback="NA"}
#'if unmatched entries should receive a \code{NA} value.
#'@param sub.fallback The level of the \code{sub} field used when no match is
#'found.  The default is \code{"N"} for normal. Specify
#'\code{sub.fallback="NA"} if unmatched entries should receive a \code{NA}
#'value.
#'@return For \code{y2z()}: A vector with \code{length(y)} elements containing
#'the standard deviation score.  For \code{z2y()}: A vector with
#'\code{length(z)} elements containing quantiles.
#'@author Stef van Buuren, 2010
#'@seealso \code{\link{y2z}}
#'@keywords distribution
#'@examples
#'
#'
#'boys <- boys7482
#'
#'# quantile at SD=0 of age 2 years, 
#'# height Dutch boys
#'z2y(z=0, x=2)
#'
#'# same for Dutch girls
#'z2y(z=0, x=2, sex="F")
#'
#'# quantile at SD=c(-1,0,1) of age 2 years, BMI Dutch boys
#'z2y(z=c(-1,0,+1), x=2, ref=nl4.bmi)
#'
#'# 0SD line (P50) in kg of weight for age in 5-10 year, Dutch boys
#'z2y(z=rep(0,6), x=5:10, ref=nl4.wgt)
#'
#'# 95th percentile (P95), age 10 years, wfa, Dutch boys
#'z2y(z=qnorm(0.95), x=10, ref=nl4.wgt)
#'
#'# table of P3, P10, P50, P90, P97 of weight for 5-10 year old dutch boys
#'# age per year
#'age <- 5:10
#'p <- c(0.03,0.1,0.5,0.9,0.97)
#'z <- rep(qnorm(p), length(age))
#'x <- rep(age, each=length(p))
#'w <- matrix(z2y(z, x=x, sex="M", ref=nl4.wgt), ncol=length(p),
#'  byrow=TRUE)
#'dimnames(w) <- list(age, p)
#'round(w,1)
#'
#'# standard set of Z-scores of weight for all tabulated ages, boys & girls
#'# and three etnicities
#'sds <- c(-2.5, -2, -1, 0, 1, 2, 2.5)
#'age <- nl4.wgt$x
#'z <- rep(sds, times=length(age))
#'x <- rep(age, each=length(sds))
#'sex <- rep(c("M","F"), each=length(z)/2)
#'w <- z2y(z=z, x=x, sex=sex, ref=nl4.wgt)
#'w <- matrix(w, ncol=length(sds), byrow=TRUE)
#'dimnames(w) <- list(age, sds)
#'data.frame(sub=nl4.wgt$sub,sex=nl4.wgt$sex,round(w,2), row.names=NULL)
#'
#'# P85 of BMI in 5-8 year old Dutch boys and girls
#'e <- expand.grid(age=5:8, sex=c("M","F"))
#'w <- z2y(z=rep(qnorm(0.85),nrow(e)), x=e$age, sex=e$sex, ref=nl4.bmi)
#'w <- matrix(w, nrow=2, byrow=TRUE)
#'dimnames(w) <- list(c("boys","girls"),5:8)
#'w
#'
#'# data transformation of height z-scores to cm-scale
#'z <- c(-1.83, 0.09, 2.33, 0.81, -1.20)
#'x <- c(8.33,  0.23, 19.2, 24.3, 10)
#'sex <- c("M", "M", "F", "M", "F")
#'round(z2y(z=z, x=x, sex=sex, ref=nl4.hgt), 1)
#'
#'# interpolate published height standard 
#'# to daily values, days 0-31, boys
#'# on centiles -2SD, 0SD and +2SD 
#'days <- 0:31
#'sds  <- c(-2, 0, +2)
#'z    <- rep(sds, length(days))
#'x    <- rep(round(days/365.25,4), each=length(sds))
#'w    <- z2y(z, x, sex="M", ref=nl4.hgt)
#'w    <- matrix(w, ncol=length(sds), byrow=TRUE)
#'dimnames(w) <- list(days, sds)
#'w
#'
#'

z2y <- function(z   = c(-2, 0, 2), 
				x   = 1, 
				sex = "M", 
				sub = "N",
				ref = get("nl4.hgt"),
				dist = "LMS",
				dec  = 3,
				sex.fallback = "M",
				sub.fallback = "N")
{
	
	z2y.grp <- function(z, x, ref, dist = "LMS") {
		
		if (dist=="NO") {
			check.names(df=ref, needed=c("x","mean","sd"))
			mean <- approx(x=ref[,"x"], y=ref[,"mean"], xout=x)$y
			sd   <- approx(x=ref[,"x"], y=ref[,"sd"],   xout=x)$y 
			return(mean + z*sd)   
		}
		if (dist=="LMS") {
			check.names(df=ref, needed=c("x","L","M","S"))
			L <- approx(x=ref[,"x"], y=ref[,"L"], xout=x)$y
			M <- approx(x=ref[,"x"], y=ref[,"M"], xout=x)$y
			S <- approx(x=ref[,"x"], y=ref[,"S"], xout=x)$y
			return(ifelse(L>0.01|L<(-0.01),M*(1+L*S*z)^(1/L),M*exp(S*z)))
		}
		if (dist=="BCCG") {
			check.names(df=ref, needed=c("x","nu","mu","sigma"))
			nu    <- approx(x=ref[,"x"], y=ref[,"nu"], xout=x)$y
			mu    <- approx(x=ref[,"x"], y=ref[,"mu"], xout=x)$y
			sigma <- approx(x=ref[,"x"], y=ref[,"sigma"], xout=x)$y
			return(qBCCG(pnorm(z), mu=mu, sigma=sigma, nu=nu))
		}
		if (dist=="BCPE") {
			check.names(df=ref, needed=c("x","nu","mu","sigma","tau"))
			mu    <- approx(x=ref[,"x"], y=ref[,"mu"], xout=x)$y
			sigma <- approx(x=ref[,"x"], y=ref[,"sigma"], xout=x)$y
			nu    <- approx(x=ref[,"x"], y=ref[,"nu"], xout=x)$y
			tau   <- approx(x=ref[,"x"], y=ref[,"tau"], xout=x)$y
			return(qBCPE(pnorm(z), mu=mu, sigma=sigma, nu=nu, tau=tau))
		}
		if (dist=="BCT") {
			check.names(df=ref, needed=c("x","nu","mu","sigma","tau"))
			mu    <- approx(x=ref[,"x"], y=ref[,"mu"], xout=x)$y
			sigma <- approx(x=ref[,"x"], y=ref[,"sigma"], xout=x)$y
			nu    <- approx(x=ref[,"x"], y=ref[,"nu"], xout=x)$y
			tau   <- approx(x=ref[,"x"], y=ref[,"tau"], xout=x)$y
			return(qBCT(pnorm(z), mu=mu, sigma=sigma, nu=nu, tau=tau))
		}
		
		stop(paste("Reference type", dist, "not implemented."))
	}
	
	if (!is.data.frame(ref)) stop("'ref' should be a data frame.")
	n <- length(z)
	if (n < 1)         stop("'z' must have 1 or more values")
	if(!is.vector(z))  stop("'z' must be a numeric vector")
	if(!is.numeric(z)) stop("'z' must be a numeric vector")
	
	x   <- rep(x,   length.out=length(z))
	sex <- rep(sex, length.out=length(z))
	sub <- rep(sub, length.out=length(z))
	dist <- match.arg(dist, choices=c("NO","LMS","BCCG","BCPE","BCT"))
	
	# available levels in ref: sex, sub
	lev.sex <- levels(ref$sex[, drop=TRUE])
	lev.sub <- levels(ref$sub[, drop=TRUE])
	
	# replace nomatching levels
	idx <- is.na(match(sub, lev.sub))
	if (any(idx)) {
		sub[idx] <- sub.fallback
		warning("Entries (n=",sum(idx),") replaced by '",sub.fallback,"'",sep="")
	}
	idx <- is.na(match(sex, lev.sex))
	if (any(idx)) {
		sex[idx] <- sex.fallback
		warning("Entries (n=",sum(idx),") replaced by '",sex.fallback,"'",sep="")
	}
	
	refs <- with(ref,split(ref, f=list(sub, sex)))
	xs <- split(x,list(sub, sex))
	zs <- split(z,list(sub, sex))
	ys <- vector("list",length(zs))
	names(ys) <- names(zs)
	
	for(i in 1:length(zs)) {
		name <- names(zs)[i]
		if(is.null(refs[[name]])) ys[[name]] <- rep(NA,length=length(zs[[name]]))
		else ys[[name]] <- z2y.grp(z=zs[[name]], x=xs[[name]], ref=refs[[name]], 
								   dist=dist)
	}
	
	y <- unsplit(ys,f=list(sub,sex))
	names(y) <- names(z)
	return(round(y,dec))
}





#'Converts measurements to standard deviation scores (SDS)
#'
#'Converts measurements into age- and sex-conditional standard deviation score
#'(SDS) using an external reference.
#'
#'Functions \code{z2y()} and \code{y2z()} are the inverse of each other.
#'
#'The argument \code{dist} determines the statistical distribution.  The
#'possibilities are as follows: \describe{ \item{list("\"NO\"")}{\code{ref}
#'should contain columns \code{mean} and \code{sd}, containing the mean and the
#'standard deviation in the external reference population.}
#'\item{list("\"LMS\"")}{\code{ref} should contain columns \code{L}, \code{S}
#'and \code{M} containing the LMS parameters.}
#'\item{list("\"BCCG\"")}{\code{ref} should contain columns \code{mu},
#'\code{sigma} and \code{nu} containing the Box-Cox Cole-Green parameters.}
#'\item{list("\"BCPE\"")}{\code{ref} should contain columns \code{mu},
#'\code{sigma}, \code{nu} and \code{tau} containing the Box-Cox Power
#'Exponential parameters.} \item{list("\"BCT\"")}{\code{ref} should contain
#'columns \code{mu}, \code{sigma}, \code{nu} and \code{tau} containing the
#'Box-Cox T distribution parameters.} }
#'
#'@aliases y2z
#'@param y A numerical vector containing the outcome measurements.  The length
#'\code{length(y)} determines the size of the output vector.
#'@param x A vector containing the values of the numerical covariate (typically
#'decimal age or height) at which conversion is desired.  Values are replicated
#'to match \code{length(y)}.
#'@param sex A character vector indicating whether the male (\code{"M"}) of
#'female (\code{"F"})reference should be used.  Values are replicated to match
#'\code{length(y)}.
#'@param sub A character vector indicating the level of the \code{sub} field of
#'the reference standard defined in \code{ref}
#'@param ref A data frame containing a factor \code{sex}, a numerical variable
#'\code{age} containing the tabulated decimal point ages, and two or more
#'numerical variables with reference values. See details.
#'@param dist A string identifying the type of distribution. Values values are:
#'\code{"NO"}, \code{"BCCG"}, \code{"LMS"}, \code{"BCPE"} and \code{"BCT"}.
#'The default is \code{"LMS"}.
#'@param dec A scalar value indicating the number of decimals used to round the
#'value.
#'@param sex.fallback The level of the \code{sex} field used when no match is
#'found.  The default is \code{"M"} for males. Specify \code{sex.fallback="NA"}
#'if unmatched entries should receive a \code{NA} value.
#'@param sub.fallback The level of the \code{sub} field used when no match is
#'found.  The default is \code{"N"} for normal. Specify
#'\code{sub.fallback="NA"} if unmatched entries should receive a \code{NA}
#'value.
#'@param tail.adjust Logical. If \code{TRUE} then the WHO method for 
#'tail adjustment is applied. The default is \code{FALSE}.
#'@return For \code{y2z()}: A vector with \code{length(y)} elements containing
#'the standard deviation score.  For \code{z2y()}: A vector with
#'\code{length(z)} elements containing quantiles.
#'@author Stef van Buuren, 2010
#'@seealso \code{\link{z2y}}
#'@keywords distribution
#'@examples
#'
#'
#'boys <- boys7482
#'
#'# SDS of height 115 cm at age 5 years, 
#'# relative to Dutch boys reference
#'y2z(y=115, x=5)
#'
#'# same relative to Dutch girls
#'y2z(y=115, x=5, sex="F")
#'
#'# SDS of IOTF BMI cut-off value for overweight (boys 2-18) 
#'# relative to Dutch boys reference
#'cutoff <- c(
#'18.41, 18.15, 17.89, 17.72, 17.55, 17.49, 17.42, 17.49, 17.55, 17.74,
#'17.92, 18.18, 18.44, 18.77, 19.10, 19.47, 19.84, 20.20, 20.55, 20.89, 
#'21.22, 21.57, 21.91, 22.27, 22.62, 22.96, 23.29, 23.60, 23.90, 24.18, 
#'24.46, 24.73, 25.00)
#'age <- seq(2, 18, by=0.5)
#'(z <- y2z(y=cutoff, x=age, sex="M", ref=nl4.bmi))
#'
#'# apply inverse transformation to check calculations
#'round(z2y(z, age, ref=nl4.bmi), 2)
#'cutoff
#'
#'# calculate percentiles of weight 12 kg at 2 years (boys, girls)
#'100*round(pnorm(y2z(y=c(12,12), x=2, sex=c("M","F"), ref=nl4.wgt)),2)
#'
#'# # percentage of children lighter than 15kg at ages 2-5
#'e <- expand.grid(age=2:5, sex=c("M","F"))
#'z <- y2z(y=rep(15,nrow(e)), x=e$age, sex=e$sex, ref=nl4.wgt)
#'w <- matrix(100*round(pnorm(z),2), nrow=2, byrow=TRUE)
#'dimnames(w) <- list(c("boys","girls"),2:5)
#'w
#'
#'# analysis in Z scale
#'hgt.z <- y2z(y=boys$hgt, x=boys$age, sex="M", ref=nl4.hgt)
#'wgt.z <- y2z(y=boys$wgt, x=boys$age, sex="M", ref=nl4.wgt)
#'plot(hgt.z, wgt.z, col="blue")
#'
#'
#'# z2y
#'
#'# quantile at SD=0 of age 2 years, 
#'# height Dutch boys
#'z2y(z=0, x=2)
#'
#'# same for Dutch girls
#'z2y(z=0, x=2, sex="F")
#'
#'# quantile at SD=c(-1,0,1) of age 2 years, BMI Dutch boys
#'z2y(z=c(-1,0,+1), x=2, ref=nl4.bmi)
#'
#'# 0SD line (P50) in kg of weight for age in 5-10 year, Dutch boys
#'z2y(z=rep(0,6), x=5:10, ref=nl4.wgt)
#'
#'# 95th percentile (P95), age 10 years, wfa, Dutch boys
#'z2y(z=qnorm(0.95), x=10, ref=nl4.wgt)
#'
#'# table of P3, P10, P50, P90, P97 of weight for 5-10 year old dutch boys
#'# age per year
#'age <- 5:10
#'p <- c(0.03,0.1,0.5,0.9,0.97)
#'z <- rep(qnorm(p), length(age))
#'x <- rep(age, each=length(p))
#'w <- matrix(z2y(z, x=x, sex="M", ref=nl4.wgt), ncol=length(p),
#'  byrow=TRUE)
#'dimnames(w) <- list(age, p)
#'round(w,1)
#'
#'# standard set of Z-scores of weight for all tabulated ages, boys & girls
#'# and three etnicities
#'sds <- c(-2.5, -2, -1, 0, 1, 2, 2.5)
#'age <- nl4.wgt$x
#'z <- rep(sds, times=length(age))
#'x <- rep(age, each=length(sds))
#'sex <- rep(c("M","F"), each=length(z)/2)
#'w <- z2y(z=z, x=x, sex=sex, ref=nl4.wgt)
#'w <- matrix(w, ncol=length(sds), byrow=TRUE)
#'dimnames(w) <- list(age, sds)
#'data.frame(sub=nl4.wgt$sub,sex=nl4.wgt$sex,round(w,2), row.names=NULL)
#'
#'# P85 of BMI in 5-8 year old Dutch boys and girls
#'e <- expand.grid(age=5:8, sex=c("M","F"))
#'w <- z2y(z=rep(qnorm(0.85),nrow(e)), x=e$age, sex=e$sex, ref=nl4.bmi)
#'w <- matrix(w, nrow=2, byrow=TRUE)
#'dimnames(w) <- list(c("boys","girls"),5:8)
#'w
#'
#'# data transformation of height z-scores to cm-scale
#'z <- c(-1.83, 0.09, 2.33, 0.81, -1.20)
#'x <- c(8.33,  0.23, 19.2, 24.3, 10)
#'sex <- c("M", "M", "F", "M", "F")
#'round(z2y(z=z, x=x, sex=sex, ref=nl4.hgt), 1)
#'
#'# interpolate published height standard 
#'# to daily values, days 0-31, boys
#'# on centiles -2SD, 0SD and +2SD 
#'days <- 0:31
#'sds  <- c(-2, 0, +2)
#'z    <- rep(sds, length(days))
#'x    <- rep(round(days/365.25,4), each=length(sds))
#'w    <- z2y(z, x, sex="M", ref=nl4.hgt)
#'w    <- matrix(w, ncol=length(sds), byrow=TRUE)
#'dimnames(w) <- list(days, sds)
#'w
#'
#'
y2z <- function(y   = c(75, 80, 85), 
				x   = 1, 
				sex = "M", 
				sub = "N",
				ref = get("nl4.hgt"),
				dist = "LMS",
				dec  = 3,
				sex.fallback = "M",
				sub.fallback = "N",
				tail.adjust = FALSE)
{
	
	y2z.grp <- function(y, x, ref, dist = "LMS", dec = 3, tail.adjust = FALSE){
		
		if (dist=="NO") {
			check.names(df=ref, needed=c("x","mean","sd"))
			mean <- approx(x=ref[,"x"], y=ref[,"mean"], xout=x)$y
			sd   <- approx(x=ref[,"x"], y=ref[,"sd"],   xout=x)$y 
			return((y-mean)/sd)   
		}
		if (dist=="LMS") {
			check.names(df=ref, needed=c("x","L","M","S"))
			L <- approx(x=ref[,"x"], y=ref[,"L"], xout=x)$y
			M <- approx(x=ref[,"x"], y=ref[,"M"], xout=x)$y
			S <- approx(x=ref[,"x"], y=ref[,"S"], xout=x)$y
			z <- ifelse(L>0.01 | L<(-0.01), (((y/M)^L)-1)/(L*S), log(y/M)/S)
			if (tail.adjust) z <- adjust.tail(y, z, L, M, S) 
			return(z)
		}
		if (dist=="BCCG") {
			check.names(df=ref, needed=c("x","nu","mu","sigma"))
			nu    <- approx(x=ref[,"x"], y=ref[,"nu"], xout=x)$y
			mu    <- approx(x=ref[,"x"], y=ref[,"mu"], xout=x)$y
			sigma <- approx(x=ref[,"x"], y=ref[,"sigma"], xout=x)$y
			return(qnorm(pBCCG(y, mu=mu, sigma=sigma, nu=nu)))
		}
		if (dist=="BCPE") {
			check.names(df=ref, needed=c("x","nu","mu","sigma","tau"))
			mu    <- approx(x=ref[,"x"], y=ref[,"mu"], xout=x)$y
			sigma <- approx(x=ref[,"x"], y=ref[,"sigma"], xout=x)$y
			nu    <- approx(x=ref[,"x"], y=ref[,"nu"], xout=x)$y
			tau   <- approx(x=ref[,"x"], y=ref[,"tau"], xout=x)$y
			return(qnorm(pBCPE(y, mu=mu, sigma=sigma, nu=nu, tau=tau)))
		}
		if (dist=="BCT") {
			check.names(df=ref, needed=c("x","nu","mu","sigma","tau"))
			mu    <- approx(x=ref[,"x"], y=ref[,"mu"], xout=x)$y
			sigma <- approx(x=ref[,"x"], y=ref[,"sigma"], xout=x)$y
			nu    <- approx(x=ref[,"x"], y=ref[,"nu"], xout=x)$y
			tau   <- approx(x=ref[,"x"], y=ref[,"tau"], xout=x)$y
			return(qnorm(pBCT(y, mu=mu, sigma=sigma, nu=nu, tau=tau)))
		}
		if (dist=="BCCG") {
			check.names(df=ref, needed=c("x","mu","sigma","nu"))
			mu    <- approx(x=ref[,"x"], y=ref[,"mu"], xout=x)$y
			sigma <- approx(x=ref[,"x"], y=ref[,"sigma"], xout=x)$y
			nu    <- approx(x=ref[,"x"], y=ref[,"nu"], xout=x)$y
			return()
			#return(ifelse(L>0.01|L<(-0.01),(((y/M)^L)-1)/(L*S),log(y/M)/S))
		}
		stop(paste("Reference type", dist, "not implemented."))
	}	
	
	if (!is.data.frame(ref)) stop("'ref' should be a data frame.")
	n <- length(y)
	if (n < 1)         stop("'y' must have 1 or more values")
	if(!is.vector(y))  stop("'y' must be a numeric vector")
	if(!is.numeric(y)) stop("'y' must be a numeric vector")
	
	x   <- rep(x,   length.out=length(y))
	sex <- rep(sex, length.out=length(y))
	sub <- rep(sub, length.out=length(y))
	dist <- match.arg(dist, choices=c("NO","LMS","BCCG","BCPE","BCT"))
	
	# available levels in ref: sex, sub
	lev.sex <- levels(ref$sex[, drop=TRUE])
	lev.sub <- levels(ref$sub[, drop=TRUE])
	
	# replace nomatching levels
	idx <- is.na(match(sub, lev.sub))
	if (any(idx)) {
		sub[idx] <- sub.fallback
		warning("Entries (n=",sum(idx),") replaced by '",sub.fallback,"'",sep="")
	}
	idx <- is.na(match(sex, lev.sex))
	if (any(idx)) {
		sex[idx] <- sex.fallback
		warning("Entries (n=",sum(idx),") replaced by '",sex.fallback,"'",sep="")
	}
	
	refs <- with(ref,split(ref, f=list(sub, sex)), drop = TRUE)
	xs <- split(x,list(sub, sex), drop = TRUE)
	ys <- split(y,list(sub, sex), drop = TRUE)
	zs <- vector("list",length(ys))
	names(zs) <- names(ys)
	
	for(i in 1:length(ys)) {
		name <- names(ys)[i]
		if(is.null(refs[[name]])) ys[[name]] <- rep(NA,length=length(ys[[name]]))
		else zs[[name]] <- y2z.grp(y=ys[[name]], x=xs[[name]], ref=refs[[name]], 
								   dist=dist, tail.adjust = tail.adjust)
	}
	
	z <- unsplit(zs,f=list(sub,sex))		
	names(z) <- names(y)
	return(round(z, dec))
}

check.names <- function(df, needed){
	if (missing(df)) stop("required argument 'df' not found")
	if (missing(needed)) stop("required argument 'needed' not found")
	
	notfound <- is.na(match(needed, names(df)))
	if (any(notfound)) stop("Not found: ",paste(needed[notfound],collapse=", "))
}

adjust.tail <- function(y, z, L, M, S){
	idx <- !is.na(z) & z > 3
	if (any(idx)) {
		sd3 <- ifelse(L>0.01|L<(-0.01), M*(1+L*S*3)^(1/L), M*exp(S*3))
		sd2 <- ifelse(L>0.01|L<(-0.01), M*(1+L*S*2)^(1/L), M*exp(S*2))
		z[idx] <- (3 + (y - sd3)/(sd3 - sd2))[idx]
	}
	idx <- !is.na(z) & z < (-3)
	if (any(idx)) {
		sd3 <- ifelse(L>0.01|L<(-0.01), M*(1+L*S*(-3))^(1/L), M*exp(S*(-3)))
		sd2 <- ifelse(L>0.01|L<(-0.01), M*(1+L*S*(-2))^(1/L), M*exp(S*(-2)))
		z[idx] <- (-3 + (y - sd3)/(sd2 - sd3))[idx]
	}
	return(z)
}

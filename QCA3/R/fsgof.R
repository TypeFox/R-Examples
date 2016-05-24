## Programme written by Ronggui HUANG (2010)
## Goodness-of-Fit Tests and Descriptive Measures in Fuzzy-Set Analysis
## Eliason S. & Stryker R. 2009. Sociological Methods & Research 38:102-146. 

zTransform <- function(x,damp=.01){
      x2 <- x
      x2[x < damp] <- damp
      x2[x > (1-damp)] <- 1-damp
	zx <- qnorm(x2)
	zx
}

Dnull <- function(y, damp=.01){
	yZ <- zTransform(y,damp=damp)
	df <- length(y) - 1
	ssd <- sum((yZ -mean(yZ))^2)
	msd <- ssd/df
	ans <- list(ssd=ssd,df=df,msd=msd)
	ans
}

fsgof.nec <- function(y, x, damp=.01, error=.05) {
	dn <- sprintf("%s and %s", deparse(substitute(y)),deparse(substitute(x)))
	yZ <- zTransform(y,damp=damp)
	xZ <- zTransform(x,damp=damp)
	d <- as.numeric(y > x)
	ssd <- sum(d*(yZ - xZ)^2)
	df1 <- sum(y > x)
	df2 <- length(y)
	para <- c(df1,df2)
	attr(para,"names") <- c("num df","denom df")
	msd <- ssd/df1
	emsd <- error^2*4
	estimate <- c(ssd,msd,emsd)
	attr(estimate,"names") <- c("SSD","MSD","EMSD") ## EMSD=EXPECTED MEAN SQUARE DISTANCE
	F <- msd/emsd
	attr(F,"names") <- "F"
	pval <- pf(F,df1, df2, lower.tail=FALSE)
	ans <- list(estimate=estimate,statistic=F,parameter=para,p.value=pval,
	method="Test of Causual Necessity.",data.name=dn)
	class(ans) <- "htest"
	ans
}

fsgof.suff <- function(y, x, damp=.01, error=.05) {
	dn <- sprintf("%s and %s", deparse(substitute(y)),deparse(substitute(x)))
	xZ <- zTransform(x,damp=damp)
	yZ <- zTransform(y,damp=damp)
	d <- as.numeric(y > x)
	ssd <- sum((1-d)*(yZ - xZ)^2)
	df1 <- sum(y < x)
	df2 <- length(y)
	para <- c(df1,df2)
	attr(para,"names") <- c("num df","denom df")
	msd <- ssd/df1
	emsd <- error^2*4
	estimate <- c(ssd,msd,emsd)
	attr(estimate,"names") <- c("SSD","MSD","EMSD") ## EMSD=EXPECTED MEAN SQUARE DISTANCE
	F <- msd/emsd
	attr(F,"names") <- "F"
	pval <- pf(F,df1, df2, lower.tail=FALSE)
	ans <- list(estimate=estimate,
	statistic=F,parameter=para,p.value=pval,
	method="Test of Causual Sufficiency.",data.name=dn)
	class(ans) <- "htest"
	ans
}

fsgof.suffnec <- function(y, x, damp=.01, error=.05) {
	dn <- sprintf("%s and %s", deparse(substitute(y)),deparse(substitute(x)))
	xZ <- zTransform(x,damp=damp)
	yZ <- zTransform(y,damp=damp)
	d <- as.numeric(y > x)
	ssd <- sum(d*(yZ - xZ)^2) + sum((1-d)*(yZ - xZ)^2)
	df1 <- df2 <- length(y)
	para <- c(df1,df2)
	attr(para,"names") <- c("num df","denom df")
	msd <- ssd/df1
	emsd <- error^2*4
	estimate <- c(ssd,msd,emsd)
	attr(estimate,"names") <- c("SSD","MSD","EMSD") ## EMSD=EXPECTED MEAN SQUARE DISTANCE
	F <- msd/emsd
	attr(F,"names") <- "F"
	pval <- pf(F,df1, df2, lower.tail=FALSE)
	ans <- list(estimate=estimate,
	statistic=F,parameter=para,p.value=pval,
	method="Test of Causual Sufficiency and Necessity.",data.name=dn)
	class(ans) <- "htest"
	ans
}

fsgof.test <- function(y, x, type=c("suffnec","suff","nec"), damp=.01, error=.05){
	type <- match.arg(type)
	ans <- switch(type,
		"suffnec"=fsgof.suffnec(y, x, damp=damp, error=error),
 		"nec"=fsgof.nec(y, x, damp=damp, error=error),
		"suff"=fsgof.suff(y, x, damp=damp, error=error)
	)
	ans
}

fsgof.descr <- function(y, x, damp=.01){
	yZ <- zTransform(y,damp=damp)
	xZ <- zTransform(x,damp=damp)
	d <- as.numeric(y > x)
	Dnull <- sum((yZ -mean(yZ))^2)
	Dnec <- sum(d*(yZ - xZ)^2)
	Dsuff <- sum((1-d)*(yZ - xZ)^2)
	RDnull <- Dnull/(Dnull+Dnec+Dsuff)
	RDnec <- Dnec/(Dnull+Dnec+Dsuff)
	RDsuff <- Dsuff/(Dnull+Dnec+Dsuff)
	RDsuffnec <- (Dnec+Dsuff)/(Dnull+Dnec+Dsuff)
	Rnull <- 1-RDnull
	Rnec<- 1-RDnec
	Rsuff <- 1-RDsuff
	Rsuffnec <- 1-RDsuffnec
	relative.distance <- data.frame(null=RDnull,nec=RDnec,suffnec=RDsuffnec)
	relative.consistency <- data.frame(null=Rnull,nec=Rnec,suffnec=Rsuffnec)
	ans <- rbind(relative.distance,relative.consistency)
	row.names(ans) <- c("distance","consistency")
	ans
}

##y=c(0.8,0.5,0.5,0.4); x=c(0.7,0.3,0.4,0.45);x2 =c(0.66,0.55,0.22,0.52)
#fsgof.nec(y,x)
#fsgof.suff(y,x)
#fsgof.suffnec(y,x)
#fsgof.descr(y,x)
#high order test
#fsgof.suffnec(y,pmin(x,x2))

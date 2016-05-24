

# Lui, Mayer, Eckhardt: CI for the risk ratio under cluster sampling based on the beta-binomial model
# Stat i Med. 2000, 19, 2933-2942.
# Inflation factor due to intraclass correlation
# page 2934, section 2.1, paragraph 1

Amrho <- function(mij, rhoi, llod=NULL)
{
misum <- sum(mij)
sum( mij*(1 + (mij - 1)*rhoi)/misum )
}

# page 2934/2935:

estrhoi <- function(yij, mij)
{
ni <- length(mij)
misum <- sum(mij)
BMSi <- (sum((yij^2)/mij) - (sum(yij)^2)/misum)/(ni - 1)
WMSi <- (sum(yij) - sum((yij^2)/mij))/(sum(mij-1))
mistar <- ((misum^2) - sum(mij^2))/((ni-1)*misum)
estimate <- (BMSi-WMSi)/(BMSi + (mistar-1)*WMSi)
if(is.na(estimate)){return(0)}else{return(estimate)}
}


# Asymptotic interval on the log-scale (Following Katz et al.) 
# Eq.(2), page 2935,

# iccpool: 	= FALSE: estimate the the intra-class- correlation separately per group
#		= TRUE: pool the estimates for the the intra-class- correlation between groups

# resbin: 	= FALSE: the overdispersion parameter is estimated from the data without a lower bound, thus underdispersion may be estimated
#		= TRUE: lower limit for the overdispersion parameter is set to 1,
#		 the overdispersion is estimated from the data, but the variance can not fall below the binomial variance


RRasylogBB <- function(y1, y2, m1, m2, conf.level=0.95, alternative="two.sided", iccpool=FALSE, resbin=FALSE)
{

# sparse data adjustment (p.2936, paragraph following Eq(6))

y1sum <- sum(y1)
y2sum <- sum(y2)

m1sum <- sum(m1)
m2sum <- sum(m2)

if(y1sum==0 | y1sum==m1sum){pi1hat <- (y1sum+0.5)/(m1sum+1)}else{pi1hat <- y1sum/m1sum}
if(y2sum==0 | y2sum==m2sum){pi2hat <- (y2sum+0.5)/(m2sum+1)}else{pi2hat <- y2sum/m2sum}


# estimate intraclass correlation

rho1hatI <- estrhoi(yij=y1, mij=m1)
rho2hatI <- estrhoi(yij=y2, mij=m2)

# common estimate for the parameter rho (pooled over the two groups)
if(iccpool){rho1hat <- rho2hat <- (rho1hatI*sum(m1) + rho2hatI*sum(m2))/(sum(m1)+sum(m2))}else{rho1hat<-rho1hatI; rho2hat<-rho2hatI}
if(resbin){
varlogRR <- ((1-pi1hat)*max(Amrho(mij=m1, rhoi=rho1hat),1)/(sum(m1)*pi1hat) + (1-pi2hat)*max(Amrho(mij=m2, rhoi=rho2hat),1)/(sum(m2)*pi2hat))
}else{
varlogRR <- ((1-pi1hat)*Amrho(mij=m1, rhoi=rho1hat)/(sum(m1)*pi1hat) + (1-pi2hat)*Amrho(mij=m2, rhoi=rho2hat)/(sum(m2)*pi2hat))
}

estimate <- pi1hat/pi2hat

switch(alternative,

"two.sided"={
quant <- qnorm(p = 1 - (1-conf.level)/2)
lower <- exp(log(estimate) - quant*sqrt(varlogRR))
upper <- exp(log(estimate) + quant*sqrt(varlogRR))
},

"less"={
quant <- qnorm(p = conf.level )
lower <- 0
upper <- exp(log(estimate) + quant*sqrt(varlogRR))
},

"greater"={
quant <- qnorm(p = conf.level)
lower <- exp(log(estimate) - quant*sqrt(varlogRR))
upper <- Inf
}
)

ICCest <- c(rho1hat, rho2hat)
out <- list(conf.int=c(lower=lower, upper=upper), estimate=estimate, esticc=ICCest)

return(out)
}

###############


###########################################


RRFiellerBaileyBB <- function(y1, y2, m1, m2, conf.level=0.95, alternative="two.sided", iccpool=FALSE, resbin=FALSE)
{

# Fieller-type with skewcorrection accord. to Bailey
# p 2936, second paragraph and Eq.(5) ff

y1sum <- sum(y1)
y2sum <- sum(y2)

m1sum <- sum(m1)
m2sum <- sum(m2)

if(y1sum==0 | y1sum==m1sum){pi1hat <- (y1sum+0.5)/(m1sum+1)}else{pi1hat <- y1sum/m1sum}
if(y2sum==0 | y2sum==m2sum){pi2hat <- (y2sum+0.5)/(m2sum+1)}else{pi2hat <- y2sum/m2sum}


rho1hatI <- estrhoi(yij=y1, mij=m1)
rho2hatI <- estrhoi(yij=y2, mij=m2)

if(iccpool){rho1hat <- rho2hat <- (rho1hatI*sum(m1) + rho2hatI*sum(m2))/(sum(m1)+sum(m2))}else{rho1hat<-rho1hatI; rho2hat<-rho2hatI}

estimate <- pi1hat/pi2hat

switch(alternative,

"two.sided"={

quant <- qnorm(p = 1 - (1-conf.level)/2)

if(resbin){
Aterm <- pi2hat^(2/3) - (quant^2) * (1-pi2hat) * max(Amrho(mij=m2, rhoi=rho2hat),1)/(9*m2sum*pi2hat^(1/3))
Bterm <- (pi1hat*pi2hat)^(1/3)
Cterm <- pi1hat^(2/3) - (quant^2) * (1-pi1hat) * max(Amrho(mij=m1, rhoi=rho1hat),1)/(9*m1sum*pi1hat^(1/3))
}else{
Aterm <- pi2hat^(2/3) - (quant^2) * (1-pi2hat) * Amrho(mij=m2, rhoi=rho2hat)/(9*m2sum*pi2hat^(1/3))
Bterm <- (pi1hat*pi2hat)^(1/3)
Cterm <- pi1hat^(2/3) - (quant^2) * (1-pi1hat) * Amrho(mij=m1, rhoi=rho1hat)/(9*m1sum*pi1hat^(1/3))
}


if((!any(is.nan(c(Aterm, Bterm, Cterm)))) && Aterm>0 & (Bterm^2 - Aterm*Cterm)>0)
{
lower <- max(0, ((Bterm - sqrt(Bterm^2 - Aterm*Cterm))/Aterm)^3)
upper <- ((Bterm + sqrt(Bterm^2 - Aterm*Cterm))/Aterm)^3
}else{lower<-0; upper<-Inf} 

},

"less"={

quant <- qnorm(p = conf.level)

if(resbin){
Aterm <- pi2hat^(2/3) - (quant^2) * (1-pi2hat) * max(Amrho(mij=m2, rhoi=rho2hat),1)/(9*m2sum*pi2hat^(1/3))
Bterm <- (pi1hat*pi2hat)^(1/3)
Cterm <- pi1hat^(2/3) - (quant^2) * (1-pi1hat) * max(Amrho(mij=m1, rhoi=rho1hat),1)/(9*m1sum*pi1hat^(1/3))
}else{
Aterm <- pi2hat^(2/3) - (quant^2) * (1-pi2hat) * Amrho(mij=m2, rhoi=rho2hat)/(9*m2sum*pi2hat^(1/3))
Bterm <- (pi1hat*pi2hat)^(1/3)
Cterm <- pi1hat^(2/3) - (quant^2) * (1-pi1hat) * Amrho(mij=m1, rhoi=rho1hat)/(9*m1sum*pi1hat^(1/3))
}

if((!any(is.nan(c(Aterm, Bterm, Cterm)))) && Aterm>0 & (Bterm^2 - Aterm*Cterm)>0)
{
lower <- 0
upper <- ((Bterm + sqrt(Bterm^2 - Aterm*Cterm))/Aterm)^3
}else{lower<-0; upper<-Inf} 

},

"greater"={

quant <- qnorm(p = conf.level)
if(resbin){
Aterm <- pi2hat^(2/3) - (quant^2) * (1-pi2hat) * max(Amrho(mij=m2, rhoi=rho2hat),1)/(9*m2sum*pi2hat^(1/3))
Bterm <- (pi1hat*pi2hat)^(1/3)
Cterm <- pi1hat^(2/3) - (quant^2) * (1-pi1hat) * max(Amrho(mij=m1, rhoi=rho1hat),1)/(9*m1sum*pi1hat^(1/3))
}else{
Aterm <- pi2hat^(2/3) - (quant^2) * (1-pi2hat) * Amrho(mij=m2, rhoi=rho2hat)/(9*m2sum*pi2hat^(1/3))
Bterm <- (pi1hat*pi2hat)^(1/3)
Cterm <- pi1hat^(2/3) - (quant^2) * (1-pi1hat) * Amrho(mij=m1, rhoi=rho1hat)/(9*m1sum*pi1hat^(1/3))
}

if( (!any(is.nan(c(Aterm, Bterm, Cterm)))) && Aterm>0 & (Bterm^2 - Aterm*Cterm)>0)
{
lower <- max(0, ((Bterm - sqrt(Bterm^2 - Aterm*Cterm))/Aterm)^3)
upper <- Inf
}else{lower<-0; upper<-Inf} 

}
)

ICCest <- c(rho1hat, rho2hat)
out <- list(conf.int=c(lower=lower, upper=upper), estimate=estimate, esticc=ICCest)

return(out)

}



# II) Implementing 
# Zaihra, T and Paul, S 2010: Interval Estimation of Some Epidemiological
# Measures of Association. The International Journal of Biostatistics. 6 (1), Article 35.

# Eq. (1)

variOD <- function(mij, yij, add=TRUE, resbin=FALSE)
{
ni <- length(mij)
misum <- sum(mij)
if(add){pihat <- (sum(yij)+0.5)/(misum+1)}else{pihat <- sum(yij)/misum}

rij <- yij - pihat*mij 

vi <- ni/(ni-1)*sum((rij^2)/misum^2)

if(resbin){
vibin<-pihat*(1-pihat)/misum
if(vi < vibin){vi<-vibin}
}

vSi <- ni/(ni-1)*vi

return(c(vi, vSi))

}


RRasylogOD <- function(y1, y2, m1, m2, conf.level=0.95, alternative="two.sided", varmethod="res", resbin=FALSE)
{

# sparse data adjustment (p.2936, paragraph following Eq(6))

y1sum <- sum(y1)
y2sum <- sum(y2)

m1sum <- sum(m1)
m2sum <- sum(m2)

if(y1sum==0 | y1sum==m1sum){pi1hat <- (y1sum+0.5)/(m1sum+1)}else{pi1hat <- y1sum/m1sum}
if(y2sum==0 | y2sum==m2sum){pi2hat <- (y2sum+0.5)/(m2sum+1)}else{pi2hat <- y2sum/m2sum}


vi1hat <- variOD(mij=m1, yij=y1, resbin=resbin)
vi2hat <- variOD(mij=m2, yij=y2, resbin=resbin)


switch(varmethod,
"res"={varlogRRhat <- vi1hat[1] / pi1hat^2 + vi2hat[1] / pi2hat^2; ODest <- c(vi1hat[1], vi2hat[1])},
"san"={varlogRRhat <- vi1hat[2] / pi1hat^2 + vi2hat[2] / pi2hat^2; ODest <- c(vi1hat[2], vi2hat[2])})

estimate <- pi1hat/pi2hat

switch(alternative,

"two.sided"={
quant <- qnorm(p = 1 - (1-conf.level)/2)
lower <- estimate * exp( - quant*sqrt(varlogRRhat))
upper <- estimate * exp( + quant*sqrt(varlogRRhat))
},

"less"={
quant <- qnorm(p = conf.level )
lower <- 0
upper <- estimate * exp( + quant*sqrt(varlogRRhat))
},

"greater"={
quant <- qnorm(p = conf.level)
lower <- estimate * exp( - quant*sqrt(varlogRRhat))
upper <- Inf
}
)

# varmethod="res" corresponds to method MR2, using the residuals variance estimator
# varmethod="san" corresponds to method MS2, using the sandwich variance estimator

out <- list(conf.int=c(lower=lower, upper=upper), estimate=estimate, estOD=ODest)

return(out)
}


RRFiellerBaileyOD <- function(y1, y2, m1, m2, conf.level=0.95, alternative="two.sided", varmethod="res", resbin=FALSE)
{

# varmethod="res" corresponds to method MR4, using the residuals variance estimator
# varmethod="san" corresponds to method MS4, using the sandwich variance estimator

y1sum <- sum(y1)
y2sum <- sum(y2)

m1sum <- sum(m1)
m2sum <- sum(m2)

if(y1sum==0 | y1sum==m1sum){pi1hat <- (y1sum+0.5)/(m1sum+1)}else{pi1hat <- y1sum/m1sum}
if(y2sum==0 | y2sum==m2sum){pi2hat <- (y2sum+0.5)/(m2sum+1)}else{pi2hat <- y2sum/m2sum}

estimate <- pi1hat/pi2hat


switch(varmethod,
"res"={varpi1hat <- variOD(mij=m1, yij=y1, resbin=resbin)[1]; varpi2hat <- variOD(mij=m2, yij=y2, resbin=resbin)[1]; ODest <- c(varpi1hat, varpi2hat)},
"san"={varpi1hat <- variOD(mij=m1, yij=y1, resbin=resbin)[2]; varpi2hat <- variOD(mij=m2, yij=y2, resbin=resbin)[2]; ODest <- c(varpi1hat, varpi2hat)})


# varmethod="res" corresponds to method MR4, using the residuals variance estimator
# varmethod="san" corresponds to method MS4, using the sandwich variance estimator



switch(alternative,

"two.sided"={

quant <- qnorm(p = 1 - (1-conf.level)/2)

Aterm <- pi2hat^(2/3) - (quant^2) * varpi2hat/(9*pi2hat^(4/3))
Bterm <- (pi1hat*pi2hat)^(1/3)
Cterm <- pi1hat^(2/3) - (quant^2) * varpi1hat/(9*pi1hat^(4/3))

if((!any(is.nan(c(Aterm, Bterm, Cterm)))) && Aterm>0 & (Bterm^2 - Aterm*Cterm)>0)
{
lower <- max(0, ((Bterm - sqrt(Bterm^2 - Aterm*Cterm))/Aterm)^3)
upper <- ((Bterm + sqrt(Bterm^2 - Aterm*Cterm))/Aterm)^3
}else{lower<-0; upper<-Inf} 

},

"less"={

quant <- qnorm(p = conf.level)

Aterm <- pi2hat^(2/3) - (quant^2) * varpi2hat/(9*pi2hat^(4/3))
Bterm <- (pi1hat*pi2hat)^(1/3)
Cterm <- pi1hat^(2/3) - (quant^2) * varpi1hat/(9*pi1hat^(4/3))

if((!any(is.nan(c(Aterm, Bterm, Cterm))))&& Aterm>0 && (Bterm^2 - Aterm*Cterm)>0)
{
lower <- 0
upper <- ((Bterm + sqrt(Bterm^2 - Aterm*Cterm))/Aterm)^3
}else{lower<-0; upper<-Inf} 

},

"greater"={

quant <- qnorm(p = conf.level)

Aterm <- pi2hat^(2/3) - (quant^2) * varpi2hat/(9*pi2hat^(4/3))
Bterm <- (pi1hat*pi2hat)^(1/3)
Cterm <- pi1hat^(2/3) - (quant^2) * varpi1hat/(9*pi1hat^(4/3))

if((!any(is.nan(c(Aterm, Bterm, Cterm)))) && Aterm>0 && (Bterm^2 - Aterm*Cterm)>0)
{
lower <- max(0, ((Bterm - sqrt(Bterm^2 - Aterm*Cterm))/Aterm)^3)
upper <- Inf
}else{lower<-0; upper<-Inf} 

}
)

out <- list(conf.int=c(lower=lower, upper=upper), estimate=estimate, estOD=ODest)
return(out)

}



# x and y must be two data sets or matrices with two columns
# containing integers
# first column is treated as success
# second is treated as failure
# row sums of matrix will be the pool sizes (m)
# numer of rows: clusters

ODbinChecksamples <- function(x, y)
{
if(!is.data.frame(x) & !is.matrix(x)){stop("x must be a data.frame or a matrix")}
if(!is.data.frame(y) & !is.matrix(y)){stop("y must be a data.frame or a matrix")}

x <- as.matrix(x)
y <- as.matrix(y)

if(ncol(x)!=2){stop("x must have exactly 2 columns")}
if(ncol(y)!=2){stop("y must have exactly 2 columns")}

if(is.data.frame(x)){x <- as.matrix(x)}
if(is.data.frame(y)){y <- as.matrix(y)}

if(!(is.numeric(x)|is.integer(x))){stop("x must have integer (or whole number) entries")}
if(!(is.numeric(y)|is.integer(y))){stop("y must have integer (or whole number) entries")}

is.wholenumber <- function(x, tol = .Machine$double.eps^0.5){abs(x - round(x)) < tol}

if(!all(is.wholenumber(x))){warning("At least one entry in sample x is not a whole number!")}
if(!all(is.wholenumber(y))){warning("At least one entry in sample y is not a whole number!")}

if(is.null(colnames(x))){colnames(x)<-c("c1","c2")}
if(is.null(colnames(y))){colnames(y)<-c("c1","c2")}

cnx <- colnames(x); cny <- colnames(y)
pdefx <- paste("Proportion is defined as ", cnx[1], "/(", paste(cnx, collapse="+"), ")", sep="") 
pdefy <- paste("Proportion is defined as ", cny[1], "/(", paste(cny, collapse="+"), ")", sep="") 

y1 <- x[,1]
m1 <- rowSums(x)

y2 <- y[,1]
m2 <- rowSums(y)

return(list(y1=y1, y2=y2, m1=m1, m2=m2, sf1 = x[,1:2], sf2 = y[,1:2], pdef1=pdefx, pdef2=pdefy))
}


ODbin.ratio <- function(x, y, conf.level=0.95, alternative="two.sided", CImethod=c("FOD", "LOD"), varmethod=c("res", "san"), resbin=FALSE) 
{

xylist <- ODbinChecksamples(x=x, y=y)

if(!is.numeric(conf.level) || length(conf.level)>1 || c(conf.level<=0 | conf.level>=1)){stop("'conf.level' must be a single numeric value between 0 and 1")}
alternative <- match.arg(alternative, choices=c("two.sided", "less", "greater"))
CImethod <- match.arg(CImethod)
varmethod <- match.arg(varmethod)
if(!is.logical(resbin)){stop("Argument 'resbin' must be TRUE or FALSE")}


switch(CImethod, 
"LOD"={
OUT <- RRasylogOD(y1=xylist$y1, y2=xylist$y2, m1=xylist$m1, m2=xylist$m2,
 conf.level=conf.level, alternative=alternative, varmethod=varmethod, resbin=resbin)
mmeth <- "Delta-method, log-scale"
},
"FOD"={OUT <- RRFiellerBaileyOD(y1=xylist$y1, y2=xylist$y2, m1=xylist$m1, m2=xylist$m2,
 conf.level=conf.level, alternative=alternative, varmethod=varmethod, resbin=resbin)
mmeth <- "Fieller-Bailey-method"
})

switch(varmethod, "res"={mvar<-"residual variance"}, "san"={mvar<-"sandwich variance"})
if(resbin){mresbin <- ", variance limited to binomial variance"}else{mresbin<-""}

METH <- paste(mmeth,", overdispersed binomial, ", mvar, mresbin, ".", sep="")

attr(OUT$conf.int, which="methodname") <- METH

OUT$pdefx <- xylist$pdefx
OUT$pdefy <- xylist$pdefy

return(OUT)

}



Betabin.ratio <- function(x, y, conf.level=0.95, alternative="two.sided", CImethod=c("FBB", "LBB"), iccpool=FALSE, resbin=FALSE) 
{

xylist <- ODbinChecksamples(x=x, y=y)

if(!is.numeric(conf.level) || length(conf.level)>1 || c(conf.level<=0 | conf.level>=1)){stop("'conf.level' must be a single numeric value between 0 and 1")}
alternative <- match.arg(alternative, choices=c("two.sided", "less", "greater"))
CImethod <- match.arg(CImethod)

if(!is.logical(iccpool)){stop("Argument 'iccpool' must be TRUE or FALSE")}
if(!is.logical(resbin)){stop("Argument 'resbin' must be TRUE or FALSE")}

switch(CImethod, 
"LBB"={
OUT <- RRasylogBB(y1=xylist$y1, y2=xylist$y2, m1=xylist$m1, m2=xylist$m2,
 conf.level=conf.level, alternative=alternative, iccpool=iccpool, resbin=resbin)
mmeth <- "Delta-method, log-scale"
},
"FBB"={OUT <- RRFiellerBaileyBB(y1=xylist$y1, y2=xylist$y2, m1=xylist$m1, m2=xylist$m2,
 conf.level=conf.level, alternative=alternative, iccpool=iccpool, resbin=resbin)
mmeth <- "Fieller-Bailey-method"
})

if(iccpool){mvar<-"icc estimated pooled over samples"}else{mvar<-"icc estimated separately for each sample"}
if(resbin){mresbin <- ", variance limited to binomial variance"}else{mresbin<-""}

METH <- paste(mmeth,", betabinomial, ", mvar, mresbin, ".", sep="")

attr(OUT$conf.int, which="methodname") <- METH

OUT$pdefx <- xylist$pdefx
OUT$pdefy <- xylist$pdefy

return(OUT)
}


############################################################

# MOVER-R method for quasibinomial GLM



MOVERR <- function(theta0, ci0, theta1, ci1, alternative="two.sided")
{
#METHOD <- "Method of variance estimates recovery (Donner, Zou, 2012)"

estimate <- (theta1)/(theta0)

# Eq. (9), Donner and Zou, 2012, Stat Methods Med Res 2012 21: 347-359.

switch(alternative,
"two.sided"={
lower <- (theta1*theta0 - sqrt( (theta1*theta0)^2 - ci1[1]*ci0[2]*(2*theta1 - ci1[1])*(2*theta0-ci0[2])))/(ci0[2]*(2*theta0 - ci0[2]))
if(theta0==0){ upper<-Inf}else{
upper <- (theta1*theta0 + sqrt( (theta1*theta0)^2 - ci1[2]*ci0[1]*(2*theta1 - ci1[2])*(2*theta0-ci0[1])))/(ci0[1]*(2*theta0 - ci0[1]))}
},

"greater"={
if(theta1==0){ lower <- 0}else{
lower <- (theta1*theta0 - sqrt( (theta1*theta0)^2 - ci1[1]*ci0[2]*(2*theta1 - ci1[1])*(2*theta0-ci0[2])))/(ci0[2]*(2*theta0 - ci0[2]))}
upper <- Inf
},

"less"={
lower <- 0
if(theta0==0){ upper <- Inf}else{
upper <- (theta1*theta0 + sqrt( (theta1*theta0)^2 - ci1[2]*ci0[1]*(2*theta1 - ci1[2])*(2*theta0-ci0[1])))/(ci0[1]*(2*theta0 - ci0[1]))}
})

conf.int<-c(lower=lower, upper=upper)

return(list(conf.int=conf.int,
estimate=estimate,
theta=c(theta1=theta1, theta0=theta0),
ci1=ci1,
ci0=ci0)) 
}

#####################################

# alternative methods for profiling


# signed root deviance over a pre-specified grid
# one factorial binomial glm with logit link

profilebin1 <- function(succ, fail, trt, grid=seq(-20,20,0.5) )
{
FIT1<-glm(cbind(succ, fail) ~ 0 + trt, data=trt, family=binomial(link="logit"))
DEV1<-deviance(FIT1)
n<-table(trt)
PARA<-coef(FIT1)
srdlist<-list()

# profiling along a fixed grid

for(i in seq(along.with=PARA)){
SRDi<-numeric(length=length(grid))
for(j in seq(along.with=grid)){
offj<-PARA
offj[i]<-grid[j]
offjv<-rep(offj, n)
sig<-sign(grid[j]-PARA[i])
FITO<-glm(cbind(succ, fail) ~0+offset(offjv), data=trt, family=binomial())
SRDi[j]<-sqrt(deviance(FITO)-DEV1)*sig
}
srdlist[[i]]<-SRDi}

names(srdlist)<-names(PARA)
return(list(srdlist=srdlist, grid=grid, dfres=FIT1$df.residual, estimate=PARA))
}

# signed root deviance over a pre-specified grid
# one factorial quasibinomial glm with logit link


profileqbin1<-function(succ, fail, trt, grid=seq(-20,20,0.5) )
{
DAT<-data.frame(succ=succ, fail=fail, trt=trt)
FIT1<-glm(cbind(succ, fail) ~ 0 + trt, data=DAT, family=quasibinomial())
SUMM<-summary(FIT1)
DEV1<-deviance(FIT1)
DISP<-SUMM$dispersion
n<-table(trt)
PARA<-coef(FIT1)
srdlist<-list()

for(i in seq(along.with=PARA)){
SRDi<-numeric(length=length(grid))
for(j in seq(along.with=grid)){
offj<-PARA
offj[i]<-grid[j]
offjv<-rep(offj, n)
sig<-sign(grid[j]-PARA[i])
DATij <- DAT
DATij$offjv
FITO <- glm(cbind(succ, fail) ~0+offset(offjv), data=DATij, family=quasibinomial())
SRDi[j] <- sig*sqrt((deviance(FITO)-DEV1)/DISP)
}
srdlist[[i]]<-SRDi}
names(srdlist)<-names(PARA)
return(list(srdlist=srdlist, grid=grid, dfres=FIT1$df.residual, disp=DISP, estimate=PARA))
}


profileci <- function(srdlist, grid, quant)
{
npar<-length(srdlist)
CI<-matrix(NA, ncol=length(quant), nrow=npar)

for(i in 1:npar)
{SPL<-spline(x=grid, y=srdlist[[i]])
CI[i,] <- approx(x=SPL$y, y=SPL$x, xout=quant)$y}

return(CI)
}


# profilebin1, profileqbin1, profileci

etabinci <- function(succ, fail, trt, grid=NULL, alternative="two.sided", conf.level=0.95)
{
if(is.null(grid)){grid<-c(seq(-20,-5, 0.5), seq(-4.9,4.9,0.1), seq(5,20,0.5))}

alternative<-match.arg(alternative, choices=c("two.sided","less","greater"))

PROF <- profilebin1(succ=succ, fail=fail, trt=trt, grid=grid )

switch(alternative,
"two.sided"={
QUANT <- qnorm(p=c((1-conf.level)/2, 1-(1-conf.level)/2))
CI <- profileci(srdlist=PROF$srdlist, grid=PROF$grid, quant=QUANT)
},
"less"={
QUANT <- qnorm(p= 1-(1-conf.level))
CIu <- profileci(srdlist=PROF$srdlist, grid=PROF$grid, quant=QUANT)
CI <- cbind(0, CIu)
},
"greater"={
QUANT <- qnorm(p= (1-conf.level))
CIl <- profileci(srdlist=PROF$srdlist, grid=PROF$grid, quant=QUANT)
CI <- cbind(CIl, Inf)
})

colnames(CI)<-c("lower","upper")

return(list(conf.int=CI, estimate=PROF$estimate, quant=QUANT, conf.level=conf.level, alternative=alternative, profile=PROF))
}


etaqbinci <- function(succ, fail, trt, grid=NULL, alternative="two.sided", conf.level=0.95)
{
if(is.null(grid)){grid<-c(seq(-20,-5, 0.5), seq(-4.9,4.9,0.1), seq(5,20,0.5))}

alternative<-match.arg(alternative, choices=c("two.sided","less","greater"))

PROF <- profileqbin1(succ=succ, fail=fail, trt=trt, grid=grid )

switch(alternative,
"two.sided"={
QUANT <- qt(p=c((1-conf.level)/2, 1-(1-conf.level)/2), df=PROF$dfres)
CI <- profileci(srdlist=PROF$srdlist, grid=PROF$grid, quant=QUANT)
},
"less"={
QUANT <- qt(p= 1-(1-conf.level), df=PROF$dfres)
CIu <- profileci(srdlist=PROF$srdlist, grid=PROF$grid, quant=QUANT)
CI <- cbind(0, CIu)
},
"greater"={
QUANT <- qt(p= (1-conf.level), df=PROF$dfres)
CIl <- profileci(srdlist=PROF$srdlist, grid=PROF$grid, quant=QUANT)
CI <- cbind(CIl, Inf)
})

colnames(CI)<-c("lower","upper")

return(list(conf.int=CI, estimate=PROF$estimate, quant=QUANT, conf.level=conf.level, alternative=alternative, profile=PROF))
}

#######################################


# adapted QBmover function: dealing with extreme cases

# mcprofile
# etabinci, etaqbinci

QBmover <- function(succ, fail, trt, conf.level=0.95, alternative="two.sided", grid=NULL)
{

trt<-as.factor(trt)
ntrt<-length(levels(trt))
if(ntrt<2){stop("At least two levels in trt needed!")}

excesss<-tapply(succ, trt, function(x){all(x[1]==x)})
excessf<-tapply(fail, trt, function(x){all(x[1]==x)})

expit<-function(x){exp(x)/(1+exp(x))}

if(is.null(grid)){grid<-c(seq(-10,-4,0.5), seq(-3.8,3.8,0.2), seq(4,10,0.5))}

if(any(excesss|excessf)){

FIT <- glm(cbind(succ, fail) ~ 0+trt, family=quasibinomial())

if(summary(FIT)$dispersion<1){
CI<-etabinci(succ=succ, fail=fail, trt=trt, grid=grid, alternative=alternative, conf.level=conf.level); resbin <- TRUE}else{
CI<-etaqbinci(succ=succ, fail=fail, trt=trt, grid=grid, alternative=alternative, conf.level=conf.level); resbin <- FALSE}

esti<-expit(CI$estimate)
oci<-expit(CI$conf.int)

oci[is.na(oci[,1]),1]<-0
oci[is.na(oci[,2]),2]<-1

}else{

# Using mcprofile

FIT <- glm(cbind(succ, fail) ~ 0+trt, family=quasibinomial())
resbin <- FALSE
if(summary(FIT)$dispersion<1){
resbin <- TRUE
FIT <- glm(cbind(succ, fail) ~ 0+trt, family=binomial())}

CM=diag(rep(1, ntrt))
GR=matrix( rep( grid, times=ntrt), ncol=ntrt)
PROF <- mcprofile(FIT, CM=CM, grid=GR)

switch(alternative,
"two.sided"={CI <- confint(PROF, adjust="none", level=conf.level)},
"less"={CI <- confint(PROF, adjust="none", level=1-(1-conf.level)*2)},
"greater"={CI <- confint(PROF, adjust="none", level=1-(1-conf.level)*2)})

esti <- expit(CI$estimate[,1])

oci <- expit(CI$confint)

oci[is.na(oci[,1]),1]<-0
oci[is.na(oci[,2]),2]<-1

# end of using mcprofile
}

MOVERCI<-list()
NAM<-character()
trtl<-levels(trt)
for(i in 1:(ntrt-1))
{
TEMP <- MOVERR(theta0=esti[1], ci0=oci[1,], theta1=esti[(i+1)], ci1=oci[(i+1),], alternative=alternative)
OUT <- as.numeric(c(TEMP$estimate, TEMP$conf.int))
names(OUT) <- c("est","lower","upper")
MOVERCI[[i]] <- OUT
NAM[i]<-paste(trtl[(i+1)], trtl[1], sep="/")
}

names(MOVERCI)<-NAM

OUT<-t(data.frame(MOVERCI))
attr(OUT, which="dispersion") <- summary(FIT)$dispersion
attr(OUT, which="resbin") <- resbin

nal<-(is.na(OUT[,2]) | is.nan(OUT[,2]))
nau<-(is.na(OUT[,3]) | is.nan(OUT[,3]))
if(any(nal)){OUT[nal,2]<-0}
if(any(nau)){OUT[nau,3]<-Inf}

return(OUT)

}

Quasibin.ratio <- function(x, y, conf.level=0.95, alternative="two.sided", grid=NULL) 
{

xylist <- ODbinChecksamples(x=x, y=y)

dxy <- rbind(data.frame(xylist$sf2, trt="y"),data.frame(xylist$sf1, trt="x"))
dxy$trt <- factor(dxy$trt, levels=c("y", "x"))

if(!is.numeric(conf.level) || length(conf.level)>1 || c(conf.level<=0 | conf.level>=1)){stop("'conf.level' must be a single numeric value between 0 and 1")}
alternative <- match.arg(alternative, choices=c("two.sided", "less", "greater"))

# Note: problem cbind(succ, fail)~ causes, that first variable in model framne is matrix with two columns, thus:

OUT <- QBmover(succ=dxy[,1], fail=dxy[,2], trt=dxy[,"trt"], conf.level=conf.level, alternative=alternative, grid=grid)

if(attr(OUT, which="resbin")){
estdisp <- 1
attr(estdisp, which="estimation") <- "dispersion restricted to be 1 (binomial variance)"
}else{
estdisp <- attr(OUT, which="dispersion")
attr(estdisp, which="estimation") <-  paste("estimated dispersion parameter= ", format(estdisp, digits=2) , sep="")}

METH <- "MOVER-R method applied to profile-deviance CIs, (quasi-)binomial assumption"

RES<-list()
RES$conf.int <- OUT[,c("lower", "upper")]
attr(RES$conf.int, which="methodname") <- METH
RES$estimate <- OUT[,"est"]
RES$pdefx <- xylist$pdefx
RES$pdefy <- xylist$pdefy
RES$estdisp <- estdisp

return(RES)
}




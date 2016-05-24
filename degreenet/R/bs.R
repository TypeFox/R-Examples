#  File degreenet/R/bs.R
#  Part of the statnet package, http://statnet.org
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) in
#    http://statnet.org/attribution
#
# Copyright 2003 Mark S. Handcock, University of California-Los Angeles
# Copyright 2007 The statnet Development Team
######################################################################
#
# Bootstrap CI
#
bootstrapzipf <- function(x,cutoff=1,
                          m=200,alpha=0.95){
aaa <- adpmle(x=x,cutoff=cutoff)
bmles <- rep(0,length=m)
for(i in seq(along=bmles)){
 xsamp <- sample(x,size=length(x),replace=TRUE)
 bmles[i] <- adpmle(x=xsamp,cutoff=cutoff)[3]
 print(i)
}
#
c(quantile(bmles,c(0.5*(1-alpha),0.5,0.5+0.5*alpha),na.rm=TRUE),aaa)
}
bootstrapdp <- function(x,cutoff=1,cutabove=1000,
                          m=200,alpha=0.95,guess=3.31,hellinger=FALSE,
                          mle.meth="adpmle"){
if(mle.meth!="adpmlef"){
 aaa <- adpmle(x=x,cutoff=cutoff,guess=guess)$theta
}else{
 aaa <- adpmlef(x=x,cutoff=cutoff,guess=guess)$theta
}
bmles <- rep(0,length=m)
for(i in seq(along=bmles)){
 xsamp <- sample(x,size=length(x),replace=TRUE)
 if(mle.meth!="adpmlef"){
  bmles[i] <- adpmle(x=xsamp,cutoff=cutoff,guess=guess)$theta
 }else{
  bmles[i] <- adpmlef(x=xsamp,cutoff=cutoff,guess=guess)$theta
 }
}
#
bmles <- bmles[!is.na(bmles)]
c(quantile(bmles,c(0.5*(1-alpha),0.5,0.5+0.5*alpha),na.rm=TRUE),aaa)
}
#
# Bootstrap CI for Yule
#
bootstrapyule <- function(x,cutoff=1,cutabove=1000,
                          m=200,alpha=0.95,guess=3.31,hellinger=FALSE,
                          mle.meth="ayulemle"){
if(mle.meth!="ayulemlef"){
 aaa <- ayulemle(x=x,cutoff=cutoff,guess=guess,hellinger=hellinger)$theta
}else{
 aaa <- ayulemlef(x=x,cutoff=cutoff,guess=guess,hellinger=hellinger)$theta
}
bmles <- rep(0,length=m)
for(i in seq(along=bmles)){
 xsamp <- sample(x,size=length(x),replace=TRUE)
 if(mle.meth!="ayulemlef"){
  bmles[i] <- ayulemle(x=xsamp,cutoff=cutoff,guess=guess,hellinger=hellinger)$theta
 }else{
  bmles[i] <- ayulemlef(x=xsamp,cutoff=cutoff,guess=guess,hellinger=hellinger)$theta
 }
}
#
bmles <- bmles[!is.na(bmles)]
c(quantile(bmles,c(0.5*(1-alpha),0.5,0.5+0.5*alpha),na.rm=TRUE),aaa)
}
#
# Bootstrap CI for Waring
#
bootstrapwar <- function(x,cutoff=1,cutabove=1000,
                          m=200,alpha=0.95,guess=c(3.31,0.1),
                          file="none", conc=FALSE){
mle <- awarmle(x=x,cutoff=cutoff,guess=guess)$theta
bmles <- rep(0,length=m)
for(i in seq(along=bmles)){
 xsamp <- sample(x,size=length(x),replace=TRUE)
 aaa <- awarmle(x=xsamp,cutoff=cutoff,guess=guess)
 if(conc){
  bmles[i] <- aaa$conc
 }else{
  bmles[i] <- aaa$theta[1]
 }
}
#
if(!missing(file)){
 save(mle,bmles,file=file)
}
c(quantile(bmles,c(0.5*(1-alpha),(1-alpha),0.5,alpha,0.5+0.5*alpha),na.rm=TRUE),mle)
}
#
# Bootstrap CI for Waring Conc
#
bootstrapwarconc <- function(x,cutoff=1,cutabove=1000,
                          m=200,alpha=0.95,guess=c(3.31,0.1),
                          file="none"){
mle <- awarmle(x=x,cutoff=cutoff,guess=guess,conc=TRUE)
cmle <- mle$conc
mle <- mle$theta
bmles <- rep(0,length=m)
bconc <- rep(0,length=m)
for(i in seq(along=bmles)){
 xsamp <- sample(x,size=length(x),replace=TRUE)
 tbs <- awarmle(x=xsamp,cutoff=cutoff,guess=mle,conc=TRUE)
 if(is.na(tbs$theta[1])){bmles[i] <- bmles[i-1]; bconc[i] <- bconc[i-1]}
 else{
  bmles[i] <- tbs$theta[1]
  bconc[i] <- tbs$conc
 }
}
#
if(!missing(file)){
 save(cmle,mle,bconc,bmles,file=file)
}
mle <- c(quantile(bconc,c(0.5*(1-alpha),(1-alpha),0.5,alpha,0.5+0.5*alpha),na.rm=TRUE),cmle,
  mean(bmles < 3.0,na.rm=TRUE))
names(mle)[7] <- "prob < 3"
mle
}
#
# Bootstrap CI for Discrete Pareto Conc
#
bootstrapdpconc <- function(x,cutoff=1,cutabove=1000,
                          m=200,alpha=0.95,guess=3.31,
                          file="none"){
mle <- adpmle(x=x,cutoff=cutoff,guess=guess)
cmle <- mle$conc
mle <- mle$theta
bmles <- rep(0,length=m)
bconc <- rep(0,length=m)
for(i in seq(along=bmles)){
 xsamp <- sample(x,size=length(x),replace=TRUE)
 tbs <- adpmle(x=xsamp,cutoff=cutoff,guess=mle)
 if(is.na(tbs$theta[1])){bmles[i] <- bmles[i-1]; bconc[i] <- bconc[i-1]}
 else{
  bmles[i] <- tbs$theta
  bconc[i] <- tbs$conc
 }
}
#
if(!missing(file)){
 save(cmle,mle,bconc,bmles,file=file)
}
mle <- c(quantile(bconc,
  c(0.5*(1-alpha),(1-alpha),0.5,alpha,0.5+0.5*alpha),na.rm=TRUE),cmle,
  mean(bmles < 3.0,na.rm=TRUE))
names(mle)[7] <- "prob < 3"
mle
}
#
# Bootstrap CI for Yule Conc
#
bootstrapyuleconc <- function(x,cutoff=1,cutabove=1000,
                          m=200,alpha=0.95,guess=3.31,
                          file="none"){
mle <- ayulemle(x=x,cutoff=cutoff,guess=guess,conc=TRUE)
cmle <- mle$conc
mle <- mle$theta
bmles <- rep(0,length=m)
bconc <- rep(0,length=m)
for(i in seq(along=bmles)){
 xsamp <- sample(x,size=length(x),replace=TRUE)
 tbs <- ayulemle(x=xsamp,cutoff=cutoff,guess=mle,conc=TRUE)
 if(is.na(tbs$theta[1])){bmles[i] <- bmles[i-1]; bconc[i] <- bconc[i-1]}
 else{
  bmles[i] <- tbs$theta
  bconc[i] <- tbs$conc
 }
}
#
if(!missing(file)){
 save(cmle,mle,bconc,bmles,file=file)
}
mle <- c(quantile(bconc,
  c(0.5*(1-alpha),(1-alpha),0.5,alpha,0.5+0.5*alpha),na.rm=TRUE),cmle,
  mean(bmles < 3.0,na.rm=TRUE))
names(mle)[7] <- "prob < 3"
mle
}
#
# Bootstrap CI for Negative Binomial
#
bootstrapnb <- function(x,cutoff=1,cutabove=1000,
                          m=200,alpha=0.95,guess=c(5,0.2),
                          file="none"){
if(missing(guess)){
  mle <- anbmle(x=x,cutoff=cutoff)
  guess <- mle$theta
}else{
  mle <- anbmle(x=x,cutoff=cutoff,guess=guess)
}
bmles <- matrix(0,nrow=m,ncol=2)
for(i in seq(along=bmles[,1])){
 xsamp <- sample(x,size=length(x),replace=TRUE)
 aaa <- anbmle(x=xsamp,cutoff=cutoff,guess=guess)
 bmles[i,] <- aaa$theta
}
#
if(!missing(file)){
 save(mle$theta,bmles,file=file)
}
c(quantile(bmles[,1],c(0.5*(1-alpha),(1-alpha),0.5,alpha,0.5+0.5*alpha),na.rm=TRUE),
  quantile(bmles[,2],c(0.5*(1-alpha),(1-alpha),0.5,alpha,0.5+0.5*alpha),na.rm=TRUE),mle$theta)
}
#
# Bootstrap CI for Negative Binomial Concentration
# Used to be bootstrapnb
#
bootstrapnbconc <- function(x,cutoff=1,cutabove=1000,
                          m=200,alpha=0.95,guess=c(5,0.2),
                          file="none"){
if(missing(guess)){
  mle <- anbmle(x=x,cutoff=cutoff,conc=TRUE)
  guess <- mle$theta
}else{
  mle <- anbmle(x=x,cutoff=cutoff,guess=guess,conc=TRUE)
}
bmles <- rep(0,m)
for(i in seq(along=bmles)){
 xsamp <- sample(x,size=length(x),replace=TRUE)
 aaa <- anbmle(x=xsamp,cutoff=cutoff,guess=guess,conc=TRUE)
 bmles[i] <- aaa$conc
}
#
if(!missing(file)){
 save(mle$theta,bmles,file=file)
}
c(quantile(bmles,c(0.5*(1-alpha),(1-alpha),0.5,alpha,0.5+0.5*alpha),na.rm=TRUE),mle$theta)
}
#
# Bootstrap CI for Negative Binomial Yule Concentration
#
bootstrapnbyconc <- function(x,cutoff=1,cutabove=1000,
                          m=200,alpha=0.95,guess=c(3.5,50,0.1),
                          file="none"){
if(missing(guess)){
  mle <- anbymle(x=x,cutoff=cutoff,conc=TRUE)
  guess <- mle$theta
}else{
  mle <- anbymle(x=x,cutoff=cutoff,guess=guess,conc=TRUE)
}
cmle <- mle$conc
mle <- mle$theta
if(missing(guess)){guess <- mle$theta}
bconc <- rep(0,length=m)
for(i in seq(along=bconc)){
 xsamp <- sample(x,size=length(x),replace=TRUE)
 while(max(xsamp)<cutoff){
  xsamp <- sample(x,size=length(x),replace=TRUE)
 }
 tmles <- try(anbymle(x=xsamp,cutoff=cutoff,guess=mle,conc=TRUE)$conc)
 while(inherits(tmles,"try-error")){
  xsamp <- sample(x,size=length(x),replace=TRUE)
  while(max(xsamp)<cutoff){
   xsamp <- sample(x,size=length(x),replace=TRUE)
  }
  tmles <- try(anbymle(x=xsamp,cutoff=cutoff,guess=mle,conc=TRUE)$conc)
 }
 bconc[i] <- tmles
}
#
if(!missing(file)){
 save(cmle,mle,bconc,file=file)
}
c(quantile(bconc,c(0.5*(1-alpha),(1-alpha),0.5,alpha,0.5+0.5*alpha),na.rm=TRUE),cmle)
}
#
# Bootstrap CI for PowerExponential
#
#bootstrappe <- function(x,cutoff=1,cutabove=1000,
#                          m=200,alpha=0.95,xr=1:10000){
#aaa <- apemle(x=x,cutoff=cutoff,xr=xr)
#bmles <- matrix(0,ncol=2,nrow=m)
#for(i in seq(along=bmles[,1])){
# xsamp <- sample(x,size=length(x),replace=TRUE)
# bmles[i,] <- apemle(x=xsamp,cutoff=cutoff,xr=xr)$theta
#}
##
#bbb <- list(c(quantile(bmles[,1],c(0.5*(1-alpha),0.5,0.5+0.5*alpha),na.rm=TRUE),
#  quantile(bmles[,2],c(0.5*(1-alpha),0.5,0.5+0.5*alpha),na.rm=TRUE)),
#  cor=cor(bmles),
#  mle=aaa$theta)
##
#bbb
#}
bootstrappe <- function(x,cutoff=1,cutabove=1000,
                          m=200,alpha=0.95,guess=c(3.31,0.1),
                          file="none", conc=FALSE){
mle <- apemle(x=x,cutoff=cutoff,guess=guess)$theta
bmles <- rep(0,length=m)
for(i in seq(along=bmles)){
 xsamp <- sample(x,size=length(x),replace=TRUE)
 aaa <- apemle(x=xsamp,cutoff=cutoff,guess=guess)
 if(conc){
  bmles[i] <- aaa$conc
 }else{
  bmles[i] <- aaa$theta[1]
 }
}
#
if(!missing(file)){
 save(mle,bmles,file=file)
}
c(quantile(bmles,c(0.5*(1-alpha),(1-alpha),0.5,alpha,0.5+0.5*alpha),na.rm=TRUE),mle)
}
#
# Bootstrap CI for Waring Conc
#
bootstrappeconc <- function(x,cutoff=1,cutabove=1000,
                          m=200,alpha=0.95,guess=c(3.31,0.1),
                          file="none"){
mle <- apemle(x=x,cutoff=cutoff,guess=guess,conc=TRUE)
cmle <- mle$conc
mle <- mle$theta
bmles <- rep(0,length=m)
bconc <- rep(0,length=m)
for(i in seq(along=bmles)){
 xsamp <- sample(x,size=length(x),replace=TRUE)
 tbs <- apemle(x=xsamp,cutoff=cutoff,guess=mle,conc=TRUE)
 if(is.na(tbs$theta[1])){bmles[i] <- bmles[i-1]; bconc[i] <- bconc[i-1]}
 else{
  bmles[i] <- tbs$theta[1]
  bconc[i] <- tbs$conc
 }
}
#
if(!missing(file)){
 save(cmle,mle,bconc,bmles,file=file)
}
mle <- c(quantile(bconc,c(0.5*(1-alpha),(1-alpha),0.5,alpha,0.5+0.5*alpha),na.rm=TRUE),cmle,
  mean(bmles < 3.0,na.rm=TRUE))
names(mle)[7] <- "prob < 3"
mle
}
##
## Bootstrap CI for Power - Poisson Mixture
##
#bootstrapmix <- function(x,cutoff=1,
#                          m=200,alpha=0.95,xr=1:10000){
#aaa <- amixmle(x=x,cutoff=cutoff,xr=xr)
#bmles <- matrix(0,ncol=3,nrow=m)
#for(i in seq(along=bmles[,1])){
# xsamp <- sample(x,size=length(x),replace=TRUE)
# bmles[i,] <- amixmle(x=xsamp,cutoff=cutoff,xr=xr)$theta
#}
##
#bbb <- list(c(quantile(bmles[,1],c(0.5*(1-alpha),0.5,0.5+0.5*alpha),na.rm=TRUE),
#  quantile(bmles[,2],c(0.5*(1-alpha),0.5,0.5+0.5*alpha),na.rm=TRUE),
#  quantile(bmles[,3],c(0.5*(1-alpha),0.5,0.5+0.5*alpha),na.rm=TRUE)),
#  cor=cor(bmles),
#  mle=aaa$theta)
##
#bbb
#}
#
#bootstrapcor <- function(x,m=200){
#bmles <- matrix(0,nrow=m,ncol=max(x)+1)
#for(i in 1:m){
# xsamp <- sample(x,size=length(x),replace=TRUE)
# for(j in (0:max(x))){
#  bmles[i,j+1] <- mean(xsamp>=j)
# }
#}
##
#bmles
#}

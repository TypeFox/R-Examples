boxcoxnc <-
function(data,method="all",lam=seq(-2,2,0.01),plotit = TRUE,rep=30,p.method="BY")
{


data<-as.numeric(data)

if (is.na(min(data))==TRUE) stop("Data include NA")
if (min(data)<=0) stop("Data must include positive values")


if (method=="all") {


# Shapiro Wilk - in nortest
sw<-NULL
for (i in 1:length(lam)){
if (round(lam[i],2)!=0) sw<-rbind(sw,c(lam[i],shapiro.test((data**(lam[i])-1)/(lam[i]))$statistic))
if (round(lam[i],2)==0) sw<-rbind(sw,c(lam[i],shapiro.test(log(data))$statistic))
}
swlam<-sw[which.max(sw[,2]),1]

if (swlam!=0) data.sw<-((data^swlam)-1)/swlam
if (swlam==0) data.sw<-log(data)

sw.pvalue_data.sw<-p.adjust(c(shapiro.test(data.sw)$p.value,sf.test(data.sw)$p.value,jarque.bera.test(data.sw)$p.value),p.method)[1]
sf.pvalue_data.sw<-p.adjust(c(shapiro.test(data.sw)$p.value,sf.test(data.sw)$p.value,jarque.bera.test(data.sw)$p.value),p.method)[2]
jb.pvalue_data.sw<-p.adjust(c(shapiro.test(data.sw)$p.value,sf.test(data.sw)$p.value,jarque.bera.test(data.sw)$p.value),p.method)[3]



# Anderson Darling - in nortest
ad<-NULL
for (i in 1:length(lam)){
if (round(lam[i],2)!=0) ad<-rbind(ad,c(lam[i],ad.test((data**(lam[i])-1)/(lam[i]))$statistic))
if (round(lam[i],2)==0) ad<-rbind(ad,c(lam[i],ad.test(log(data))$statistic))
}
adlam<-ad[which.min(ad[,2]),1]

if (adlam!=0) data.ad<-((data^adlam)-1)/adlam
if (adlam==0) data.ad<-log(data)


sw.pvalue_data.ad<-p.adjust(c(shapiro.test(data.ad)$p.value,sf.test(data.ad)$p.value,jarque.bera.test(data.ad)$p.value),p.method)[1]
sf.pvalue_data.ad<-p.adjust(c(shapiro.test(data.ad)$p.value,sf.test(data.ad)$p.value,jarque.bera.test(data.ad)$p.value),p.method)[2]
jb.pvalue_data.ad<-p.adjust(c(shapiro.test(data.ad)$p.value,sf.test(data.ad)$p.value,jarque.bera.test(data.ad)$p.value),p.method)[3]


# Cramer Von-Mises - in nortest
cvm<-NULL
for (i in 1:length(lam)){
if (round(lam[i],2)!=0) cvm<-rbind(cvm,c(lam[i],cvm.test((data**(lam[i])-1)/(lam[i]))$statistic))
if (round(lam[i],2)==0) cvm<-rbind(cvm,c(lam[i],cvm.test(log(data))$statistic))
}
cvmlam<-cvm[which.min(cvm[,2]),1]


if (cvmlam!=0) data.cvm<-((data^cvmlam)-1)/cvmlam
if (cvmlam==0) data.cvm<-log(data)

sw.pvalue_data.cvm<-p.adjust(c(shapiro.test(data.cvm)$p.value,sf.test(data.cvm)$p.value,jarque.bera.test(data.cvm)$p.value),p.method)[1]
sf.pvalue_data.cvm<-p.adjust(c(shapiro.test(data.cvm)$p.value,sf.test(data.cvm)$p.value,jarque.bera.test(data.cvm)$p.value),p.method)[2]
jb.pvalue_data.cvm<-p.adjust(c(shapiro.test(data.cvm)$p.value,sf.test(data.cvm)$p.value,jarque.bera.test(data.cvm)$p.value),p.method)[3]



# Pearson Test - in nortest
pt<-NULL
for (i in 1:length(lam)){
if (round(lam[i],2)!=0) pt<-rbind(pt,c(lam[i],pearson.test((data**(lam[i])-1)/(lam[i]))$statistic))
if (round(lam[i],2)==0) pt<-rbind(pt,c(lam[i],pearson.test(log(data))$statistic))
}
ptlam<-pt[which.min(pt[,2]),1]


if (ptlam!=0) data.pt<-((data^ptlam)-1)/ptlam
if (ptlam==0) data.pt<-log(data)

sw.pvalue_data.pt<-p.adjust(c(shapiro.test(data.pt)$p.value,sf.test(data.pt)$p.value,jarque.bera.test(data.pt)$p.value),p.method)[1]
sf.pvalue_data.pt<-p.adjust(c(shapiro.test(data.pt)$p.value,sf.test(data.pt)$p.value,jarque.bera.test(data.pt)$p.value),p.method)[2]
jb.pvalue_data.pt<-p.adjust(c(shapiro.test(data.pt)$p.value,sf.test(data.pt)$p.value,jarque.bera.test(data.pt)$p.value),p.method)[3]



#Shapiro Francia - in nortest
sf<-NULL
for (i in 1:length(lam)){
if (round(lam[i],2)!=0) sf<-rbind(sf,c(lam[i],sf.test((data**(lam[i])-1)/(lam[i]))$statistic))
if (round(lam[i],2)==0) sf<-rbind(sf,c(lam[i],sf.test(log(data))$statistic))
}
sflam<-sf[which.max(sf[,2]),1]


if (sflam!=0) data.sf<-((data^sflam)-1)/sflam
if (sflam==0) data.sf<-log(data)

sw.pvalue_data.sf<-p.adjust(c(shapiro.test(data.sf)$p.value,sf.test(data.sf)$p.value,jarque.bera.test(data.sf)$p.value),p.method)[1]
sf.pvalue_data.sf<-p.adjust(c(shapiro.test(data.sf)$p.value,sf.test(data.sf)$p.value,jarque.bera.test(data.sf)$p.value),p.method)[2]
jb.pvalue_data.sf<-p.adjust(c(shapiro.test(data.sf)$p.value,sf.test(data.sf)$p.value,jarque.bera.test(data.sf)$p.value),p.method)[3]



# Lilliefors - Kolmogorov Smirnov - in nortest
lt<-NULL
for (i in 1:length(lam)){
if (round(lam[i],2)!=0) lt<-rbind(lt,c(lam[i],lillie.test((data**(lam[i])-1)/(lam[i]))$statistic))
if (round(lam[i],2)==0) lt<-rbind(lt,c(lam[i],lillie.test(log(data))$statistic))
}
ltlam<-lt[which.min(lt[,2]),1]


if (ltlam!=0) data.lt<-((data^ltlam)-1)/ltlam
if (ltlam==0) data.lt<-log(data)

sw.pvalue_data.lt<-p.adjust(c(shapiro.test(data.lt)$p.value,sf.test(data.lt)$p.value,jarque.bera.test(data.lt)$p.value),p.method)[1]
sf.pvalue_data.lt<-p.adjust(c(shapiro.test(data.lt)$p.value,sf.test(data.lt)$p.value,jarque.bera.test(data.lt)$p.value),p.method)[2]
jb.pvalue_data.lt<-p.adjust(c(shapiro.test(data.lt)$p.value,sf.test(data.lt)$p.value,jarque.bera.test(data.lt)$p.value),p.method)[3]



# Jarque-Bera in tseries
jb<-NULL
for (i in 1:length(lam)){
if (round(lam[i],2)!=0) jb<-rbind(jb,c(lam[i],jarque.bera.test((data**(lam[i])-1)/(lam[i]))$statistic))
if (round(lam[i],2)==0) jb<-rbind(jb,c(lam[i],jarque.bera.test(log(data))$statistic))
}
jblam<-jb[which.min(jb[,2]),1]


if (jblam!=0) data.jb<-((data^jblam)-1)/jblam
if (jblam==0) data.jb<-log(data)

sw.pvalue_data.jb<-p.adjust(c(shapiro.test(data.jb)$p.value,sf.test(data.jb)$p.value,jarque.bera.test(data.jb)$p.value),p.method)[1]
sf.pvalue_data.jb<-p.adjust(c(shapiro.test(data.jb)$p.value,sf.test(data.jb)$p.value,jarque.bera.test(data.jb)$p.value),p.method)[2]
jb.pvalue_data.jb<-p.adjust(c(shapiro.test(data.jb)$p.value,sf.test(data.jb)$p.value,jarque.bera.test(data.jb)$p.value),p.method)[3]




aclam1<-NULL
for (q in 1:rep) {
ac<-rnorm(length(data),0,100)
lm1<-glm(data~ac,family=gaussian)

bc1<-boxcox(lm1,lam,plotit=FALSE)
aclam<-bc1$x[which.max(bc1$y)]
aclam1<-cbind(aclam1,aclam)
}



if (plotit==TRUE){
par(mfrow=c(2,4))
plot(sw[,1],sw[,2],ylab="test statistic",xlab=expression(lambda), main="Shapiro-Wilk")
abline(v=swlam,lty=2)
plot(ad[,1],ad[,2],ylab="test statistic",xlab=expression(lambda),main="Anderson-Darling")
abline(v=adlam,lty=2)
plot(cvm[,1],cvm[,2],ylab="test statistic",xlab=expression(lambda),main="Cramer-von Mises")
abline(v=cvmlam,lty=2)
plot(pt[,1],pt[,2],ylab="test statistic",xlab=expression(lambda),main="Pearson Chi-Square")
abline(v=ptlam,lty=2)
plot(sf[,1],sf[,2],ylab="test statistic",xlab=expression(lambda),main="Shapiro-Francia")
abline(v=sflam,lty=2)
plot(lt[,1],lt[,2],ylab="test statistic",xlab=expression(lambda),main="Lilliefors")
abline(v=ltlam,lty=2)
plot(jb[,1],jb[,2],ylab="test statistic",xlab=expression(lambda),main="Jarque-Bera")
abline(v=jblam,lty=2)
boxcox(lm1,lam,plotit)

if (plotit==TRUE){
title("Artificial Covariate")
}

}



aclam1<-as.numeric(aclam1)
aclam2<-mean(aclam1)

if (swlam==max(lam)) stop("Enlarge the range of the lambda in a positive direction")
if (swlam==min(lam)) stop("Enlarge the range of the lambda in a negative direction")
if (adlam==max(lam)) stop("Enlarge the range of the lambda in a positive direction")
if (adlam==min(lam)) stop("Enlarge the range of the lambda in a negative direction")
if (cvmlam==max(lam)) stop("Enlarge the range of the lambda in a positive direction")
if (cvmlam==min(lam)) stop("Enlarge the range of the lambda in a negative direction")
if (ptlam==max(lam)) stop("Enlarge the range of the lambda in a positive direction")
if (ptlam==min(lam)) stop("Enlarge the range of the lambda in a negative direction")
if (sflam==max(lam)) stop("Enlarge the range of the lambda in a positive direction")
if (sflam==min(lam)) stop("Enlarge the range of the lambda in a negative direction")
if (ltlam==max(lam)) stop("Enlarge the range of the lambda in a positive direction")
if (ltlam==min(lam)) stop("Enlarge the range of the lambda in a negative direction")
if (jblam==max(lam)) stop("Enlarge the range of the lambda in a positive direction")
if (jblam==min(lam)) stop("Enlarge the range of the lambda in a negative direction")
if (aclam==max(lam)) stop("Enlarge the range of the lambda in a positive direction")
if (aclam==min(lam)) stop("Enlarge the range of the lambda in a negative direction")





if (aclam2!=0) data.ac<-((data^aclam2)-1)/aclam2
if (aclam2==0) data.ac<-log(data)

sw.pvalue_data.ac<-p.adjust(c(shapiro.test(data.ac)$p.value,sf.test(data.ac)$p.value,jarque.bera.test(data.ac)$p.value),p.method)[1]
sf.pvalue_data.ac<-p.adjust(c(shapiro.test(data.ac)$p.value,sf.test(data.ac)$p.value,jarque.bera.test(data.ac)$p.value),p.method)[2]
jb.pvalue_data.ac<-p.adjust(c(shapiro.test(data.ac)$p.value,sf.test(data.ac)$p.value,jarque.bera.test(data.ac)$p.value),p.method)[3]




result<-matrix(0,4,8)
colnames(result)<-c("sw","ad","cvm","pt","sf","lt","jb","ac")
rownames(result)<-c("lambda.hat","sw.pvalue","sf.pvalue","jb.pvalue")

result[1,]<-c(swlam,adlam,cvmlam,ptlam,sflam,ltlam,jblam,aclam2)
result[2,]<-c(sw.pvalue_data.sw,sw.pvalue_data.ad,sw.pvalue_data.cvm,sw.pvalue_data.pt,sw.pvalue_data.sf,sw.pvalue_data.lt,sw.pvalue_data.jb,sw.pvalue_data.ac)
result[3,]<-c(sf.pvalue_data.sw,sf.pvalue_data.ad,sf.pvalue_data.cvm,sf.pvalue_data.pt,sf.pvalue_data.sf,sf.pvalue_data.lt,sf.pvalue_data.jb,sf.pvalue_data.ac)
result[4,]<-c(jb.pvalue_data.sw,jb.pvalue_data.ad,jb.pvalue_data.cvm,jb.pvalue_data.pt,jb.pvalue_data.sf,jb.pvalue_data.lt,jb.pvalue_data.jb,jb.pvalue_data.ac)


out<-list()
out$title<-"Implementation of Box-Cox Power Transformation when No Covariate Is Available"
out$method="All"
out$date<-date()
out$result<-result
out
}

else if (method=="sw") {
sw<-NULL
for (i in 1:length(lam)){
if (round(lam[i],2)!=0) sw<-rbind(sw,c(lam[i],shapiro.test((data**(lam[i])-1)/(lam[i]))$statistic))
if (round(lam[i],2)==0) sw<-rbind(sw,c(lam[i],shapiro.test(log(data))$statistic))
}

swlam<-sw[which.max(sw[,2]),1]

if (plotit==TRUE){
plot(sw[,1],sw[,2],ylab="test statistic",xlab=expression(lambda), main="Shapiro-Wilk")
abline(v=swlam,lty=2)
}


if (swlam==max(lam)) stop("Enlarge the range of the lambda in a positive direction")
if (swlam==min(lam)) stop("Enlarge the range of the lambda in a negative direction")



if (swlam!=0) data.sw<-((data^swlam)-1)/swlam
if (swlam==0) data.sw<-log(data)


sw.pvalue_data.sw<-p.adjust(c(shapiro.test(data.sw)$p.value,sf.test(data.sw)$p.value,jarque.bera.test(data.sw)$p.value),p.method)[1]
sf.pvalue_data.sw<-p.adjust(c(shapiro.test(data.sw)$p.value,sf.test(data.sw)$p.value,jarque.bera.test(data.sw)$p.value),p.method)[2]
jb.pvalue_data.sw<-p.adjust(c(shapiro.test(data.sw)$p.value,sf.test(data.sw)$p.value,jarque.bera.test(data.sw)$p.value),p.method)[3]


result<-matrix(0,4,1)
colnames(result)<-c("sw");rownames(result)<-c("lambda.hat","sw.pvalue","sf.pvalue","jb.pvalue")
result[,1]<-c(swlam,sw.pvalue_data.sw,sf.pvalue_data.sw,jb.pvalue_data.sw)

out<-list()
out$title<-"Implementation of Box-Cox Power Transformation when No Covariate Is Available"
out$method="Shapiro-Wilk"
out$date<-date()
out$result<-result
out
}

else if (method=="ad") {

ad<-NULL
for (i in 1:length(lam)){
if (round(lam[i],2)!=0) ad<-rbind(ad,c(lam[i],ad.test((data**(lam[i])-1)/(lam[i]))$statistic))
if (round(lam[i],2)==0) ad<-rbind(ad,c(lam[i],ad.test(log(data))$statistic))
}

adlam<-ad[which.min(ad[,2]),1]
adstat<-ad[which.min(ad[,2]),2]

if (plotit==TRUE){
plot(ad[,1],ad[,2],ylab="test statistic",xlab=expression(lambda),main="Anderson-Darling")
abline(v=adlam,lty=2)
}


if (adlam==max(lam)) stop("Enlarge the range of the lambda in a positive direction")
if (adlam==min(lam)) stop("Enlarge the range of the lambda in a negative direction")


if (adlam!=0) data.ad<-((data^adlam)-1)/adlam
if (adlam==0) data.ad<-log(data)

sw.pvalue_data.ad<-p.adjust(c(shapiro.test(data.ad)$p.value,sf.test(data.ad)$p.value,jarque.bera.test(data.ad)$p.value),p.method)[1]
sf.pvalue_data.ad<-p.adjust(c(shapiro.test(data.ad)$p.value,sf.test(data.ad)$p.value,jarque.bera.test(data.ad)$p.value),p.method)[2]
jb.pvalue_data.ad<-p.adjust(c(shapiro.test(data.ad)$p.value,sf.test(data.ad)$p.value,jarque.bera.test(data.ad)$p.value),p.method)[3]


result<-matrix(0,4,1)
colnames(result)<-c("ad");rownames(result)<-c("lambda.hat","sw.pvalue","sf.pvalue","jb.pvalue")
result[,1]<-c(adlam,sw.pvalue_data.ad,sf.pvalue_data.ad,jb.pvalue_data.ad)

out<-list()
out$title<-"Implementation of Box-Cox Power Transformation when No Covariate Is Available"
out$method="Anderson-Darling"
out$date<-date()
out$result<-result
out
}

else if (method=="cvm") {

# Cramer Von-Mises - in nortest
cvm<-NULL
for (i in 1:length(lam)){
if (round(lam[i],2)!=0) cvm<-rbind(cvm,c(lam[i],cvm.test((data**(lam[i])-1)/(lam[i]))$statistic))
if (round(lam[i],2)==0) cvm<-rbind(cvm,c(lam[i],cvm.test(log(data))$statistic))
}

cvmlam<-cvm[which.min(cvm[,2]),1]
cvmstat<-cvm[which.min(cvm[,2]),2]


if (plotit==TRUE){
plot(cvm[,1],cvm[,2],ylab="test statistic",xlab=expression(lambda),main="Cramer-von Mises")
abline(v=cvmlam,lty=2)
}


if (cvmlam==max(lam)) stop("Enlarge the range of the lambda in a positive direction")
if (cvmlam==min(lam)) stop("Enlarge the range of the lambda in a negative direction")


if (cvmlam!=0) data.cvm<-((data^cvmlam)-1)/cvmlam
if (cvmlam==0) data.cvm<-log(data)

sw.pvalue_data.cvm<-p.adjust(c(shapiro.test(data.cvm)$p.value,sf.test(data.cvm)$p.value,jarque.bera.test(data.cvm)$p.value),p.method)[1]
sf.pvalue_data.cvm<-p.adjust(c(shapiro.test(data.cvm)$p.value,sf.test(data.cvm)$p.value,jarque.bera.test(data.cvm)$p.value),p.method)[2]
jb.pvalue_data.cvm<-p.adjust(c(shapiro.test(data.cvm)$p.value,sf.test(data.cvm)$p.value,jarque.bera.test(data.cvm)$p.value),p.method)[3]

result<-matrix(0,4,1)
colnames(result)<-c("cvm");rownames(result)<-c("lambda.hat","sw.pvalue","sf.pvalue","jb.pvalue")
result[,1]<-c(cvmlam,sw.pvalue_data.cvm,sf.pvalue_data.cvm,jb.pvalue_data.cvm)


out<-list()
out$title<-"Implementation of Box-Cox Power Transformation when No Covariate Is Available"
out$method="Cramer-von Mises"
out$date<-date()
out$result<-result
out
}

else if (method=="pt") {

pt<-NULL
for (i in 1:length(lam)){
if (round(lam[i],2)!=0) pt<-rbind(pt,c(lam[i],pearson.test((data**(lam[i])-1)/(lam[i]))$statistic))
if (round(lam[i],2)==0) pt<-rbind(pt,c(lam[i],pearson.test(log(data))$statistic))
}

ptlam<-pt[which.min(pt[,2]),1]
ptstat<-pt[which.min(pt[,2]),2]

if (plotit==TRUE){
plot(pt[,1],pt[,2],ylab="test statistic",xlab=expression(lambda),main="Pearson Chi-Square")
abline(v=ptlam,lty=2)
}


if (ptlam==max(lam)) stop("Enlarge the range of the lambda in a positive direction")
if (ptlam==min(lam)) stop("Enlarge the range of the lambda in a negative direction")


if (ptlam!=0) data.pt<-((data^ptlam)-1)/ptlam
if (ptlam==0) data.pt<-log(data)

sw.pvalue_data.pt<-p.adjust(c(shapiro.test(data.pt)$p.value,sf.test(data.pt)$p.value,jarque.bera.test(data.pt)$p.value),p.method)[1]
sf.pvalue_data.pt<-p.adjust(c(shapiro.test(data.pt)$p.value,sf.test(data.pt)$p.value,jarque.bera.test(data.pt)$p.value),p.method)[2]
jb.pvalue_data.pt<-p.adjust(c(shapiro.test(data.pt)$p.value,sf.test(data.pt)$p.value,jarque.bera.test(data.pt)$p.value),p.method)[3]

result<-matrix(0,4,1)
colnames(result)<-c("pt");rownames(result)<-c("lambda.hat","sw.pvalue","sf.pvalue","jb.pvalue")
result[,1]<-c(ptlam,sw.pvalue_data.pt,sf.pvalue_data.pt,jb.pvalue_data.pt)



out<-list()
out$title<-"Implementation of Box-Cox Power Transformation when No Covariate Is Available"
out$method="Pearson Chi-Square"
out$date<-date()
out$result<-result
out
}

else if (method=="sf") {

sf<-NULL
for (i in 1:length(lam)){
if (round(lam[i],2)!=0) sf<-rbind(sf,c(lam[i],sf.test((data**(lam[i])-1)/(lam[i]))$statistic))
if (round(lam[i],2)==0) sf<-rbind(sf,c(lam[i],sf.test(log(data))$statistic))
}

sflam<-sf[which.max(sf[,2]),1]
sfstat<-sf[which.max(sf[,2]),2]

if (plotit==TRUE){
plot(sf[,1],sf[,2],ylab="test statistic",xlab=expression(lambda),main="Shapiro-Francia")
abline(v=sflam,lty=2)
}


if (sflam==max(lam)) stop("Enlarge the range of the lambda in a positive direction")
if (sflam==min(lam)) stop("Enlarge the range of the lambda in a negative direction")


if (sflam!=0) data.sf<-((data^sflam)-1)/sflam
if (sflam==0) data.sf<-log(data)


sw.pvalue_data.sf<-p.adjust(c(shapiro.test(data.sf)$p.value,sf.test(data.sf)$p.value,jarque.bera.test(data.sf)$p.value),p.method)[1]
sf.pvalue_data.sf<-p.adjust(c(shapiro.test(data.sf)$p.value,sf.test(data.sf)$p.value,jarque.bera.test(data.sf)$p.value),p.method)[2]
jb.pvalue_data.sf<-p.adjust(c(shapiro.test(data.sf)$p.value,sf.test(data.sf)$p.value,jarque.bera.test(data.sf)$p.value),p.method)[3]

result<-matrix(0,4,1)
colnames(result)<-c("sf");rownames(result)<-c("lambda.hat","sw.pvalue","sf.pvalue","jb.pvalue")
result[,1]<-c(sflam,sw.pvalue_data.sf,sf.pvalue_data.sf,jb.pvalue_data.sf)


out<-list()
out$title<-"Implementation of Box-Cox Power Transformation when No Covariate Is Available"
out$method="Shapiro-Francia"
out$date<-date()
out$result<-result
out
}

else if (method=="lt") {

lt<-NULL
for (i in 1:length(lam)){
if (round(lam[i],2)!=0) lt<-rbind(lt,c(lam[i],lillie.test((data**(lam[i])-1)/(lam[i]))$statistic))
if (round(lam[i],2)==0) lt<-rbind(lt,c(lam[i],lillie.test(log(data))$statistic))
}

ltlam<-lt[which.min(lt[,2]),1]
ltstat<-lt[which.min(lt[,2]),2]

if (plotit==TRUE){
plot(lt[,1],lt[,2],ylab="test statistic",xlab=expression(lambda),main="Lilliefors")
abline(v=ltlam,lty=2)
}

if (ltlam==max(lam)) stop("Enlarge the range of the lambda in a positive direction")
if (ltlam==min(lam)) stop("Enlarge the range of the lambda in a negative direction")


if (ltlam!=0) data.lt<-((data^ltlam)-1)/ltlam
if (ltlam==0) data.lt<-log(data)

sw.pvalue_data.lt<-p.adjust(c(shapiro.test(data.lt)$p.value,sf.test(data.lt)$p.value,jarque.bera.test(data.lt)$p.value),p.method)[1]
sf.pvalue_data.lt<-p.adjust(c(shapiro.test(data.lt)$p.value,sf.test(data.lt)$p.value,jarque.bera.test(data.lt)$p.value),p.method)[2]
jb.pvalue_data.lt<-p.adjust(c(shapiro.test(data.lt)$p.value,sf.test(data.lt)$p.value,jarque.bera.test(data.lt)$p.value),p.method)[3]

result<-matrix(0,4,1)
colnames(result)<-c("lt");rownames(result)<-c("lambda.hat","sw.pvalue","sf.pvalue","jb.pvalue")
result[,1]<-c(ltlam,sw.pvalue_data.lt,sf.pvalue_data.lt,jb.pvalue_data.lt)


out<-list()
out$title<-"Implementation of Box-Cox Power Transformation when No Covariate Is Available"
out$method="Lilliefors"
out$date<-date()
out$result<-result
out
}

else if (method=="jb") {

jb<-NULL
for (i in 1:length(lam)){
if (round(lam[i],2)!=0) jb<-rbind(jb,c(lam[i],jarque.bera.test((data**(lam[i])-1)/(lam[i]))$statistic))
if (round(lam[i],2)==0) jb<-rbind(jb,c(lam[i],jarque.bera.test(log(data))$statistic))
}


jblam<-jb[which.min(jb[,2]),1]
jbstat<-jb[which.min(jb[,2]),2]

if (plotit==TRUE){
plot(jb[,1],jb[,2],ylab="test statistic",xlab=expression(lambda),main="Jarque-Bera")
abline(v=jblam,lty=2)
}

if (jblam==max(lam)) stop("Enlarge the range of the lambda in a positive direction")
if (jblam==min(lam)) stop("Enlarge the range of the lambda in a negative direction")


if (jblam!=0) data.jb<-((data^jblam)-1)/jblam
if (jblam==0) data.jb<-log(data)

sw.pvalue_data.jb<-p.adjust(c(shapiro.test(data.jb)$p.value,sf.test(data.jb)$p.value,jarque.bera.test(data.jb)$p.value),p.method)[1]
sf.pvalue_data.jb<-p.adjust(c(shapiro.test(data.jb)$p.value,sf.test(data.jb)$p.value,jarque.bera.test(data.jb)$p.value),p.method)[2]
jb.pvalue_data.jb<-p.adjust(c(shapiro.test(data.jb)$p.value,sf.test(data.jb)$p.value,jarque.bera.test(data.jb)$p.value),p.method)[3]

result<-matrix(0,4,1)
colnames(result)<-c("jb");rownames(result)<-c("lambda.hat","sw.pvalue","sf.pvalue","jb.pvalue")
result[,1]<-c(jblam,sw.pvalue_data.jb,sf.pvalue_data.jb,jb.pvalue_data.jb)


out<-list()
out$title<-"Implementation of Box-Cox Power Transformation when No Covariate Is Available"
out$method="Jarque-Bera"
out$date<-date()
out$result<-result
out
}

else if (method=="ac") {
aclam1<-NULL

for (t in 1:rep) {

ac<-rnorm(length(data),0,100)
lm1<-glm(data~ac,family=gaussian)
bc1<-boxcox(lm1,lam,plotit=FALSE)
aclam<-bc1$x[which.max(bc1$y)]
aclam1<-cbind(aclam1,aclam)
}

boxcox(lm1,lam,plotit)
if (plotit==TRUE){
title("Artificial Covariate")
}
aclam1<-as.numeric(aclam1)
aclam2<-mean(aclam1)

if (aclam2==max(lam)) stop("Enlarge the range of the lambda in a positive direction")
if (aclam2==min(lam)) stop("Enlarge the range of the lambda in a negative direction")


if (aclam2!=0) data.ac<-((data^aclam2)-1)/aclam2
if (aclam2==0) data.ac<-log(data)

sw.pvalue_data.ac<-p.adjust(c(shapiro.test(data.ac)$p.value,sf.test(data.ac)$p.value,jarque.bera.test(data.ac)$p.value),p.method)[1]
sf.pvalue_data.ac<-p.adjust(c(shapiro.test(data.ac)$p.value,sf.test(data.ac)$p.value,jarque.bera.test(data.ac)$p.value),p.method)[2]
jb.pvalue_data.ac<-p.adjust(c(shapiro.test(data.ac)$p.value,sf.test(data.ac)$p.value,jarque.bera.test(data.ac)$p.value),p.method)[3]

result<-matrix(0,4,1)
colnames(result)<-c("ac");rownames(result)<-c("lambda.hat","sw.pvalue","sf.pvalue","jb.pvalue")
result[,1]<-c(aclam2,sw.pvalue_data.ac,sf.pvalue_data.ac,jb.pvalue_data.ac)

out<-list()
out$title<-"Implementation of Box-Cox Power Transformation when No Covariate Is Available"
out$method="Artifical Covariate"
out$date<-date()
out$result<-result
out
}

}

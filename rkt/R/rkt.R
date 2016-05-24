rkt <- function(date,y,block,cv,correct=F,rep="e"){

#
# preparing output
#

ans<-list(sl=NA,S=NA,B=NA,varS=NA,sl.corrected=NA,varS.corrected=NA,partial.S=NA,partial.sl=NA,partial.varS=NA,partial.sl.corrected=NA,partial.varS.corrected=NA,tau=NA)
oldClass(ans)<-"rkt"

#
# testing the input parameters
#

if (missing(date)|missing(y)) {
cat ("Error: first mandatory vector should be positive: year or year+fraction\n")
cat ("       second mandatory vector should be numerical: measured data (ties allowed)\n")
cat ("       third optional vector should be positive integer: season, month, site or a unique code for season and site\n")
cat ("       fourth optional vector should be numerical: covariable\n") 
return(ans)
}

if (missing(block)) block<-replicate(length(date),1)
a<-search()
ok=0
if(!is.numeric(date)|!is.numeric(block)|!is.numeric(y)|length(date)<4|length(y)!=length(block)|length(y)!=length(date)|min(date)<0 ) {
cat ("Error: first mandatory vector should be positive: year or year+fraction\n")
cat ("       second mandatory vector should be numerical: measured data (ties allowed)\n")
cat ("       third optional vector should be positive integer: season, month, site or a unique code for season and site\n")
cat ("       fourth optional vector should be numerical: covariable\n") 
cat ("       all with the same length, at least 4\n") 
return(ans)
 }

if (!missing(cv)) {
if (!is.numeric(cv)|length(cv)!=length(date)) {
cat ("Error: fourth vector should be numerical: covariable\n") 
cat ("       with the same length as the others\n") 
return(ans)
}
}
block<-trunc(block)
if (min(block)<1) {
cat ("Error: only positive values accepted for blocks\n")
return(ans)
 }
block<-trunc(block)

#
# if correction for correlation between block is requested, only one datum per year per block is allowed
#

if (correct) date<-trunc(date)

#
# replacing values sharing the same date with median or average
#

if (rep=="a"){
for (i in 1:length(date)){
date_ties<-y[date==date[i] & block==block[i] & !is.na(y)]
if (!missing(cv)){
cv_ties<-cv[date==date[i] & block==block[i] & !is.na(y)]
}
if (length(date_ties)>1) {
y[date==date[i]& block==block[i]]<-NA
y[i]<-mean(date_ties[!is.na(date_ties)])
if (!missing(cv)){
cv[date==date[i]& block==block[i]]<-NA
cv[i]<-mean(cv_ties[!is.na(date_ties)])
}
}
}
} 
if (rep=="m"){
for (i in 1:length(date)){
date_ties<-y[date==date[i] & block==block[i] & !is.na(y)]
if (!missing(cv)){
cv_ties<-cv[date==date[i] & block==block[i] & !is.na(y)]
}
if (length(date_ties)>1) {
y[date==date[i]& block==block[i]]<-NA
y[i]<-median(date_ties[!is.na(date_ties)])
if (!missing(cv)){
cv[date==date[i]& block==block[i]]<-NA
cv[i]<-median(cv_ties[!is.na(date_ties)])
}
}
}
}

#
# eliminating records with missing values for y 
# from date,block,y and cv 
#

date=date[!is.na(y)]
block=block[!is.na(y)]
if (!missing(cv)) cv=cv[!is.na(y)]
y=y[!is.na(y)]

#
# making blocks consecutive integers from 1 to maxblock
#

for (i in max(block):1){
if (length(block[block==i])==0){
block[block>=i]<-block[block>=i]-1
}
}

#
# initializing variables
#

tau<-NA
maxblock<-max(block)
minyear=min(date)
ny<-max(date)-minyear+1

S<-0
Scv<-0
varS<-0

pc<-NA
varSc<-NA

p.S<-NA
p.varS<-NA
p.p<-NA

varScv<-0

p.pc<-NA
p.varSc<-NA

sen<-date[0]
varmix<-0

#
# SKT/RKT and Sen on y (see Hirsch et al. 1982)
#

kd=split(date,block)
ky=split(y,block)
for (b in 1:maxblock) {
kkd=kd[[b]]
kky=ky[[b]]
if (length(kkd)>3) {
for (i1 in 2:length(kkd)){
for (i2 in 1:(i1-1)) {
sen<-c(sen,(kky[i2]-kky[i1])/(kkd[i2]-kkd[i1]))
S<-S+sign1((kky[i2]-kky[i1])/(kkd[i2]-kkd[i1]))
}
 }
nd<-length(kky[!is.na(kky)])
tn<-length(kky[is.na(kky)])
SumTies<-0
for (i1 in 1:length(kky)) {
if (!is.na(kky[i1])) {
ties<-length(kky[kky==kky[i1]])-tn
SumTies<-SumTies+(ties-1)*(2*ties+5)
 }
}
varS<-varS+nd*(nd-1)*(2*nd+5)/18-SumTies/18
} else cat ("1 block with less than 4 points ignored\n")
}
B<-median(sen)
ncomp<-length(sen)

#
# two equal dates gives delta/0 = Inf or 0/0 = NaN
#

if (length(sen[sen==Inf])>0 |  length(sen[is.nan(sen)])>0) {
cat("At least two equal dates in the same block\n")
return(ans)
} 

#
# if all block were too short, sen is empty and SKT/RKT not performed
#

if (length(sen)==0) 
{
cat("Less than 4 dates in all blocks\n")
return(ans)
}

if (S==0)
{
p=1
}
else
{
zeta<-(S-sign(S))/sqrt(varS)
p<-(1-pnorm(abs(zeta)))*2
}

#
# performing SKT/RKT on covariable for partial test (see Liniseller and Grinvall 2002)
#

if (!missing(cv)) {
kcv=split(cv,block)
for (b in 1:maxblock) {
kkd=kd[[b]]
kkcv=kcv[[b]]
kky=ky[[b]]
if (length(kkd)>3) {
for (i1 in 2:length(kkd)){
for (i2 in 1:(i1-1)) {
Scv<-Scv+sign1((kkcv[i2]-kkcv[i1])/(kkd[i2]-kkd[i1]))
}
 }
nd<-length(kkcv[!is.na(kkcv)])
tn<-length(kkcv[is.na(kkcv)])
SumTies<-0
for (i1 in 1:length(kkcv)) {
if (!is.na(kkcv[i1])) {
ties<-length(kkcv[kkcv==kkcv[i1]])-tn
SumTies<-SumTies+(ties-1)*(2*ties+5)
 }
}
varScv<-varScv+nd*(nd-1)*(2*nd+5)/18-SumTies/18
} else cat ("1 block with less than 4 points ignored\n")
Ry<-rank(kky,na.last="keep",ties.method="average")
Ry[is.na(Ry)]<-mean(Ry[!is.na(Ry)])
Rcv<-rank(kkcv,na.last="keep",ties.method="average")
Rcv[is.na(Rcv)]<-mean(Rcv[!is.na(Rcv)])
Kmix<-0
for (i in 2:length(kky)) {
for (j in 1:(i-1)) Kmix<-Kmix+sign1((kkcv[i]-kkcv[j])*(kky[i]-kky[j]))
}
Gmix<-0
for (i in 1:length(kky)) Gmix<-Gmix+Ry[i]*Rcv[i]
varmix<-varmix+(Kmix+4*Gmix-length(kkd)*(length(kky[!is.na(kky)])+1)*(length(kkcv[!is.na(kkcv)])+1))/3
}
p.varS=varS-varmix*varmix/varScv
p.S<-S-varmix/varScv*Scv
if (p.S==0)
{
p.p=1
}
else
{
zeta<-(p.S-sign(p.S))/sqrt(p.varS)
p.p<-(1-pnorm(abs(zeta)))*2
}
}

#
# correction for correlation among blocks (see Hirsch & Slack 1984) for at least 10 years
#

if (correct  & abs(S)>0 & ny>9 & maxblock>1) {
X<-rep(NA,times=(ny*maxblock))
dim(X)<- c(ny,maxblock)
Rks<-X
varSc<-0
for (i in 1:length(date)) X[date[i]-minyear+1,block[i]]<-y[i]
for (b in 1:maxblock){ 
R<-X[,b]
R<-rank(R,na.last="keep",ties.method="average")
R[is.na(R)]<-mean(R[!is.na(R)])
Rks[,b]<-R
}
for (g in 1:maxblock) {
for (h in 1:maxblock){
K<-0
for (i in 2:ny) {
for (j in 1:(i-1)) K<-K+sign1((X[j,g]-X[i,g])*(X[j,h]-X[i,h]))
}
G<-0
for (i in 1:ny) G<-G+Rks[i,g]*Rks[i,h] 

xg<-X[,g]
xh<-X[,h]
varSc<-varSc+(K+4*G-length(xg)*(length(xg[!is.na(xg)])+1)*(length(xh[!is.na(xh)])+1))/3

}
}
zeta<-(S-sign(S))/sqrt(varSc)
pc<-(1-pnorm(abs(zeta)))*2

#
# partial test with correction for correlation among blocks
#

if (!missing(cv)){
XC<-rep(NA,times=(ny*maxblock))
dim(XC)<- c(ny,maxblock)
Rks<-X
varSccv<-0
for (i in 1:length(date)) XC[date[i]-minyear+1,block[i]]<-cv[i]
for (i in 1:maxblock){ 
R<-XC[,i]
R<-rank(R,na.last="keep",ties.method="average")
R[is.na(R)]<-mean(R[!is.na(R)])
Rks[,i]<-R
}
# 
# recalculating varmix for intrablock variance
#
varmix<-0
for (g in 1:maxblock){
for (h in 1:maxblock){
kkcv=XC[,g]
kky=X[,h]
Ry<-rank(kky,na.last="keep",ties.method="average")
Ry[is.na(Ry)]<-mean(Ry[!is.na(Ry)])
Rcv<-rank(kkcv,na.last="keep",ties.method="average")
Rcv[is.na(Rcv)]<-mean(Rcv[!is.na(Rcv)])
Kmix<-0
for (i in 2:length(kky)) {
for (j in 1:(i-1)) Kmix<-Kmix+sign1((kkcv[i]-kkcv[j])*(kky[i]-kky[j]))
}
Gmix<-0
for (i in 1:length(kky)) Gmix<-Gmix+Ry[i]*Rcv[i]
varmix<-varmix+(Kmix+4*Gmix-length(kky)*(length(kky[!is.na(kky)])+1)*(length(kkcv[!is.na(kkcv)])+1))/3
}
}
#
#
#
for (g in 1:maxblock) {
for (h in 1:maxblock){
K<-0
for (i in 2:ny) {
for (j in 1:(i-1)) K<-K+sign1((XC[j,g]-XC[i,g])*(XC[j,h]-XC[i,h]))
}
G<-0
for (i in 1:ny) G<-G+Rks[i,g]*Rks[i,h] 
xg<-X[,g]
xh<-X[,h]
varSccv<-varSccv+(K+4*G-length(xg)*(length(xg[!is.na(xg)])+1)*(length(xh[!is.na(xh)])+1))/3
}
}
p.varSc<-varSc-varmix*varmix/varSccv
zeta <- (p.S-sign(p.S))/sqrt(p.varSc)
p.pc <- (1-pnorm(abs(zeta)))*2
}
}
if (ncomp > 0) tau<-S/ncomp
ans=list(sl=p,S=S,B=B,varS=varS,sl.corrected=pc,varS.corrected=varSc,partial.S=p.S,partial.sl=p.p,partial.varS=p.varS, partial.sl.corrected=p.pc,partial.varS.corrected=p.varSc,tau=tau)
oldClass(ans)<-"rkt"
return(ans)

}

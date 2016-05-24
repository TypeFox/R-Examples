Steps <-
function(m,n,stage="all"){
lst<-NULL
V<-NULL
WW<-NULL
W<-NULL

A<-BIB(m)
s<-A$BIB
while (m >= 2) {
d<-dim(s)[2]
reso<-Resolvable(1,s)
U<-reso$RBIB
UD<-Uniform(U)
if (d > 3) {
ss<-Gen(n,s)
s<-ss$BIB2
WW<-list(WW,ss)}
W<-c(W,reso)
V<-c(V,UD)
m<-m-1}

BIB<-FALSE
BIBg<-FALSE
RBIB<-FALSE
UDs<-FALSE
if (length(stage) == 1 && stage == "all") {
stage <- c("S1", "S2", "S3", "S4")}
for (i in 1:length(stage)) {
stage_ <- stage[i]
        switch(stage_, S1 = {
            BIB <- TRUE
        }, S2 = {
            BIBg <- TRUE
        }, S3 = {
            RBIB <- TRUE
        }, S4 = {
            UDs <- TRUE
        })
}
if (BIB==TRUE) {
lst<-c(lst,BIB1=A)}
if (BIBg==TRUE){
lst<-c(lst,BIBg=WW)}
if (RBIB==TRUE) {
lst<-c(lst,Resolvables=W)}
if (UDs==TRUE){
lst<-c(lst,V)}
return(lst)}

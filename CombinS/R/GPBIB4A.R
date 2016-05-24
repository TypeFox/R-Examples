GPBIB4A <-
function(n,l,s,w){
if (s<2 & l<2 & n<2) {"n,l,s should be great than 1"}
else {
V<-n*l
reso<-(n*l)%%(2*s)
bbo<-ifelse(reso==0,"Yes","No")

A<-NULL;mat<-NULL;lamda<-NULL
for (i in 1:w){
A[[i]]<-matrix(1:V, ncol=l, byrow=TRUE)
z<-(i-1)*V
A[[i]]<-A[[i]]+z}

######################################################
##Op : Fonction interne pour les autres fonctions :
Op<-function(mt,s){
n<-dim(mt)[1];l<-dim(mt)[2];V<-n*l
a<-combn(1:l, s);b<-combn(1:n, 2)
v<-dim(b)[2];vv<-dim(a)[2]
MAT<-NULL;A<-1;y<-1
while(A<=vv) {
for (k in 1:v){
s<-a[,A];ss<-b[,k]
MAT[[y]]<-as.vector(t(mt[ss,s]))
y<-y+1}
A<-A+1}
return(Reduce("rbind",MAT))}
######################################################
######################################################

for (i in 1:w){
mw<-A[[i]]
mat[[i]]<-Op(mw,s)}

BIB<-Reduce("rbind",mat)
T<-BIB[1,1];R<-length(which(T==BIB))
lamda[1]<-(n-1)*choose(l-2,s-2);lamda[2]<-choose(l-1,s-1);lamda[3]<-choose(l-2,s-2);lamda[4]<-0
return(list(PBIB=BIB,Type="Generalized rectangular right angular (4) (GPBIB_4) design",V=w*V,B=dim(BIB)[1],R=R,K=2*s,lamda=lamda, Resolvable=bbo))}}

GPBIB7B <-
function(n,l,s,w){
if (s<3 & l<2 & n<2) {"n,l,s should be great than 1 and s great than 2"}
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
Bp<-NULL
BB<-NULL
for (j in 1:w) {
AA<-A[[j]]
AB<-A[-j]
M<-length(AB)


for (m in 1:M){
AS<-AB[[m]]

for (k in 1:l) {
AS1<-AS
AS2<-AS
co1<-AA[,k]
co2<-AA[k,]

AS1[,k]<-AA[,k]
AS2[k,]<-AA[k,]


X1<-Op(AS1,s)
y1<-dim(X1)[1]
for (x in y1:1){
if (length(intersect(X1[x,],co1))==0) {
X1<-X1[-x,]}}
Bp<-rbind(Bp,X1)

X2<-Op(AS2,s)
y2<-dim(X2)[1]
for (z in y2:1){
if (length(intersect(X2[z,],co2))==0) {
X2<-X2[-z,]}}

BB<-rbind(BB,X2)

}}}
PBIB<-rbind(Bp,BB)
T <- PBIB[1, 1]
R <- length(which(T == PBIB))
lamda[1] <- s*(n-1)*(w-1)*choose(l-2,s-2)
lamda[2] <- s*(w-1)*choose(l-1,s-1)
lamda[3] <- (l - 2) * choose(l - 3, s - 3)*(w-1)
lamda[4] <- 0
lamda[5] <- 2 * (n - 1) * choose(l - 2, s - 2)
lamda[6] <- 2 * choose(l - 1, s - 1)
lamda[7] <- 4 * choose(l - 2, s - 2)
return(list(PBIB = PBIB, Type = "Generalized rectangular right angular (7) (PBIB_7) design", V = w * V, B = dim(PBIB)[1], R = R, K = 2 * s, lamda = lamda, Resolvable=bbo))
}}

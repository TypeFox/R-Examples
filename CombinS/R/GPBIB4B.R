GPBIB4B <-
function(n,l,s,w){
if (s<3 & l<2 & n<2) {"n,l should be great than 1 and s great than 2"}
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
for (j in 1:w) {
AA<-A[[j]]
AB<-A[-j]
M<-length(AB)
vec<-NULL

for (m in 1:M){
AS<-AB[[m]]

for (k in 1:l) {
co<-AA[,k]
vec[[1]]<-AA[,k]

for (p in 2:n) {
vec[[p]]<-c(co[-1],co[1])
co<-c(co[-1],co[1])}
N<-length(vec)

for (p in 1:N) {

mt<-cbind(vec[[p]],AS)
X<-Op(mt,s)
y<-dim(X)[1]
nn<-vec[[p]]
for (x in y:1){
if (any(X[x,1]==nn)==FALSE){
X<-X[-x,]}}
Bp<-rbind(Bp,X)}}}}
PBIB<-Bp
T <- PBIB[1, 1]
R <- length(which(T == PBIB))
lamda[1] <- l * n * (n - 1) * choose(l - 2, s - 3)*(w-1)
lamda[2] <- n * (choose(l, s - 1) + (l * choose(l - 1,s - 2)))*(w-1)
lamda[3] <- n * l * choose(l - 2, s - 3)*(w-1)
lamda[4] <- 4 * (n - 1) * choose(l - 1, s - 2)
return(list(PBIB = PBIB, Type = "Generalized rectangular right angular (4) (GPBIB_4) design (lamda not equal to 0)", V = w * V, B = dim(PBIB)[1], R = R, K = 2 * s, lamda = lamda, Resolvable=bbo))
}}

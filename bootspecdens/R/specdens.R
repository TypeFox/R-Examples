`specdens` <-
function(data,h,B, level){
######
bootsample<- numeric(B)
q.oben<- numeric(6)

d<- dim(data)
m<- d[1]
n<- d[2]
N<- n*m

#######determine periodogram matrix

J <- function(lambda){
g=rep(0,m)
for(j in 1:n){
g <- g+(data[,j]*exp(-1i*lambda*j))
}
g/sqrt(2*pi*n)
}

I <- function(lambda){##periodogram matrix
J(lambda) %*% t(Conj(J(lambda)))
}

####### determine kernel estimator of the spectral density matrix
Kernel<- function(lambda){##Parzen kernel
y<-rep(0,length(lambda))
for(i in 1:length(lambda)){
if(lambda[i]==0)y[i]<-1
else y[i]<- (3/(8*pi)*(sin(lambda[i]/4)/(lambda[i]/4))^4)
}
return(y)
}

fHat <- function(lambda){##kernel estimator of the spectral density matrix
g=matrix(rep(0,m*m), ncol=m)
for(j in -floor((n-1)/2):(floor(n/2)) ){
g <- g+( Kernel((lambda-(2*pi*j/n))/h)*I(2*pi*j/n))
}
g/(n*h)
}
####
wHat <- function(lambda){##pooled kernel estimator
sum(diag(fHat(lambda)))/m
}

DiagfHat<- function(lambda){
diag(fHat(lambda))
}

Dhat <- function(lambda){
diag(wHat(lambda)*(1/DiagfHat(lambda)))
}

fHat_T<- function(lambda){
sqrt(Dhat(lambda))%*%fHat(lambda)%*%sqrt(Dhat(lambda))
} 

fHat_T_j<- fHat_T(2*pi*(-floor((n-1)/2))/n)
for(j in (-floor((n-1)/2)+1):(floor(n/2)) ){
fHat_T_j<- matrix(c(fHat_T_j, fHat_T(2*pi*j/n)),nrow=m)
}

fHat_T_sqrt_j<- svd(fHat_T_j[1:m,1:m])$u %*% sqrt(diag(svd(fHat_T_j[1:m,1:m])$d)) %*% svd(fHat_T_j[1:m,1:m])$v
for(j in seq(m+1,(m*n-1), by=m) ){
fHat_T_sqrt_j<- matrix(c(fHat_T_sqrt_j, svd(fHat_T_j[1:m,j:(j+(m-1))])$u %*% sqrt(diag(svd(fHat_T_j[1:m,j:(j+(m-1))])$d)) %*% svd(fHat_T_j[1:m,j:(j+(m-1))])$v), nrow=m)
}

fHat_Ttilde_j<-1/2*matrix(c(Re(fHat_T_j[1:m,1:m][c(1:m)]),Im(fHat_T_j[1:m,1:m])[c(1:m)],Re(fHat_T_j[1:m,1:m])[c((m+1):m^2)], Im(fHat_T_j[1:m,1:m])[c((m+1):m^2)],-Im(fHat_T_j[1:m,1:m])[c(1:m)], Re(fHat_T_j[1:m,1:m])[c(1:m)], -Im(fHat_T_j[1:m,1:m])[c((m+1):m^2)],Re(fHat_T_j[1:m,1:m])[c((m+1):m^2)]),ncol=2*m)
for(j in seq(m+1,(m*n-1), by=m) ){
fHat_Ttilde_j<- matrix(c(fHat_Ttilde_j, 1/2*matrix(c(Re(fHat_T_j[1:m,j:(j+(m-1))][c(1:m)]),Im(fHat_T_j)[1:m,j:(j+(m-1))][c(1:m)],Re(fHat_T_j)[1:m,j:(j+(m-1))][c((m+1):m^2)], Im(fHat_T_j)[1:m,j:(j+(m-1))][c((m+1):m^2)],-Im(fHat_T_j)[1:m,j:(j+(m-1))][c(1:m)], Re(fHat_T_j)[1:m,j:(j+(m-1))][c(1:m)], -Im(fHat_T_j)[1:m,j:(j+(m-1))][c((m+1):m^2)],Re(fHat_T_j)[1:m,j:(j+(m-1))][c((m+1):m^2)]),ncol=2*m)), nrow=2*m)
}

fHat_Ttilde_sqrt_j<- svd(fHat_Ttilde_j[1:(2*m),1:(2*m)])$u %*% sqrt(diag(svd(fHat_Ttilde_j[1:(2*m),1:(2*m)])$d)) %*% svd(fHat_Ttilde_j[1:(2*m),1:(2*m)])$v
for(j in seq((2*m+1),((2*m)*n-3),by=(2*m)) ){
fHat_Ttilde_sqrt_j <-matrix(c(fHat_Ttilde_sqrt_j, svd(fHat_Ttilde_j[1:(2*m),j:(j+(2*m-1))])$u %*% sqrt(diag(svd(fHat_Ttilde_j[1:(2*m),j:(j+(2*m-1))])$d)) %*% svd(fHat_Ttilde_j[1:(2*m),j:(j+(2*m-1))])$v), nrow=(2*m))
}

######## standardisation factors of the test statistic T_n

g_rs <- function(r,s){
delta <- function(r,s){
ifelse(r==s,1,0)
}
delta(r,s)-(1/m)
}

kappaHatsq_s1s2 <- function(s1,s2, lambda){
g=0
for(j in -floor((n-1)/2):(floor(n/2)-1)){
g <- g+( Kernel((lambda-(2*pi*j/n))/h)*I(2*pi*j/n)[s1,s2])
}
(Mod(g/(n*h)))^2 * (diag(fHat(lambda))[s1]* diag(fHat(lambda))[s2])^(-1)
}

######determine tau0Hatsq

c=0####sum over r1, r2; Riemann-sum over lambda; sum over s1, s2 in g_r1s1*g_r2,s2*kappaHatsq_s1s2
for(i in 1:m){
#
k=60
bb <- function(r1,r2){
#
a<- function(r1,r2,lambda){
p=0
for(i in 1:m){
for(j in 1:m){
p<- p+(g_rs(r1,i)*g_rs(r2,j)*kappaHatsq_s1s2(i,j,lambda))
}
}
p
}
#
g=0
for(i in -k:k){
g<- g+(a(r1,r2,pi*i/k))^2
}
g*(pi/k)
}
#
for(j in 1:m){
c<- c+bb(i,j)
}
}


####

d<- function(x,y){##kernel(x+y)
Kernel(x+y)
}

k=40
g=0##Riemann-sum over Kernel(x)*Kernel(x+y)
for(i in -k:k){
#
k=40
e<- function(y){
g=0
for(i in -k:k){
g<- g+Kernel(30*i/k)*d(30*i/k,y)
}
g*(30/k)
}
#
g<- g+(e(30*i/k))^2
}
f<- g*(30/k)


tau0Hatsq<- 1/(2*pi^2)*f*c

#######determine my_n
pp=0##Riemann-sum over Kern(x)^2
k=50
for(i in -k:k){
pp<- pp+(Kernel(i*10/k))^2
}
pp<- pp*(10/k)


ggg=0##sum over r
for(i in 1:m){
#
zzz<- function(r){##sum over s1
#
zz<- function(r,s1){##sum over s2
#
k=100    
z <- function(s1,s2){##Riemann-sum over kappadachquad_s1s2
g=0
for(i in -k:k){
g<- g+kappaHatsq_s1s2(s1,s2,pi*i/k)
}
g*(pi/k)
}
#
g=0
for(i in 1:m){
g<- g+g_rs(r,s1)*g_rs(r,i)*z(s1,i)
}
g
}
#
g=0
for(i in 1:m){
g<- g+zz(r,i)
}
g
}
#
ggg<- ggg+zzz(i)
}

my_n <- 1/(2*pi*sqrt(h))*pp*ggg

########determine test statistic T_n
g=0##sum over r
for(r in 1:m){
#
jj<- function(r){##Riemnann-sum over lambda
g=0
for(i in -floor((n-1)/2):floor(n/2)){
g<- g+((fHat(2*pi*i/n)[r,r]/wHat(2*pi*i/n))-1)^2
}
g*(2*pi/n)
}
#
g<- g+jj(r)
}

T_n <- g/m

Z_n<- (N*sqrt(h)*T_n-my_n)/sqrt(tau0Hatsq)##standardized test statistic


############bootstrap loop
for(b in 1:B){
#####
yy<- mvrnorm(1, rep(0,m), diag(1,m))##random N(0,I) vector

ifelse(floor(n/2)==(n/2),XX_j0<- fHat_T_sqrt_j[1:m, ((n*m)/2-(m-1)):(n*m/2)] %*% yy ,XX_j0<- fHat_T_sqrt_j[1:m, (((n*m)-(m-2))/2):((((n*m)-(m-2))/2)+(m-1))] %*% yy)

XX_jpi<- fHat_T_sqrt_j[1:m,(m*n-(m-1)):(m*n)] %*% yy

Istar_n_j0<- XX_j0 %*% t(XX_j0)

Istar_n_jpi<- XX_jpi %*% t(XX_jpi)

####

y<- mvrnorm(1, rep(0,2*m), diag(1,2*m))

X_j<- fHat_Ttilde_sqrt_j[1:(2*m),1:(2*m)] %*% y
for(j in seq((2*m+1),((2*m)*n-3), by=(2*m)) ){
X_j <- matrix(c(X_j, fHat_Ttilde_sqrt_j[1:(2*m),j:(j+(2*m-1))] %*% y), nrow=(2*m))
}

Z_j<- X_j[1:(2*m),1][c(1:m)]+X_j[1:(2*m),1][c((m+1):(2*m))]*1i
for(j in 2:n ){
Z_j<- c(Z_j, X_j[1:(2*m),j][c(1:m)]+X_j[1:(2*m),j][c((m+1):(2*m))]*1i)
}

Istar_n_j<- Z_j[1:m] %*% Conj(t(Z_j[1:m]))##complex Wishart distributed
for(j in seq(m+1,(m*n-1), by=m) ){
Istar_n_j <- matrix(c(Istar_n_j, Z_j[j:(j+(m-1))] %*% Conj(t(Z_j[j:(j+(m-1))]))), nrow=m)
}

###############
if(floor(n/2)==(n/2)){
g=0
for(j in 1:(n/2-1)){
g <- g+( Kernel(((2*pi*(-floor((n-1)/2))/n)-(2*pi*(j-(n/2))/n))/h)*Conj(diag(Istar_n_j[1:m,(m*n-m*j-(m-1)):(m*n-m*j)])))
}
g<- g+Kernel((2*pi*(-floor((n-1)/2))/n)/h)*diag(Istar_n_j0)
for(j in (n/2+1):(n-1)){
g <- g+( Kernel(((2*pi*(-floor((n-1)/2))/n)-(2*pi*(j-(n/2))/n))/h)*diag(Istar_n_j[1:m,(m*j-(m-1)):(m*j)]))
}
g<- g+Kernel(((2*pi*(-floor((n-1)/2))/n)-(2*pi*floor(n/2)/n))/h)*diag(Istar_n_j[1:m,(m*n-(m-1)):(m*n)])
fHatStar_j<- g/(n*h)

for(i in (-floor((n-1)/2)+1):(floor(n/2)) ){
g=0
for(j in 1:(n/2-1)){
g <- g+( Kernel(((2*pi*i/n)-(2*pi*(j-(n/2))/n))/h)*Conj(diag(Istar_n_j[1:m,(m*n-m*j-(m-1)):(m*n-m*j)])))
}
g<- g+Kernel((2*pi*i/n)/h)*diag(Istar_n_j0)
for(j in (n/2+1):(n-1) ){
g <- g+( Kernel(((2*pi*i/n)-(2*pi*(j-(n/2))/n))/h)*diag(Istar_n_j[1:m,(m*j-(m-1)):(m*j)]))
}
g<- g+Kernel(((2*pi*i/n)-(2*pi*floor(n/2)/n))/h)*diag(Istar_n_j[1:m,(m*n-(m-1)):(m*n)])
fHatStar_j <- c(fHatStar_j,g/(n*h))
}
} else{
g=0
g<- g+Kernel(((2*pi*(-floor((n-1)/2))/n)-(2*pi*floor(-(n-1)/2)/n))/h)*diag(Istar_n_j[1:m,(m*n-(m-1)):(m*n)])
for(j in 1:(((n-1)/2)-1)){
g <- g+( Kernel(((2*pi*(-floor((n-1)/2))/n)-(2*pi*(j-((n-1)/2))/n))/h)*Conj(diag(Istar_n_j[1:m,(m*n-m*j-(m-1)):(m*n-2*j)])))
}
g<- g+Kernel((2*pi*(-floor((n-1)/2))/n)/h)*diag(Istar_n_j0)
for(j in ((n+1)/2+1):(n-1)){
g <- g+( Kernel(((2*pi*(-floor((n-1)/2))/n)-(2*pi*(j-((n+1)/2))/n))/h)*diag(Istar_n_j[1:m,(m*j-(m-1)):(m*j)]))
}
g<- g+Kernel(((2*pi*(-floor((n-1)/2))/n)-(2*pi*floor(n/2)/n))/h)*diag(Istar_n_j[1:m,(m*n-(m-1)):(m*n)])
fHatStar_j<- g/(n*h)

for(i in (-floor((n-1)/2)+1):(floor(n/2)) ){
g=0
g<- g+Kernel(((2*pi*i/n)-(2*pi*floor(-(n-1)/2)/n))/h)*diag(Istar_n_j[1:m,(m*n-(m-1)):(m*n)])
for(j in 1:(((n-1)/2)-1)){
g <- g+( Kernel(((2*pi*i/n)-(2*pi*(j-((n-1)/2))/n))/h)*Conj(diag(Istar_n_j[1:m,(m*n-m*j-(m-1)):(m*n-m*j)])))
}
g<- g+Kernel((2*pi*i/n)/h)*diag(Istar_n_j0)
for(j in ((n+1)/2+1):(n-1)){
g <- g+( Kernel(((2*pi*i/n)-(2*pi*(j-((n+1)/2))/n))/h)*diag(Istar_n_j[1:m,(m*j-(m-1)):(m*j)]))
}
g<- g+Kernel(((2*pi*i/n)-(2*pi*floor(n/2)/n))/h)*diag(Istar_n_j[1:m,(m*n-(m-1)):(m*n)])
fHatStar_j <- c(fHatStar_j,g/(n*h))
}
}
####

wHatStar_j<- 1/m*sum(fHatStar_j[1:m])
for(j in seq((m+1),(m*n-1), by=m) ){
wHatStar_j<- c(wHatStar_j, (1/m*sum(fHatStar_j[j:(j+(m-1))])))
}

######statistic Tstar for determine quantile

p=0
for(r in 1:m){
g=0
for(j in 1:n){
g<- g+((fHatStar_j[(m*j-(m-1)):(m*j)][r]/wHatStar_j[j])-1)^2
}
g<- g*(2*pi/n)
p<- p+g
}

Tstar_n<-p/m

N<- m*n

Zstar_n <- (N*sqrt(h)*Tstar_n-my_n)/sqrt(tau0Hatsq)##standadized statistic

bootsample[b]<- Zstar_n

}

q.oben<- quantile(Re(bootsample), 1-level)

ifelse(Re((N*sqrt(h)*T_n-my_n)/sqrt(tau0Hatsq)) >= Re(q.oben), print("not equal spectral densities"), print("equal spectral desnities"))

pzaehler=0

for(i in 1:B){
	ifelse(Re((N*sqrt(h)*T_n-my_n)/sqrt(tau0Hatsq))>= Re(bootsample[i]), pzaehler<- pzaehler, pzaehler<- pzaehler+1)
	}

pvalue <- pzaehler/B
print("p-value:")
print(pvalue)
}


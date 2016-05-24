#-----------------------------------------------------------------------
#
# Statistical shape analysis routines
# written by Ian Dryden in R  (see http://cran.r-project.org) 
# (c) Ian Dryden
#     University of Nottingham
#                        2015
#
#          Version 1.1-11   Date: August 18, 2015
#
#----------------------------------------------------------------------
#
#
#
#
#
#
###########################################################################



############################################################
#
#  FUNCTIONS FOR CALCULATING NON-EUCLIDEAN MEANS AND DISTANCES
#  OF COVARIANCE MATRICES
#
############################################################
# Log Euclidean mean: Sigma_L

estLogEuclid<-function(S,weights=1){
M <- dim(S)[3]
if (length(weights)==1){
weights <- rep(1,times=M)
}
sum<-S[,,1]*0
for (j in 1:M){
eS<-eigen(S[,,j],symmetric=TRUE)
sum<-sum+
weights[j]*eS$vectors%*%diag(log(eS$values))%*%t(eS$vectors)/
            sum(weights)
}
ans<-sum
eL<-eigen(ans,symmetric=TRUE)
eL$vectors%*%diag(exp(eL$values))%*%t(eL$vectors)
}


estPowerEuclid<-function(S,weights=1,alpha=0.5){
M <- dim(S)[3]
if (length(weights)==1){
weights <- rep(1,times=M)
}
sum<-S[,,1]*0
for (j in 1:M){
eS<-eigen(S[,,j],symmetric=TRUE)
sum<-sum+
weights[j]*eS$vectors%*%diag((eS$values)**alpha)%*%t(eS$vectors)/
            sum(weights)
}
ans<-sum
eL<-eigen(ans,symmetric=TRUE)
eL$vectors%*%diag((eL$values)**(1/alpha))%*%t(eL$vectors)
}




# Riemannian (weighted mean) : Sigma_R

estLogRiem2<-function(S,weights=1){
M <- dim(S)[3]
if (length(weights)==1){
weights <- rep(1,times=M)
}
check<-9
tau<-1
Hold<-99999
mu<-estLogEuclid(S,weights)
while (check>0.0000000001){
ev<-eigen(mu,symmetric=TRUE)
logmu<- ev$vectors%*%diag(log(ev$values))%*%t(ev$vectors)
Hnew<-Re(Hessian2(S,mu,weights))
logmunew<- logmu+tau*Hnew
ev<-eigen(logmunew,symmetric=TRUE)
mu<- ev$vectors%*%diag(exp(ev$values))%*%t(ev$vectors)
check <- Re(Enorm(Hold)-Enorm(Hnew))
if (check < 0){ 
tau <- tau/2 
check <- 999999
}
Hold<-Hnew
}
mu
}

# Hessian used in calculating Sigma_R

Hessian2 <- function(S,Sigma,weights=1){
M <- dim(S)[3]
k<-dim(S)[1]
if (length(weights)==1){
weights <- rep(1,times=M)
}
ev0<-eigen(Sigma,symmetric=TRUE)
shalf <- ev0$vectors %*% (diag( sqrt((ev0$values)))) %*%t( ev0$vectors )  
sumit<-matrix(0,k,k)
for (i in 1:(M)){
ev2<- eigen(shalf%*%solve(S[,,i])%*%shalf,symmetric=TRUE)
sumit<-sumit + weights[i]*ev2$vectors %*% diag ( log((ev2$values)) ) %*% t( ev2$vectors )/sum(weights)
}
-sumit
}

# Euclidean : Sigma_E

estEuclid<-function(S,weights=1){
M <- dim(S)[3]
if (length(weights)==1){
weights <- rep(1,times=M)
}
sum<-S[,,1]*0
for (j in 1:M){
sum<-sum+S[,,j]*weights[j]/sum(weights)
}
sum
}
# Cholesky mean : Sigma_C

estCholesky<-function(S,weights=1){
M <- dim(S)[3]
if (length(weights)==1){
weights <- rep(1,times=M)
}
sum<-S[,,1]*0
for (j in 1:M){
sum<-sum+t(chol(S[,,j]))*weights[j]/sum(weights)
}
cc<-sum
cc%*%t(cc)
}


estSS<-function(S,weights=1){
M <- dim(S)[3]
k<-dim(S)[1]
H<-defh(k)
if (length(weights)==1){
weights <- rep(1,times=M)
}
Q<-array(0,c(k+1,k,M))
for (j in 1:M){
Q[,,j]<-t(H)%*%(rootmat(S[,,j]))
}
ans<-procWGPA(Q,fixcovmatrix=diag(k+1),scale=FALSE,reflect=TRUE,
                     sampleweights=weights)
H%*%ans$mshape%*%t(H%*%ans$mshape)
}


estShape<-function(S,weights=1){
M <- dim(S)[3]
k<-dim(S)[1]
H<-defh(k)
if (length(weights)==1){
weights <- rep(1,times=M)
}
Q<-array(0,c(k+1,k,M))
for (j in 1:M){
Q[,,j]<-t(H)%*%(rootmat(S[,,j]))
}
ans<-procWGPA(Q,fixcovmatrix=diag(k+1),scale=TRUE,reflect=TRUE,
                     sampleweights=weights)
H%*%ans$mshape%*%t(H%*%ans$mshape)
}



estRiemLe <- function(S,weights){
M<-dim(S)[3]
k<-dim(S)[1]
if (M!=2) print("Sorry - Calculation not implemented for M>2 yet") 
if (M == 2){
P1 <- S[,,1]
P2<- S[,,2]
detP1<- prod(eigen(P1)$values)
detP2<- prod(eigen(P2)$values)
P1 <- P1/(detP1)^(1/k)
P2 <- P2/(detP2)^(1/k)
P1inv <- solve(P1)
P12sq <- P1inv %*% P2%*%P2 %*% P1inv
tem <- eigen(P12sq,symmetric=TRUE)
A2 <- tem$vectors %*% diag( log(tem$values) ) %*%t(tem$vectors)
logPs2 <- weights[2]*A2
tem2<- eigen(logPs2,symmetric=TRUE)
Ps2 <- tem2$vectors %*% diag( exp(tem2$values) ) %*%t(tem2$vectors)
P12s <- P1 %*% Ps2 %*% P1
tem3 <- eigen(P12s,symmetric=TRUE)
P12sA <- tem3$vectors %*% diag( sqrt(tem3$values) ) %*% t( tem3$vectors )
Ptildes <- (detP1*(detP2/detP1)^weights[2])^(1/k)*P12sA
Ptildes
}
}






##########distances#################################

distRiemPennec <- function( P1, P2){
eig <- eigen(P1,symmetric=TRUE)
P1half <- eig$vectors %*% diag(sqrt(eig$values)) %*% t(eig$vectors)
P1halfinv <- solve(P1half)
AA <- P1halfinv%*%P2%*%P1halfinv
tem <- eigen(AA,symmetric=TRUE)
A2 <- tem$vectors %*% diag( log(tem$values) ) %*%t(tem$vectors)
dd <- Enorm(A2)
dd
}

distLogEuclidean <- function( P1, P2){
eig <- eigen(P1,symmetric=TRUE)
logP1<- eig$vectors%*%diag( log(eig$values) ) %*%t(eig$vectors)
tem <- eigen(P2,symmetric=TRUE) 
logP2 <- tem$vectors %*% diag( log(tem$values) ) %*%t(tem$vectors)
dd <- Enorm( logP1-logP2)
dd
}

distRiemannianLe <-
function( P1, P2){
dd <- distRiemPennec( P1%*%t(P1), P2%*%t(P2))/2
dd
}


distProcrustesSizeShape<- 
function (P1, P2) 
{
    H <- defh(dim(P1)[1])
    Q1 <- t(H) %*% rootmat(P1)
    Q2 <- t(H) %*% rootmat(P2)
    ans <- sqrt(centroid.size(Q1)^2 + centroid.size(Q2)^2 - 2 * 
        centroid.size(Q1) * centroid.size(Q2) * cos(riemdist(Q1, 
        Q2, reflect = TRUE)))
    ans
}

distProcrustesFull <- function(P1,P2){
H<-defh(dim(P1)[1])
Q1<-t(H)%*%rootmat(P1) 
Q2<-t(H)%*%rootmat(P2) 
ans<- riemdist(Q1,Q2,reflect=TRUE)
ans
}



distPowerEuclidean<-function(P1,P2,alpha=1/2){
if (alpha != 0){
eS<-eigen(P1,symmetric=TRUE)
Q1<- eS$vectors%*%diag(abs(eS$values)^alpha)%*%t(eS$vectors)
eS<-eigen(P2,symmetric=TRUE)
Q2<- eS$vectors%*%diag(abs(eS$values)^alpha)%*%t(eS$vectors)
dd <- Enorm( Q1 - Q2 )/abs(alpha)
}
if (alpha == 0){
dd <- distLogEuclidean(P1,P2)
}
dd
}

distCholesky<-function(P1,P2){
H<-defh(dim(P1)[1])
Q1<-t(H)%*%t(chol(P1))
Q2<-t(H)%*%t(chol(P2))
ans <- Enorm( Q1 - Q2)
ans
}

distEuclidean<-function(P1,P2){
ans <- Enorm( P1 - P2 )
ans
}

##################


distcov <- function( S1, S2 , method="Riemannian",alpha=1/2){
if (method=="Procrustes"){
dd <- distProcrustesSizeShape(S1,S2)
}
if (method=="ProcrustesShape"){
dd <- distProcrustesFull(S1,S2)
}
if (method=="Riemannian"){
dd <- distRiemPennec(S1,S2)
}
if (method=="Cholesky"){
dd <- distCholesky(S1,S2)
}
if (method == "Power"){
dd <- distPowerEuclidean(S1,S2,alpha)
}
if (method == "Euclidean"){
dd <- distEuclidean(S1,S2)
}
if (method == "LogEuclidean"){
dd <- distLogEuclidean(S1,S2)
}
if (method=="RiemannianLe"){
dd <- distRiemannianLe(S1,S2)
}
dd
}


estcov <-
function (S, method = "Riemannian", weights = 1, alpha = 1/2, MDSk = 2) 
{
    out <- list(mean = 0, sd = 0, pco = 0, eig = 0, dist=0)
    M <- dim(S)[3]
    if (length(weights) == 1) {
        weights <- rep(1, times = M)
    }
    if (method == "Procrustes") {
        dd <- estSS(S, weights)
    }
    if (method == "ProcrustesShape") {
        dd <- estShape(S, weights)
    }
    if (method == "Riemannian") {
        dd <- estLogRiem2(S, weights)
    }
    if (method == "Cholesky") {
        dd <- estCholesky(S, weights)
   }
    if (method == "Power") {
        dd <- estPowerEuclid(S, weights, alpha)
    }
    if (method == "Euclidean") {
        dd <- estEuclid(S, weights)
    }
    if (method == "LogEuclidean") {
        dd <- estLogEuclid(S, weights)
    }
    if (method == "RiemannianLe") {
        dd <- estRiemLe(S, weights)
    }
    out$mean <- dd
    sum <- 0
    for (i in 1:M) {
        sum <- sum + weights[i] * distcov(S[, , i], dd, method=method)^2/sum(weights)
    }
    out$sd <- sqrt(sum)
    dist <- matrix(0, M, M)
    for (i in 2:M) {
        for (j in 1:(i - 1)) {
            dist[i, j] <- distcov(S[, , i], S[, , j], method = method)
            dist[j, i] <- dist[i, j]
        }
    }
    out$dist <- dist
    if (M > MDSk){
    ans <- cmdscale(dist, k = MDSk, eig = TRUE, add = TRUE, x.ret = TRUE)
    out$pco <- ans$points
    out$eig <- ans$eig
    if (MDSk>2){
    shapes3d(out$pco[,1:min(MDSk,3)],axes3=TRUE)
    }
if (MDSk==2){
plot(out$pco,type="n",xlab="MDS1",ylab="MDS2")
text(out$pco[,1],out$pco[,2],1:length(out$pco[,1]))
}
}
    out
}

rootmat<- function( P1 ){
        eS <- eigen(P1, symmetric = TRUE)
        if (min(eS$values) < -0.001) {
print("Not positive-semi definite") 
}
else{
        Q1 <- eS$vectors %*% diag(sqrt(abs(eS$values))) %*% t(eS$vectors)
Q1
}
}



##########################

shapes.cva <- function( X , groups , scale=TRUE, ncv=2){
g <- dim( table ( groups ) )
ans <- procGPA( X , scale=scale)

if (scale==TRUE) pp <- (ans$k-1)*ans$m - (ans$m*(ans$m-1)/2) - 1
if (scale==FALSE) pp <- (ans$k-1)*ans$m - (ans$m*(ans$m-1)/2) 

pracdim <- min( pp,  ans$n - g)
out <- lda( ans$scores[,1:pracdim] , groups)
print((out))
cv <- predict(out,dimen=3)$x
if (dim(cv)[2]==1) {
cv <- cbind( cv, rnorm(dim(cv)[1])/1000 )
}
if (ncv==2){
eqscplot(cv,type="n",xlab="CV1",ylab="CV2")
text(cv, labels=groups)
}
if (ncv==3){
shapes3d(cv,color=groups,axes3=TRUE)
}
cv
}


groupstack<-function( A1, A2, A3=0, A4=0, A5=0, A6=0, A7=0, A8=0){
out<-list(x=0,groups="")
dat <- abind(A1,A2)
group <- c(rep(1,times=dim(A1)[3]),rep(2,times=dim(A2)[3]))
if (is.array(A3)){
dat <- abind(dat,A3)
group <- c(group,rep(3,times=dim(A3)[3]))
if (is.array(A4)){
dat <- abind(dat,A4)
group <- c(group,rep(4,times=dim(A4)[3]))
if (is.array(A5)){
dat <- abind(dat,A5)
group <- c(group,rep(5,times=dim(A5)[3]))
if (is.array(A6)){
dat <- abind(dat,A6)
group <- c(group,rep(6,times=dim(A6)[3]))
if (is.array(A7)){
dat <- abind(dat,A7)
group <- c(group,rep(7,times=dim(A7)[3]))
if (is.array(A8)){
dat <- abind(dat,A8)
group <- c(group,rep(8,times=dim(A8)[3]))
}
}
}
}
}
}
out$x <- dat
out$groups <- group
out
}

 



###########################


transformations<-function(Xrotated,Xoriginal){
# outputs the translations, rotations and 
# scalings for ordinary Procrustes rotation 
# of each individual in Xoriginal to the 
# Procrustes rotated individuals in Xrotated
X1<-Xrotated
X2<-Xoriginal
n<-dim(X1)[3]
m<-dim(X1)[2]
translation<-matrix(0,m,n)
scale<-rep(0,times=n)
rotation<-array(0,c(m,m,n))
for (i in 1:n){
translation[,i]<- -apply(X2[,,i]-X1[,,i],2,mean)
ans<-procOPA(X1[,,i],X2[,,i])
scale[i]<-ans$s
rotation[,,i]<-ans$R
}
out<-list(translation=0,scale=0,rotation=0)
out$translation<-translation
out$scale<-scale
out$rotation<-rotation
out
}


iglogl<- function( x ,lam, nlam){
gamma<- abs(x[1])
alpha<- gamma/ mean(1/lam[1:nlam])
ll <-  - (gamma+1) * sum( log(lam[1:nlam] )) - alpha * sum (1/lam[1:nlam]) +
              nlam * gamma *log(alpha) - nlam * lgamma (gamma) 
-ll
}


procWGPA<- function( x, fixcovmatrix=FALSE, initial="Identity", maxiterations=10, 
scale=TRUE, reflect=FALSE, prior="Exponential",diagonal=TRUE,
  sampleweights="Equal"){
X<-x
priorargument<-prior
alpha<-"not estimated"
gamma<-"not estimated"
k<-dim(X)[1]
n<-dim(X)[3]
m<-dim(X)[2]

if (initial[1]=="Identity"){
Sigmak <- diag(k)
}
else{
if (initial[1]=="Rawdata"){
tol<-0.0000000001
if (m==2){
Sigmak<- diag( diag(var( t(X[,1,])))+diag(var( t(X[,2,]))) )/2 + tol
}
if (m==3){
Sigmak<- diag( diag(var( t(X[,1,])))+diag(var( t(X[,2,])))+diag(var( t(X[,3,]))) )/3 + tol
}
}
else
{
Sigmak <- initial
}
}

mu<-procGPA(X,scale=scale)$mshape



#cat("Iteration 1 \n")

if (fixcovmatrix[1]!=FALSE){
Sigmak <- fixcovmatrix
}

ans<-procWGPA1(X,mu,metric=Sigmak,scale=scale,reflect=reflect,sampleweights=sampleweights)


if ((maxiterations>1)&&(fixcovmatrix[1]==FALSE)){

ans0<-ans
dif <- 999999
it<-1
while ((dif > 0.00001)&&(it< maxiterations)){

it <- it + 1
if (it==2){
cat("Differences in norm of Sigma estimates... \n ")
}

if (prior[1]=="Identity"){
prior <- diag(k)
}

if (prior[1]=="Inversegamma"){
lam <- eigen(ans$Sigmak)$values
nlam<- min( c( n*m-m-3, k-3 ) )
mu <- mean(1/lam[1:( nlam)])
alpha <- 1/mu
out <- nlm( iglogl, p=c(1) ,lam=lam, nlam=nlam)
#print(out)
gamma <- abs(out$estimate[1])
alpha<- gamma/ mean(1/lam[1:nlam])
newmetric <- n*m/(n*m+2*(1+gamma))* ( ans$Sigmak + (2*alpha/(n*m))*diag(k) )
#dif2<-999999
#while (dif2> 0.000001){
#old<-alpha
#lam <- eigen(newmetric)$values
#out <- nlm( iglogl, p=c(1) ,lam=lam, nlam=nlam)
#gamma <- abs(out$estimate[1])
#alpha<- gamma/ mean(1/lam[1:nlam])
#newmetric <- n*m/(n*m+2*(1+gamma))*( ans$Sigmak + (2*alpha/(n*m))*diag(k) )
#dif2<- abs(alpha- old)
#print(dif2)
#}
}


if (prior[1]=="Exponential"){
lam <- eigen(ans$Sigmak)$values
nlam<- min( c( n*m-m-2, k-2 ) )
mu <- mean(1/lam[1:( nlam)])
alpha <- 1/mu
gamma<-1
newmetric <- n*m/(n*m+2*(2))*( ans$Sigmak + (2*alpha/(n*m))*diag(k) )
#dif2<-999999
#while (dif2> 0.000001){
#old<-alpha
#newmetric <- n*m/(n*m+2*(2))*( ans$Sigmak + (2*alpha/(n*m))*diag(k) )
#lam <- eigen(newmetric)$values
#mu <- mean(1/lam[1:( nlam)])
#alpha <- 1/mu
#newmetric <- n*m/(n*m+2*(2))*( ans$Sigmak + (2*alpha/(n*m))*diag(k) )
#dif2<- abs(alpha- old)
#}
}

if (is.double(prior[1])){
newmetric <- ( ans$Sigmak + prior )
}



if (diagonal==TRUE){
newmetric<- diag(diag(newmetric))
}

if (fixcovmatrix[1]!=FALSE){
newmetric <- fixcovmatrix
}

ans2<-procWGPA1(X, ans$mshape,  metric= newmetric ,scale=scale,sampleweights=sampleweights)
plotshapes(ans2$rotated)
dif<-Enorm( (ans$Sigmak - ans2$Sigmak) )
ans<-ans2
cat(c(it," ",dif," \n"))
}
}
if ((priorargument[1]=="Exponential")||(priorargument[1]=="Inversegamma")){
ans$alpha<-alpha
ans$gamma<-gamma
}
cat(" \n")
ans
}





procWGPA1 <- function( X, mu, metric="Identity", scale=TRUE, 
          reflect=FALSE, sampleweights="Equal"){



k<-dim(X)[1]
n<-dim(X)[3]
m<-dim(X)[2]

sum <- 0
for (i in 1:n){
sum <- sum + centroid.size(X[,,i])**2
}
size1 <- sqrt(sum)





if (sampleweights[1]=="Equal"){
sampleweights<- rep(1/n, times=n)
}


if (length(sampleweights)!= n){
cat("Sample weight vector not of correct length \n")
}


if (metric[1]=="Identity"){
Sigmak <- diag(k)
}
else{
Sigmak <- metric
}

eig<-eigen(Sigmak,symmetric=TRUE)

Sighalf <- eig$vectors %*% diag (sqrt(abs(eig$values)))%*%t(eig$vectors)
Siginvhalf <- eig$vectors %*% diag(1/sqrt(abs(eig$values)))%*%t(eig$vectors)
Siginv<-eig$vectors %*% diag (1/(eig$values))%*%t(eig$vectors)

one <- matrix(rep(1,times=k),k,1)

Xstar <- X



for (i in 1:n){
Xstar[,,i]<- Xstar[,,i]- one%*%t(one)%*%Siginv %*% Xstar[,,i]/
                                     c( t(one)%*%Siginv%*%one )
Xstar[,,i] <- Siginvhalf%*%Xstar[,,i]
}

mu<- mu - one%*%t(one)%*%Siginv %*% mu /c( t(one)%*%Siginv%*%one )


ans<-procGPA(Xstar,eigen2d=FALSE)
ans2<-ans

dif3 <- 99999999
while (dif3 > 0.00001){

for (i in 1:n){
old<-mu
tem<- procOPA( Siginvhalf%*%mu , Xstar[,,i],scale=scale,reflect=reflect)
Gammai<- tem$R
betai <- tem$s
#ci <- t(one)%*% Siginvhalf %*% X[,,i] %*% Gammai*betai/k
#Yi <- Sighalf%*% ans$rotated[,,i] + Sighalf%*%one%*% ci
#Zi <- Yi - one %*% t(one)%*% Siginv %*% Yi / c( t(one)%*%Siginv%*%one )

Zi <- Sighalf%*% Xstar[,,i]%*%Gammai*betai
ans2$rotated[,,i]<-Zi
}

sum2 <- 0
for (i in 1:n){
sum2 <- sum2 + centroid.size(ans2$rotated[,,i])**2
}
size2 <- sqrt(sum2)



tem<-ans2$mshape*0
for (i in 1:n){
ans2$rotated[,,i] <- ans2$rotated[,,i]*size1/size2
tem<-tem+ ans2$rotated[,,i]*sampleweights[i]/sum(sampleweights)
}

mu<- tem
dif3<- riemdist( old,  mu)
} 





z<-ans2
z$mshape <- tem

    tan <- z$rotated[, 1, ] - z$mshape[, 1]
    for (i in 2:m) {
        tan <- rbind(tan, z$rotated[, i, ] - z$mshape[, i])
    }
    pca<-prcomp1(t(tan))
    z$tan <- tan
    npc <- 0
    for (i in 1:length(pca$sdev)) {
        if (pca$sdev[i] > 1e-07) {
            npc <- npc + 1
        }
    }
    z$scores <- pca$x
    z$rawscores <- pca$x
    z$stdscores <-pca$x
    for (i in 1:npc) {
        z$stdscores[, i] <- pca$x[, i]/pca$sdev[i]
    }
    z$pcar <- pca$rotation
    z$pcasd <- pca$sdev
    z$percent <- z$pcasd^2/sum(z$pcasd^2) * 100
    
    size <- rep(0, times = n)
    rho <- rep(0, times = n)
    x <- X
    size <- apply(x, 3, centroid.size)
    rho <- apply(x, 3, y <- function(x) {
        riemdist(x, z$mshape)
    })

    z$rho <- rho
    z$size <- size
    z$rmsrho <- sqrt(mean(rho^2))
    z$rmsd1 <- sqrt(mean(sin(rho)^2))
   
    z$k <- k
    z$m <- m
    z$n <- n


tem<- matrix(0, k,k)

for (i in 1:n){
tem <- tem +  (z$rotated[,,i] - z$mshape)%*%t((z$rotated[,,i] - z$mshape))
}
tem<-tem/(n*m)
    z$Sigmak <-  tem 


return(z)
}

testshapes <- function( A, B, resamples = 1000, replace=TRUE, scale=TRUE){
if (replace==TRUE){
out<-bootstraptest(A,B,resamples=resamples, scale=scale)
}
if (replace==FALSE){
out<-permutationtest(A,B,nperms=resamples, scale=scale)
}
out
}


testmeanshapes <- function( A, B, resamples = 1000, replace=FALSE, scale=TRUE){
if (replace==TRUE){
out<-bootstraptest(A,B,resamples=resamples, scale=scale)
}
if (replace==FALSE){
out<-permutationtest(A,B,nperms=resamples, scale=scale)
}


if (resamples > 0){

aa <- list(H = 0, H.pvalue=0, H.table.pvalue = 0, 
           G = 0, G.pvalue = 0, G.table.pvalue = 0,
            J = 0, J.pvalue = 0, J.table.pvalue = 0)

aa$H <- out$H 
aa$H.pvalue <- out$H.pvalue 
aa$H.table.pvalue <- out$H.table.pvalue 
aa$G <- out$G 
aa$G.pvalue <- out$G.pvalue 
aa$G.table.pvalue <- out$G.table.pvalue 
aa$J <- out$J 
aa$J.pvalue <- out$J.pvalue 
aa$J.table.pvalue <- out$J.table.pvalue 

 }

if (resamples == 0) {
 
aa <- list(H = 0, H.table.pvalue = 0, 
           G = 0, G.table.pvalue = 0,
            J = 0, J.table.pvalue = 0)


aa$H <- out$H 
aa$H.table.pvalue <- out$H.table.pvalue 
aa$G <- out$G 
aa$G.table.pvalue <- out$G.table.pvalue 
aa$J <- out$J 
aa$J.table.pvalue <- out$J.table.pvalue 

}

aa
}



permutationtest2<-function (A, B, nperms = 1000, scale=scale)
{
    A1 <- A
    A2 <- B
mdim <- dim(A1)[2]
    B <- nperms
    nsam1 <- dim(A1)[3]
    nsam2 <- dim(A2)[3]
    pool<-procGPA( abind (A1, A2) , scale=scale, tangentcoords = "partial", pcaoutput=FALSE)

    tempool <- pool
for (i in 1:(nsam1+nsam2)){
tempool$tan[,i] <- pool$tan[,i]/Enorm(pool$tan[,i])*pool$rho[i]
}
pool<-tempool

    permpool<-pool
    Gtem <- Goodall( pool, nsam1, nsam2)
    Htem <- Hotelling( pool, nsam1, nsam2)
    Jtem <- James( pool, nsam1, nsam2, table=TRUE)
    Ltem <- Lambdamin( pool, nsam1, nsam2)
    Gumc <- Gtem$F
    Humc <- Htem$F
    Jumc <- Jtem$Tsq
    Lumc <- Ltem$lambdamin
    Gtabpval <- Gtem$pval
    Htabpval <- Htem$pval
    Jtabpval <- Jtem$pval
    Ltabpval <- Ltem$pval
    if (B > 0) {
        Apool <- array(0, c(dim(A1)[1], dim(A1)[2], dim(A1)[3] +
            dim(A2)[3]))
        Apool[, , 1:nsam1] <- A1
        Apool[, , (nsam1 + 1):(nsam1 + nsam2)] <- A2
        out <- list(H = 0, H.pvalue = 0, H.table.pvalue = 0, J = 0, J.pvalue = 0, J.table.pvalue = 0,
            G = 0, G.pvalue = 0, G.table.pvalue = 0)
        Gu <- rep(0, times = B)
        Hu <- rep(0, times = B)
        Ju <- rep(0, times = B) 
        Lu <- rep(0, times = B)       
cat("Permutations - sampling without replacement: ")
        cat(c("No of permutations = ", B, "\n"))
        for (i in 1:B) {
            if (i/100 == trunc(i/100)){
            cat(c(i, " ")) }
            select <- sample(1:(nsam1 + nsam2))
            permpool$tan <- pool$tan[,select]
            Gu[i] <- Goodall(permpool, nsam1, nsam2 )$F
            Hu[i] <- Hotelling(permpool, nsam1, nsam2 )$F
            Ju[i] <- James(permpool, nsam1, nsam2 )$Tsq        
            Lu[i] <- Lambdamin(permpool, nsam1, nsam2 )$lambdamin      
}
        Gu <- sort(Gu)
        numbig <- length(Gu[Gumc < Gu])
        pvalG <- (1 + numbig)/(B + 1)
        Lu <- sort(Lu)
        numbig <- length(Lu[Lumc < Lu])
        pvalL <- (1 + numbig)/(B + 1)
        Ju <- sort(Ju)
        numbig <- length(Ju[Jumc < Ju])
        pvalJ <- (1 + numbig)/(B + 1)
        Hu <- sort(Hu)
        numbig <- length(Hu[Humc < Hu])
        pvalH <- (1 + numbig)/(B + 1)
        cat(" \n")
        out$Hu <- Hu
        out$Ju <- Ju
        out$Gu <- Gu
        out$Lu <- Lu
        out$H <- Humc
        out$H.pvalue <- pvalH
        out$H.table.pvalue <- Htabpval
        out$J <- Jumc
        out$J.pvalue <- pvalJ
        out$J.table.pvalue <- Jtabpval
        out$G <- Gumc
        out$G.pvalue <- pvalG
        out$G.table.pvalue <- Gtabpval
        out$L <- Lumc
        out$L.pvalue <- pvalL
        out$L.table.pvalue <- Ltabpval
 }
    if (B == 0) {
        out <- list(H = 0, H.table.pvalue = 0, G = 0, G.table.pvalue = 0)
        out$H <- Humc
        out$H.table.pvalue <- Htabpval
        out$J <- Jumc
        out$J.table.pvalue <- Jtabpval
        out$G <- Gumc
        out$G.table.pvalue <- Gtabpval
        out$L <- Lumc
        out$L.table.pvalue <- Ltabpval
}
    out
}


bootstraptest<-function (A, B, resamples = 200, scale=TRUE)
{
    A1 <- A
    A2 <- B
mdim <- dim(A1)[2]
    B <- resamples
    nsam1 <- dim(A1)[3]
    nsam2 <- dim(A2)[3]
    pool<-procGPA( abind (A1, A2) ,scale=scale , tangentcoords="partial", pcaoutput=FALSE)
    tempool <- pool
for (i in 1:(nsam1+nsam2)){
tempool$tan[,i] <- pool$tan[,i]/Enorm(pool$tan[,i])*pool$rho[i]
}
pool<-tempool
    bootpool<-pool
    Gtem <- Goodall( pool, nsam1, nsam2)
    Htem <- Hotelling( pool, nsam1, nsam2)
    Jtem <- James( pool, nsam1, nsam2, table=TRUE)
    Ltem <- Lambdamin( pool, nsam1, nsam2 )
    Gumc <- Gtem$F
    Humc <- Htem$F
    Jumc <- Jtem$Tsq
    Lumc <- Ltem$lambdamin
    Gtabpval <- Gtem$pval
    Htabpval <- Htem$pval
    Jtabpval <- Jtem$pval
    Ltabpval <- Ltem$pval
    if (B > 0) {
        Apool <- array(0, c(dim(A1)[1], dim(A1)[2], dim(A1)[3] +
            dim(A2)[3]))
        Apool[, , 1:nsam1] <- A1
        Apool[, , (nsam1 + 1):(nsam1 + nsam2)] <- A2
        out <- list(H = 0, H.pvalue = 0, H.table.pvalue = 0,J = 0, J.pvalue = 0, J.table.pvalue = 0,
            G = 0, G.pvalue = 0, G.table.pvalue = 0)
        Gu <- rep(0, times = B)
        Hu <- rep(0, times = B)
       Ju <- rep(0, times = B)
       Lu <- rep(0, times = B)
pool2<-pool
pool2$tan[,1:nsam1] <- pool$tan[,1:nsam1] - apply(pool$tan[,1:nsam1],1,mean)
pool2$tan[,(nsam1+1):(nsam1+nsam2)] <- pool$tan[,(nsam1+1):(nsam1+nsam2)] - 
                    apply(pool$tan[,(nsam1+1):(nsam1+nsam2)],1,mean)


        cat("Bootstrap - sampling with replacement within each group under H0: ")
        cat(c("No of resamples = ", B, "\n"))
        for (i in 1:B) {
            if (i/100 == trunc(i/100)){
            cat(c(i, " ")) }
            select1 <- sample(1:nsam1,replace=TRUE)
            select2 <- sample((nsam1+1):(nsam1+nsam2),replace=TRUE)
            bootpool$tan <- pool2$tan[,c(select1,select2)]
            Gu[i] <- Goodall( bootpool, nsam1, nsam2 )$F
            Hu[i] <- Hotelling( bootpool, nsam1, nsam2 )$F
            Ju[i] <- James( bootpool, nsam1, nsam2 )$Tsq
            Lu[i] <- Lambdamin( bootpool, nsam1, nsam2 )$lambdamin
       }
        Gu <- sort(Gu)
        numbig <- length(Gu[Gumc < Gu])
        pvalG <- (1 + numbig)/(B + 1)
        Ju <- sort(Ju)
        numbig <- length(Ju[Jumc < Ju])
        pvalJ <- (1 + numbig)/(B + 1)
        Hu <- sort(Hu)
        numbig <- length(Hu[Humc < Hu])
        pvalH <- (1 + numbig)/(B + 1)
        numbig <- length(Lu[Lumc < Lu])
        pvalL <- (1 + numbig)/(B + 1)
       cat(" \n")
       out$Hu <- Hu
        out$Ju <- Ju
        out$Gu <- Gu
        out$Lu <- Lu
        out$H <- Humc
        out$H.pvalue <- pvalH
        out$H.table.pvalue <- Htabpval
        out$J <- Jumc
        out$J.pvalue <- pvalJ
        out$J.table.pvalue <- Jtabpval
        out$G <- Gumc
        out$G.pvalue <- pvalG
        out$G.table.pvalue <- Gtabpval
        out$L <- Lumc
        out$L.pvalue <- pvalL
        out$L.table.pvalue <- Ltabpval
}
    if (B == 0) {
        out <- list(H = 0, H.table.pvalue = 0, G = 0, G.table.pvalue = 0, J = 0, J.table.pvalue = 0)
        out$H <- Humc
        out$H.table.pvalue <- Htabpval
        out$J <- Jumc
        out$J.table.pvalue <- Jtabpval
        out$G <- Gumc
        out$G.table.pvalue <- Gtabpval
        out$L <- Lumc
        out$L.table.pvalue <- Ltabpval
}
    out
}


 Lambdamin<-
function (pool, n1, n2, p = 0) 
{
    censiz <- centroid.size(pool$mshape)
    tan1 <- pool$tan[, 1:n1]
    tan2 <- pool$tan[, (n1 + 1):(n1 + n2)]
    kt <- dim(tan1)[1]
    n <- n1 + n2
    k <- pool$k
    m <- pool$m
    if (p == 0) {
        p <- min(k * m - (m * (m - 1))/2 - 1 - m, n1 + n2 - 2)
    }


HH <- diag(k)


mu1<-pool$mshape

if (dim(tan1)[1] == k*m - m){
 HH <- defh(k-1)

mu1<-preshape(pool$mshape)
}


    if (m == 2) {
        mu <- c(mu1[,1], mu1[, 2])      

    }
    if (m == 3) {
        mu <- c(mu1[, 1], mu1[, 2], mu1[, 
            3])

    }





    dd <- kt
    X1 <- tan1 * 0
    X2 <- tan2 * 0
    S1 <- matrix(0, dd, dd)
    S2 <- matrix(0, dd, dd)
    for (i in 1:n1) {
        X1[, i] <- (mu + tan1[, i])/Enorm(mu + tan1[, i])
        S1 <- S1 + X1[, i] %*% t(X1[, i])
    }
    for (i in 1:n2) {
        X2[, i] <- (mu + tan2[, i])/Enorm(mu + tan2[, i])
        S2 <- S2 + X2[, i] %*% t(X2[, i])
    }

sumx1 <- 0
sumx2 <- 0
for (i in 1:n1){
sumx1 <- sumx1 + X1[,i]
}
for (i in 1:n2){
sumx2 <- sumx2 + X2[,i]
}


    sum1 <- apply(X1, 1, sum)
    sum2 <- apply(X2, 1, sum)
    mean1 <- sum1/Enorm(sum1)
    mean2 <- sum2/Enorm(sum2)


    bb1 <- mean1[1:(dd - 1)]
    cc1 <- mean1[dd]
    bb2 <- mean2[1:(dd - 1)]
    cc2 <- mean2[dd]
    A1 <- cc1/abs(cc1) * diag(dd - 1) - cc1/(abs(cc1) + cc1^2) * 
        bb1 %*% t(bb1)
    M1 <- cbind(A1, -bb1)
    A1 <- cc2/abs(cc2) * diag(dd - 1) - cc2/(abs(cc2) + cc2^2) * 
        bb2 %*% t(bb2)
    M2 <- cbind(A1, -bb2)
    G1 <- matrix(0, dd - 1, dd - 1)
    G2 <- matrix(0, dd - 1, dd - 1)
    for (iu in 1:(dd - 1)) {
        for (iv in iu:(dd - 1)) {
            G1[iu, iv] <- G1[iu, iv] + t((t(M1))[, iu]) %*% S1 %*% 
                (t(M1))[, iv]
            G1[iv, iu] <- G1[iu, iv]
            G2[iu, iv] <- G2[iu, iv] + t((t(M2))[, iu]) %*% S2 %*% 
                (t(M2))[, iv]
            G2[iv, iu] <- G2[iu, iv]
        }
    }
    G1 <- G1/n1/Enorm(sumx1/n1)^2
    G2 <- G2/n2/Enorm(sumx2/n2)^2
 #   eva1 <- eigen(G1, symmetric = TRUE, EISPACK = TRUE)
 eva1 <- eigen(G1, symmetric = TRUE)
    pcar1 <- eva1$vectors[, 1:p]
    pcasd1 <- sqrt(abs(eva1$values[1:p]))
#    eva2 <- eigen(G2, symmetric = TRUE, EISPACK = TRUE)
 eva2 <- eigen(G2, symmetric = TRUE)
    pcar2 <- eva2$vectors[, 1:p]
    pcasd2 <- sqrt(abs(eva2$values[1:p]))
    if ((pcasd1[p] < 1e-06) || (pcasd2[p] < 1e-06)) {
        offset <- 1e-06
        cat("*")
        pcasd1 <- sqrt(pcasd1^2 + offset)
        pcasd2 <- sqrt(pcasd2^2 + offset)
    }
    Ahat1 <- n1 * t(M1) %*% (pcar1 %*% diag(1/pcasd1^2) %*% t(pcar1)) %*% 
        M1
    Ahat2 <- n2 * t(M2) %*% (pcar2 %*% diag(1/pcasd2^2) %*% t(pcar2)) %*% 
        M2


    Ahat <- ( Ahat1 +  Ahat2)
#    eva <- eigen(Ahat, symmetric = TRUE, EISPACK = TRUE)
eva <- eigen(Ahat, symmetric = TRUE)



    lambdamin <- eva$values[p+1]
    pval <- 1 - pchisq(lambdamin, p )
#print(lambdamin)
#print(pval)
    z <- list()
    z$pval <- pval
    z$df <- p
    z$lambdamin <- lambdamin
    return(z)
}


Goodall<-function( pool , n1, n2, p=0){
   tan1<- pool$tan[,1:n1]
   tan2<- pool$tan[,(n1+1):(n1+n2)]
   kt <- dim(tan1)[1]
   n <- n1+n2
   k<-pool$k
   m<-pool$m
   if (p == 0){
    p <- min(k * m - (m * (m - 1))/2 - 1 - m, n1 + n2 - 2)
}
    top <-  Enorm( apply(tan1,1,mean) - apply(tan2,1,mean) )**2
    bot <- sum(diag(var(t(tan1))))* (n1-1) +  sum(diag(var(t(tan2)))) * (n2-1)
    Fstat <- ((n1 + n2 - 2)/(1/n1 + 1/n2) * top)/bot
    pval <- 1 - pf(Fstat, p, (n1 + n2 - 2) * p)
    z<-list()
    z$F <- Fstat
    z$pval <- pval
    z$df1 <- p
    z$df2 <- (n1 + n2 - 2) * p
    return(z)
}

Hotelling<-function( pool , n1, n2, p=0){
   tan1<- pool$tan[,1:n1]
   tan2<-pool$tan[,(n1+1):(n1+n2)]
   kt <- dim(tan1)[1]
   n <- n1+n2
   k<-pool$k
   m<-pool$m
    S1 <- var(t(tan1))
    S2 <- var(t(tan2))
    Sw <- ((n1 - 1) * S1 + (n2 - 1) * S2)/(n1 + n2 - 2)
   if (p == 0){
    p <- min(k * m - (m * (m - 1))/2 - 1 - m, n1 + n2 - 2)
}
#    eva <- eigen(Sw, symmetric = TRUE,EISPACK=TRUE)
 eva <- eigen(Sw, symmetric = TRUE)
    pcar <- eva$vectors[, 1:p]
    pcasd <- sqrt(abs(eva$values[1:p]))
    if (pcasd[p] < 1e-06) {
        offset <- 1e-06
        cat("*")
        pcasd <- sqrt(pcasd^2 + offset)
    }
    lam <- rep(0, times =  kt )
    lam[1:p] <- 1/pcasd^2
    Suinv <- eva$vectors %*% diag(lam) %*% t(eva$vectors)
    pcax <- t(pool$tan) %*% pcar
    one1 <- matrix(1/n1, n1, 1)
    one2 <- matrix(1/n2, n2, 1)
    oneone <- rbind(one1, -one2)
    vbar <- pool$tan %*% oneone
    scores1 <- matrix(vbar, 1, kt ) %*% pcar
    scores <- scores1/pcasd
    F.partition <- ((scores[1:p]^2) * (n1 * n2 * (n1 + n2 - p -
        1)))/((n1 + n2) * (n1 + n2 - 2) * p)
    FF <- sum(F.partition)
    pval <- 1 - pf(FF, p, (n1 + n2 - p - 1))
    z<-list()
    z$F.partition <- F.partition
    z$F <- FF
    z$pval <- pval
    z$df1 <- p
    z$T.df1 <- p
    z$df2 <- (n1 + n2 - p - 1)
    mm <- n - 2
    z$T.df2 <- mm
    z$Tsq <- FF * (n1 + n2) * (n1 + n2 - 2) * p/(n1 * n2)/(n1 +
        n2 - p - 1)
    z$Tsq.partition <- F.partition * (n1 + n2) * (n1 + n2 - 2) *
        p/(n1 * n2)/(n1 + n2 - p - 1)
    return(z)
}

James<-function( pool , n1, n2, p=0, table=FALSE){
   tan1<- pool$tan[,1:n1]
   tan2<-pool$tan[,(n1+1):(n1+n2)]
   kt <- dim(tan1)[1]
   n <- n1+n2
   k<-pool$k
   m<-pool$m
    S1 <- var(t(tan1))
    S2 <- var(t(tan2))
    Sw <- S1/n1 + S2/n2
   if (p == 0){
    p <- min(k * m - (m * (m - 1))/2 - 1 - m, n1 + n2 - 2)
}
#    eva <- eigen(Sw, symmetric = TRUE,EISPACK=TRUE)
 eva <- eigen(Sw, symmetric = TRUE)
    pcar <- eva$vectors[, 1:p]
    pcasd <- sqrt(abs(eva$values[1:p]))
    if (pcasd[p] < 1e-06) {
        offset <- 1e-06
        cat("*")
        pcasd <- sqrt(pcasd^2 + offset)
    }
    lam <- rep(0, times =  kt )
    lam[1:p] <- 1/pcasd^2
    Suinv <- eva$vectors %*% diag(lam) %*% t(eva$vectors)
    pcax <- t(pool$tan) %*% pcar
    one1 <- matrix(1/n1, n1, 1)
    one2 <- matrix(1/n2, n2, 1)
    oneone <- rbind(one1, -one2)
    vbar <- pool$tan %*% oneone
#    scores1 <- matrix(vbar, 1, kt ) %*% pcar
#    scores <- scores1/pcasd
#    F.partition <- ((scores[1:p]^2) * (n1 * n2 * (n1 + n2 - p -
#        1)))/((n1 + n2) * (n1 + n2 - 2) * p)
#    FF <- sum(F.partition)
#    pval <- 1 - pf(FF, p, (n1 + n2 - p - 1))
#########
#        ginvSw<- pcar%*%diag(1/pcasd**2)%*%t(pcar)
        ginvSw <- Suinv
        pval=0
	T1<- sum(diag(( ginvSw%*%S1/n1 )))
	T2<- sum(diag(( ginvSw%*%S2/n2 )))
	T1sq<- sum(diag(( (ginvSw%*%S1/n1) %*% ginvSw%*%S1/n1)))
	T2sq<- sum(diag(( (ginvSw%*%S2/n2) %*% ginvSw%*%S2/n2)))

	Tsq<- (t(vbar)%*%( ginvSw) %*%vbar)[1,1]
if (table==TRUE){
        AA <- 1 + 1/(2*p) * (T1**2/(n1-1) + T2**2/(n2-1) )
	BB <- 1/(p*(p+2)) * (  ( T1**2/2 + T1sq )/(n1-1) + (T2**2/2 + T2sq )/(n2-1) )
	 kk<-rep(0,times=1000)
	for (i in 0:999){
	alphai<-i/1000
	kk[i+1]<- qchisq(alphai,df=p)*(AA + BB * qchisq(alphai,df=p) )
	}
        pval<- 1 - max( c(1:1000)[kk<Tsq] )/1000
}
#######
    z<-list()
    z$pval <- pval
 
    z$Tsq <- Tsq

    return(z)
}


tpsgrid<-function (TT, YY, xbegin = -999, ybegin = -999, xwidth= -999, opt = 1, ext = 0.1, ngrid = 22,cex=1,pch=20,
col=2,zslice = 0, mag=1, axes3=FALSE )
{

    k <- nrow(TT)
    m <- dim(TT)[2]

YY <- TT + (YY-TT)*mag

bb<-array(TT,c(dim(TT),1))
aa<-defplotsize2(bb)

if (xwidth== -999){
xwidth<-aa$width
}
if (xbegin== -999){
xbegin<-aa$xl
}
if (ybegin== -999){
ybegin<-aa$yl
}


if (m == 3){
zup<- max(TT[,3])
zlo<- min(TT[,3])

zpos<-zslice
for (ii in 1:length(zslice) ){
zpos[ii] <- (zup + zlo)/2 + (zup - zlo)/2* zslice[ii]
}
}



    xstart <- xbegin
    ystart <- ybegin
    ngrid <- trunc(ngrid/2) * 2
    kx <- ngrid
    ky <- ngrid - 1
    l <- kx * ky
    step <- xwidth/(kx - 1)
    r <- 0
    X <- rep(0, times = kx)
    Y2 <- rep(0, times = ky)
    for (p in 1:kx) {
        ystart <- ybegin
        xstart <- xstart + step
        for (q in 1:ky) {
            ystart <- ystart + step
            r <- r + 1
            X[r] <- xstart
            Y2[r] <- ystart
        }
    }

    TPS <- bendingenergy(TT)
    gamma11 <- TPS$gamma11
    gamma21 <- TPS$gamma21
    gamma31 <- TPS$gamma31
    W <- gamma11 %*% YY
    ta <- t(gamma21 %*% YY)
    B <- gamma31 %*% YY
    WtY <- t(W) %*% YY
    trace <- c(0)
    for (i in 1:m) {
        trace <- trace + WtY[i, i]
    }
    benergy <- (1/(8 * pi)) * trace
    l <- kx * ky
    phi <- matrix(0, l, m)
    s <- matrix(0, k, 1)


for (islice in 1:length(zslice)){
if (m == 3){  
    refc <- matrix(c(X, Y2, rep(zpos[islice],times=kx*ky) ), kx * ky, m)
}
if (m == 2){
    refc <- matrix(c(X, Y2), kx * ky, m)
}


    for (i in 1:l) {
        s <- matrix(0, k, 1)
        for (im in 1:k) {
            s[im, ] <- sigma(refc[i, ] - TT[im, ])
        }
        phi[i, ] <- ta + t(B) %*% refc[i, ] + t(W) %*% s
    }

if (m == 3){

if (opt==2){
    shapes3d(TT,color=2,axes3=axes3,rglopen=FALSE)
    shapes3d(YY,color=4,rglopen=FALSE)
for (i in 1:k){
lines3d(rbind(TT[i,],YY[i,]),col=1)
}
for (j in 1:kx){
    lines3d(refc[ ( (j-1)*ky+1) : (ky*j) ,],color=6)
}
for (j in 1:ky){
    lines3d(refc[  (0:(kx-1)*ky) + j ,],color=6)
}
}

shapes3d(TT,color=2,axes3=axes3,rglopen=FALSE)
shapes3d(YY,color=4,rglopen=FALSE)
for (i in 1:k){
lines3d(rbind(TT[i,],YY[i,]),col=1)
}
for (j in 1:kx){
    lines3d(phi[ ( (j-1)*ky+1) : (ky*j) ,],color=3)
}
for (j in 1:ky){
    lines3d(phi[  (0:(kx-1)*ky) + j ,],color=3)
}
}
}

if (m == 2){
    par(pty = "s")
    if (opt == 2) {
        par(mfrow=c(1,2))
        order <- linegrid(refc, kx, ky)
        plot(order[1:l, 1], order[1:l, 2], type = "l", xlim = c(xbegin -
            xwidth * ext, xbegin + xwidth * (1 + ext)), ylim = c(ybegin -
            (xwidth * ky)/kx * ext, ybegin + (xwidth * ky)/kx *
            (1 + ext)), xlab = " ", ylab = " ")
        lines(order[(l + 1):(2 * l), 1], order[(l + 1):(2 * l),
            2], type = "l")
        points(TT, cex = cex,pch=pch,col=col)
    }
    order <- linegrid(phi, kx, ky)
    plot(order[1:l, 1], order[1:l, 2], type = "l", xlim = c(xbegin -
        xwidth * ext, xbegin + xwidth * (1 + ext)), ylim = c(ybegin -
        (xwidth * ext * ky)/kx, ybegin + (xwidth * (1 + ext) *
        ky)/kx), xlab = " ", ylab = " ")
    lines(order[(l + 1):(2 * l), 1], order[(l + 1):(2 * l), 2],
        type = "l")
    points(YY, cex = cex,pch=pch,col=col+1)
    points(TT, cex=cex,pch=pch,col=col)
    for (i in 1:(k)){
    arrows( TT[i,1],TT[i,2] , YY[i,1],YY[i,2] ,col=col+2,length=0.1,angle=20)
    }
}


}
#



rotatexyz<-function(x,thetax,thetay,thetaz){
thetax<-thetax/180*pi
thetay<-thetay/180*pi
thetaz<-thetaz/180*pi

Rx<-matrix(c(1,0,0,0,cos(thetax),sin(thetax),0,-sin(thetax),cos(thetax)),3,3)
Ry<-matrix(c(cos(thetay),0,sin(thetay),0,1,0,-sin(thetay),0,cos(thetay)),3,3)
Rz<-matrix(c(cos(thetaz),sin(thetaz),0,-sin(thetaz),cos(thetaz),0,0,0,1),3,3)


y<-x
n<-dim(x)[3]
for (i in 1:n){
y[,,i]<-x[,,i]%*%Rx%*%Ry%*%Rz
}
y
}


rigidbody<-function(X,transx=0, transy=0, transz=0, thetax=0,thetay=0,thetaz=0){
if (is.matrix(X)){
X<-array(X,c(dim(X),1))
}
m<-dim(X)[2]
n<-dim(X)[3]
Y<-X
if (m == 2){
#xx<-as.3d(X)
if(dim(X)[3]<2){xx <- array(as.3d(X),dim=c(nrow(X),3,1))}else{xx <- as.3d(X)}
for (i in 1:n){
for (j in 1:m){
xx[j,,i]<-xx[j,,i]-c(transx,transy,transz)
} }
yy<-rotatexyz(xx,thetax,thetay,thetaz)
Y<-yy
if ( (sum(abs(yy[,3,])))< 0.000000001){
Y<-yy[,1:2,]
}
}
if (m == 3){
for (i in 1:n){
for (j in 1:m){
X[j,,i]<-X[j,,i]-c(transx,transy,transz)
} }
Y<-rotatexyz(X,thetax,thetay,thetaz)
}
Y
}





as.3d<-function(X){
k<-dim(X)[1]

if (is.matrix(X)){
X<-array( X, c(dim(X),1) )
}


n<-dim(X)[3]
if (dim(X)[2] != 2){
print("not 2 dimensional!")
}
Y <- array(0, c(k,3,n) )
Y[,1:2,] <- X

if (n==1){
Y <- Y[,,1]
}

Y
}


abind <- function( X1 , X2 ){
k<-dim(X1)[1]
m<-dim(X1)[2]

if (is.matrix(X1)){
tem<- array(0, c(k,m,1) )
tem[,,1]<-X1
X1<-tem
}

if (is.matrix(X2)){
tem<- array(0, c(k,m,1) )
tem[,,1]<-X2
X2<-tem
}

n1<-dim(X1)[3]
n2<-dim(X2)[3]
Y<-array(0,c(k,m,n1+n2))
Y[,,1:n1]<-X1
Y[,,(n1+1):(n1+n2)]<- X2
Y
}


shapes3d<-function(x,loop=0,type="p",color=2,joinline=c(1:1),axes3=FALSE,rglopen=TRUE){



if (is.matrix(x)){ 
xt <- array( 0, c(dim(x),1) )
xt[,,1]<-x
x<-xt 
}

if (is.array(x)==FALSE){
cat("Data not in right format : require an array \n")
}

if (is.array(x)==TRUE){

if (rglopen){
rgl.open()
rgl.bg(color="white")
}



if (dim(x)[2]==2){
x <- as.3d(x)
}



if (loop == 0){

k<-dim(x)[1]
sz<- centroid.size( x[,,1] )/sqrt(k) /30
plotshapes3d(x,type=type,color=color,size=sz,joinline=joinline)
if (axes3){
axes3d( color= "black")
title3d(xlab="x",ylab="y",zlab="z",color="black")
}
}
if (loop > 0){

for (i in 1:loop){

plotshapestime3d(x,type=type)

}

}
}
}

plotshapes3d<-
function (x, type = "p", rgl = TRUE, color = 2, size = 1, joinline = c(1:1))
{
    k <- dim(x)[1]
    n <- dim(x)[3]
    y <- matrix(0, k * n, 3)
    for (i in 1:n) {
        y[(i - 1) * k + (1:k), ] <- x[, , i]
    }
    if (rgl == FALSE) {
        par(mfrow = c(1, 1))
        out <- defplotsize3(x)
        xl <- out$xl
        xu <- out$xu
        yl <- out$yl
        yu <- out$yu
        zl <- out$zl
        zu <- out$zu
        scatterplot3d(y, xlim = c(xl, xu), ylim = c(yl, yu),
            zlim = c(zl, zu), xlab = "x", ylab = "y", zlab = "z",
            axis = TRUE, type = type, color = color, highlight.3d = TRUE)
    }
    if (rgl == TRUE) {
if (type=="l"){
        points3d(y, col = color, size = size)
            for (j in 1:n) {
                lines3d(x[, , j],col= 8)
            }
}
if (type=="dots"){
        points3d(y, col = color, size = size)
}
        if (type=="p"){
        spheres3d(y, col = color, radius = size)
}
        if (length(joinline) > 1) {
            for (j in 1:n) {
                lines3d(x[joinline, , j],col= 8)
            }
        }

    }
}



plotshapestime3d<-function (x,type="p")
{
par(mfrow=c(1,1))
out<-defplotsize3(x)
xl<-out$xl
xu<-out$xu
yl<-out$yl
yu<-out$yu
zl<-out$zl
zu<-out$zu

n<-dim(x)[3]
for (i in 1:n){
scatterplot3d(x[,,i],xlim=c(xl,xu),ylim=c(yl,yu),zlim=c(zl,zu),xlab="x",ylab="y",zlab="z",
  axis=TRUE,type=type,highlight.3d=TRUE)
title(i)

}






}



plotPDMnoaxis3<-function (mean, pc, sd, xl, xu, yl, yu, lineorder, i)
{
    fig <- mean + i * pc * sd
    k <- length(mean)/2
    figx <- fig[1:k]
    figy <- fig[(k + 1):(2 * k)]
    plot(figx, figy, axes = FALSE, xlab = " ", ylab = " ", ylim = c(yl,
        yu), type = "n", xlim = c(xl, xu))
    text(figx, figy, 1:k)
    lines(figx[lineorder], figy[lineorder])
    for (aa in 1:9999){aaa<-1}
}
#################################


shapepca<-function (proc, pcno = c(1, 2, 3), type = "r", mag = 1, joinline = c(1,1),project=c(1,2),scores3d=FALSE,color=2,axes3=FALSE,rglopen=TRUE,zslice=0)
{
if (scores3d==TRUE){
axes3<-TRUE
sz<-max(proc$rawscores[,max(pcno)]) - min(proc$rawscores[,max(pcno)])
spheres3d( proc$rawscores[,pcno] , radius=sz/30, col=color)
if (axes3){
axes3d()
}
}
        
m <- dim(proc$mshape)[2]
k<-dim(proc$mshape)[1]


if (scores3d==FALSE){
    if ((m == 2)) {
        out <- defplotsize2(proc$rotated,project=project)
        xl <- out$xl
        yl <- out$yl
        width <- out$width
        plotpca(proc, pcno, type, mag, xl, yl, width, joinline,project)
    }
    if ((m == 3)&&(type=="m")) {
#        plot3Dmean(proc)
#        cat("Mean shape \n")
#        for (i in 1:length(pcno)) {
#            cat("PC ", pcno[i], " \n")
#            plot3Dpca(proc, pcno[i])
#        }

        for (i in 1:length(pcno)) {
            cat("PC ", pcno[i], " \n")
plotpca3d(proc,pcno[i])
        }

    }


## correct length of tangent vector if in Helmertized space





   h <- defh(k - 1)
    zero <- matrix(0, k - 1, k)
    H <- cbind(h, zero, zero)
    H1 <- cbind(zero, h, zero)
   H2 <- cbind( zero, zero, h)
    H <- rbind(H, H1, H2)

    if (dim(proc$pcar)[1] == (3 * (k - 1))) {
        pcarot <- (t(H) %*% proc$pcar)
proc$pcar<-pcarot    
}


       if (((m == 3)&&(type!="m"))&&(type!="g")) {
if (rglopen){
rgl.open()
rgl.bg(color="white")
}


sz<- centroid.size( proc$mshape )/sqrt(k) /30
        spheres3d(proc$mshape,radius=sz, col=color)
if (axes3){
axes3d()
}

        for (i in pcno ) {
  pc<- proc$mshape + 3*mag*proc$pcasd[i]*cbind(proc$pcar[1:k,i],
           proc$pcar[(k+1):(2*k),i],proc$pcar[(2*k+1):(3*k),i])
for (j in 1:k){
lines3d(rbind(proc$mshape[j,],pc[j,]),col=i)
}
 
        }
    }
}


if ((m == 3)&&(type=="g")){
if (rglopen){
rgl.open()
rgl.bg(color="white")
}

       for (i in pcno ) {
  pc<- proc$mshape + 3*mag*proc$pcasd[i]*cbind(proc$pcar[1:k,i],
           proc$pcar[(k+1):(2*k),i],proc$pcar[(2*k+1):(3*k),i])
tpsgrid( proc$mshape, pc, zslice=zslice)
}
}




}


############

plotpca3d<-function (procreg, pcno=1)
{
par(mfrow=c(1,1))
out<-defplotsize3(procreg$rotated)
xl<-out$xl
xu<-out$xu
yl<-out$yl
yu<-out$yu
zl<-out$zl
zu<-out$zu

    k <- dim(procreg$mshape)[1]
subx<-1:k
suby<-(k+1):(2*k)
subz<-(2*k+1):(3*k)

evec<-cbind(procreg$pcar[subx,pcno],procreg$pcar[suby,pcno],procreg$pcar[subz,pcno])


    for (j in 1:10) {
        for (ii in -12:12) {
            mag <- ii/4
scatterplot3d(procreg$mshape+mag*evec*procreg$pcasd[pcno],xlim=c(xl,xu),ylim=c(yl,yu),zlim=c(zl,zu),xlab="x",ylab="y",zlab="z",
  axis=TRUE,highlight.3d=TRUE)
title(pcno)
}
      for (ii in -11:11) {
            mag <- -ii/4
scatterplot3d(procreg$mshape+mag*evec*procreg$pcasd[pcno],xlim=c(xl,xu),ylim=c(yl,yu),zlim=c(zl,zu),xlab="x",ylab="y",zlab="z",
  axis=TRUE)
title(pcno)
        }
    }
}


##############################
Hotelling2Djames<-
function (A, B)
{
    z <- list(Tsq = 0, pval = 0)
    n1 <- dim(A)[3]
    n2 <- dim(B)[3]
    n <- n1 + n2
    k <- dim(A)[1]
    m <- dim(B)[2]
    if (m != 2) {
        print("Data not two dimensional")
        return(z)
    }
    else {
        pool <- array(0, c(k, m, n))
        pool[, , 1:n1] <- A
        pool[, , (n1 + 1):n] <- B
        poolpr <- procrustes2d(pool, 1, 2)
        S1 <- var(t(poolpr$tan[, 1:n1]))
        S2 <- var(t(poolpr$tan[, (n1 + 1):(n1 + n2)]))
        gamma <- realtocomplex(preshape(poolpr$mshape))
        Sw <- ( S1/n1 + S2/n2)
        p <- 2 * k - 4
#	TT<-eigen(Sw,symmetric=TRUE,EISPACK=TRUE)
TT<-eigen(Sw,symmetric=TRUE)
        pcar <- TT$vectors[, 1:p]
        pcasd <- sqrt(abs(TT$values[1:p]))
####### add small offset if defecient in rank
if (pcasd[p] < 0.000001)
{offset<- 0.000001
cat("*")
pcasd <- sqrt(pcasd**2 + offset**2)
}
#######################################
        pcax <- t(poolpr$tan) %*% pcar
        h <- defh(k - 1)
        zero <- matrix(0, k - 1, k)
        H <- cbind(h, zero)
        H1 <- cbind(zero, h)
        H <- rbind(H, H1)
        meanxy <- t(H) %*% V(gamma)
        realrot <- t(H) %*% pcar
        one1 <- matrix(1/n1, n1, 1)
        one2 <- matrix(1/n2, n2, 1)
        oneone <- rbind(one1, -one2)
        vbar <- poolpr$tan %*% oneone
        scores1 <- matrix(vbar, 1, (2 * k - 2)) %*% pcar
        scores <- scores1/pcasd
        F.partition <- ((scores[1:p]^2) * (n1 * n2 * (n1 + n2 -
            p - 1)))/((n1 + n2) * (n1 + n2 - 2) * p)
        FF <- sum(F.partition)
        pval <- 1 - pf(FF, p, (n1 + n2 - p - 1))
        ginvSw<- pcar%*%diag(1/pcasd**2)%*%t(pcar)
	T1<- sum(diag(( ginvSw%*%S1/n1 )))
	T2<- sum(diag(( ginvSw%*%S2/n2 )))
	T1sq<- sum(diag(( (ginvSw%*%S1/n1) %*% ginvSw%*%S1/n1)))
	T2sq<- sum(diag(( (ginvSw%*%S2/n2) %*% ginvSw%*%S2/n2)))

	Tsq<- (t(vbar)%*%( ginvSw) %*%vbar)[1,1]
        AA <- 1 + 1/(2*p) * (T1**2/(n1-1) + T2**2/(n2-1) )
	BB <- 1/(p*(p+2)) * (  ( T1**2/2 + T1sq )/(n1-1) + (T2**2/2 + T2sq )/(n2-1) )

	 kk<-rep(0,times=1000)
	for (i in 0:999){
	alphai<-i/1000
	kk[i+1]<- qchisq(alphai,df=p)*(AA + BB * qchisq(alphai,df=p) )
	}
        pval<- 1 - max( c(1:1000)[kk<Tsq] )/1000

#        z$F.partition <- F.partition
#        z$F <- FF
        z$pval <- pval
        z$Tsq <- Tsq

#        z$df1 <- p
#        z$T.df1 <- p
#        z$df2 <- (n1 + n2 - p - 1)
#        mm <- n - 2
#        z$T.df2 <- mm
#        z$Tsq <- FF * (n1 + n2) * (n1 + n2 - 2) * p/(n1 * n2)/(n1 +
#            n2 - p - 1)
#        z$Tsq.partition <- F.partition * (n1 + n2) * (n1 + n2 -
#            2) * p/(n1 * n2)/(n1 + n2 - p - 1)
        return(z)
    }
}





MGM<-function(zst){
nsam<-dim(zst)[2]
k<-dim(zst)[1]
Mhat<-matrix(0,k-1,k-2)
lamhat<-rep(0,times=(k-1))
Sighat<-matrix(0,k-2,k-2)
kk<-k*2-2
t1 <- reassqpr(preshape(zst))/nsam
#	t2 <- eigen(t1,symmetric=TRUE,EISPACK=TRUE)
t2 <- eigen(t1,symmetric=TRUE)
reagamma <- (t2$vectors[, 1]+t2$vectors[,2])/sqrt(2)
	gamma <- Vinv(reagamma)
muhat<-gamma

for (i in 1:(k-2)){
Mhat[,i]<-Vinv( t2$vectors[, 1+(2*i)])
}
for (i in 1:(k-1)){
lamhat[i]<-t2$values[(2*i)-1]
}
for (j in 2:(k-1)){
for (l in 2:(k-1)){
sum<-0
for (i in 1:nsam){
zi<-preshape(zst[,i])
sum<-sum+ st(Mhat[,j-1])%*%zi*st(zi)%*%Mhat[,l-1]*st(zi)%*%muhat*st(muhat)%*%zi
}
Sighat[j-1,l-1]<-1/(lamhat[1]-lamhat[j])/(lamhat[1]-lamhat[l])*sum/nsam
}
}
SR<-Re(Sighat)
SI<-Im(Sighat)
S1<-cbind(SR,SI)
S2<-cbind(-(SI),SR)
S<-rbind(S1,S2)

offset<-0
#es<-eigen(S,symmetric=TRUE,EISPACK=TRUE)$values
es<-eigen(S,symmetric=TRUE)$values
nn<-length(es)
if (es[nn] < 0.000001)
{offset<- 0.000001
#cat("Warning: test: small samples, lambda I added to within group covariance matrix \n")
cat("*")
}
invS<-solve(S+offset*diag(nn))

invSR<-invS[1:(k-2),1:(k-2)]
invSI<-invS[1:(k-2),(k-1):(2*k-4)]
invS<-invSR+1i*invSI
Mhat<-st(Mhat)
MGM<-st(Mhat)%*%invS%*%Mhat
MGM
}




resampletest<-function(A,B,resamples=200,replace=TRUE){  

A1<-A
A2<-B
B<-resamples

k<-dim(A1)[1]
m<-dim(A1)[2]
nmin<- min(dim(A1)[3],dim(A2)[3])
ntot<- dim(A1)[3] + dim(A2)[3]
M<- (k-1)*m - m*(m-1)/2 - 1
if (M >= ntot) {
cat("Warning: Low sample size (n1 + n2 <= p) \n")
}
if ((M >= nmin)&&(replace==TRUE)){
cat("Warning: Low sample sizes : min(n1,n2)<=p : * indicates some regularization carried out \n")
}

permutation<- !replace



   if (is.complex(A1)) {
        tem <- array(0, c(nrow(A1), 2, ncol(A1)))
        tem[, 1, ] <- Re(A1)
        tem[, 2, ] <- Im(A1)
        A1 <- tem
    }
    if (is.complex(A2)) {
        tem <- array(0, c(nrow(A2), 2, ncol(A2)))
        tem[, 1, ] <- Re(A2)
        tem[, 2, ] <- Im(A2)
        A2 <- tem
    }
    m <- dim(A1)[2]
        if (m != 2) {
        print("Data not two dimensional")
        print("Carrying out tests on Procrustes residuals")
        out <- testmeanshapes( A1, A2, resamples = resamples, replace=replace)
        return(out)
    }

    zst1<-A1[,1,]+1i*A1[,2,]
    zst2<-A2[,1,]+1i*A2[,2,]




nsam1<-dim(zst1)[2]
nsam2<-dim(zst2)[2]
k<-dim(zst1)[1]
LL<-(MGM(zst1)+MGM(zst2))*(nsam1+nsam2)
LL1<-cbind(Re(LL),Im(LL))
LL2<-cbind(-Im(LL),Re(LL))
LL<-rbind(LL1,LL2)
#Tumc<-min(eigen(LL,symmetric=TRUE,only.values=TRUE,EISPACK=TRUE)$values)
Tumc<-min(eigen(LL,symmetric=TRUE,only.values=TRUE)$values)


m1<-preshape(procrustes2d(zst1)$mshape)
m1<-m1[,1]+1i*m1[,2]
m2<-preshape(procrustes2d(zst2)$mshape)
m2<-m2[,1]+1i*m2[,2]
m0<-preshape(procrustes2d(cbind(zst1,zst2))$mshape)
m0<-m0[,1]+1i*m0[,2]
d<-length(m1)
H<-defh(k-1)
b<-m1
a<-m0
bt<-b*c((st(b)%*%a)/Mod(st(b)%*%a))
abt<-c(Re(st(bt)%*%a))
ct<-(bt-a*abt)
ct<-ct/sqrt(st(ct)%*%ct)
At<-a%*%st(ct)-ct%*%st(a)
salph<-sqrt(1-abt**2)
calph<-abt
Id<-diag(rep(1,times=d))
U1<-Id+salph*At+(calph-1)*(a%*%st(a)+ct%*%st(ct))
b<-m2
a<-m0
bt<-b*c((st(b)%*%a)/Mod(st(b)%*%a))
abt<-c(Re(st(bt)%*%a))
ct<-(bt-a*abt)
ct<-ct/sqrt(st(ct)%*%ct)
At<-a%*%st(ct)-ct%*%st(a)
salph<-sqrt(1-abt**2)
calph<-abt
Id<-diag(rep(1,times=d))
U2<-Id+salph*At+(calph-1)*(a%*%st(a)+ct%*%st(ct))
yst1<-t(H)%*%U1%*%preshape(zst1)
yst2<-t(H)%*%U2%*%preshape(zst2)
ybind<-cbind(yst1,yst2)

zr1<-array(0,c(k,2,nsam1))
zr2<-array(0,c(k,2,nsam2))
zr3<-array(0,c(k,2,nsam1+nsam2))
zr1[,1,]<-Re(zst1)
zr1[,2,]<-Im(zst1)
zr2[,1,]<-Re(zst2)
zr2[,2,]<-Im(zst2)
zr3[,1,]<-cbind(Re(zst1),Re(zst2))
zr3[,2,]<-cbind(Im(zst1),Im(zst2))
yr1<-array(0,c(k,2,nsam1))
yr2<-array(0,c(k,2,nsam2))
yr3<-array(0,c(k,2,nsam1+nsam2))
yr1[,1,]<-Re(yst1)
yr1[,2,]<-Im(yst1)
yr2[,1,]<-Re(yst2)
yr2[,2,]<-Im(yst2)
yr3[,1,]<-cbind(Re(yst1),Re(yst2))
yr3[,2,]<-cbind(Im(yst1),Im(yst2))


Gtem<-Goodall2D(zr1,zr2)
Htem<-Hotelling2D(zr1,zr2)
Jtem<-Hotelling2Djames(zr1,zr2)

Gumc<-Gtem$F
Humc<-Htem$F
Jumc<-Jtem$Tsq

Gtabpval<-Gtem$pval
Htabpval<-Htem$pval
Jtabpval<-Jtem$pval


if (B > 0){

Tu<-rep(0,times=B)
Gu<-Tu
Hu<-Tu
Ju<-Tu


cat("Resampling...")
cat(c("No of resamples = ",B,"\n"))

if (permutation){
cat("Permutations - sampling without replacement \n")
}
if (permutation==FALSE){
cat("Bootstrap - sampling with replacement \n")
}

for (i in 1:B){
cat(c(i," "))


select1<-sample(1:nsam1,replace=TRUE)
select2<-sample(1:nsam2,replace=TRUE)
zb1<-yst1[,select1]
zb2<-yst2[,select2]
zbgh1<-yr1[,,select1]
zbgh2<-yr2[,,select2]

if (permutation){
select0<-sample(c(1:(nsam1+nsam2)),(nsam1+nsam2),replace=FALSE)
select1<-select0[1:nsam1]
select2<-select0[(nsam1+1):(nsam1+nsam2)]
zb1<-zr3[,1,select1]+1i*zr3[,2,select1]
zb2<-zr3[,1,select2]+1i*zr3[,2,select2]
zbgh1<-yr3[,,select1]
zbgh2<-yr3[,,select2]
}


LL<-(MGM(zb1)+MGM(zb2))*(nsam1+nsam2)
LL1<-cbind(Re(LL),Im(LL))
LL2<-cbind(-Im(LL),Re(LL))
LL<-rbind(LL1,LL2)
#lmin<-min(eigen(LL,symmetric=TRUE,only.values=TRUE,EISPACK=TRUE)$values)
lmin<-min(eigen(LL,symmetric=TRUE,only.values=TRUE)$values)
Tu[i]<-lmin
Gu[i]<-Goodall2D(zbgh1,zbgh2)$F
Hu[i]<-Hotelling2D(zbgh1,zbgh2)$F
Ju[i]<-Hotelling2Djames(zbgh1,zbgh2)$Tsq
}

Tu<-sort(Tu)
numbig<-length(Tu[Tumc<Tu])
pvalb<-(1+numbig)/(B+1)
Gu<-sort(Gu)
numbig<-length(Gu[Gumc<Gu])
pvalG<-(1+numbig)/(B+1)
Hu<-sort(Hu)
numbig<-length(Hu[Humc<Hu])
pvalH<-(1+numbig)/(B+1)
Ju<-sort(Ju)
numbig<-length(Ju[Jumc<Ju])
pvalJ<-(1+numbig)/(B+1)




cat(" \n")
out<-list(lambda=0,lambda.pvalue=0,lambda.table.pvalue=0,H=0,H.pvalue=0,H.table.pvalue=0,J=0,J.pvalue=0,J.table.pvalue=0,
G=0,G.pvalue=0,G.table.pvalue=0)

out$lambda<-Tumc
out$lambda.pvalue<-pvalb
out$lambda.table.pvalue<-1-pchisq(Tumc,2*k-4)
out$H<-Humc
out$H.pvalue<-pvalH
out$H.table.pvalue<-Htabpval
out$J<-Jumc
out$J.pvalue<-pvalJ
out$J.table.pvalue<-Jtabpval
out$G<-Gumc
out$G.pvalue<-pvalG
out$G.table.pvalue<-Gtabpval

}

if (resamples==0){
out<-list(lambda=0,lambda.table.pvalue=0,H=0,H.table.pvalue=0,J=0,J.table.pvalue=0,G=0,G.table.pvalue=0)

out$lambda<-Tumc
out$lambda.table.pvalue<-1-pchisq(Tumc,2*k-4)
out$H<-Humc
out$H.table.pvalue<-Htabpval
out$J<-Jumc
out$J.table.pvalue<-Jtabpval
out$G<-Gumc
out$G.table.pvalue<-Gtabpval


}



out
}



prcomp1<-function (x, retx = TRUE, center = TRUE, scale. = FALSE, tol = NULL, svd=TRUE) 
{
    x <- as.matrix(x)
    x <- scale(x, center = center, scale = scale.)
if (svd==FALSE){
a<-eigen(cov(x))
r<-list(sdev=0,rotation=0,x=0)
r$sdev<-sqrt(abs(a$values))
r$rotation<-a$vectors
r$x<-x%*%a$vectors
}
else
{
    s <- svd(x,nu=0)
  if (!is.null(tol)) {
        rank <- sum(s$d > (s$d[1] * tol))
        if (rank < ncol(x)) 
            s$v <- s$v[, 1:rank, drop = FALSE]
    }
    s$d <- s$d/sqrt(max(1, nrow(x) - 1))
    dimnames(s$v) <- list(colnames(x), paste("PC", seq(len = ncol(s$v)), 
        sep = ""))
    r <- list(sdev = s$d, rotation = s$v)
    if (retx) 
        r$x <- x %*% s$v
    class(r) <- "prcomp1"
}
    r
}




# if using Splus
# prcomp1<-prcomp

defplotsize3<-function(Y)
{
    out <- list(xl = 0, yl = 0, zl=0, xu = 0, yu = 0, zu=0, width = 0)
    n <- dim(Y)[3]
        xm <- mean(Y[, 1, ])
        ym <- mean(Y[, 2, ])
        zm <- mean(Y[, 3, ])
    x <- Y
    x[, 1, ] <- Y[, 1, ] - xm
    x[, 2, ] <- Y[, 2, ] - ym
    x[, 3, ] <- Y[, 3, ] - zm
mn1<-min(x[,1,])
mn2<-min(x[,2,])
mn3<-min(x[,3,])
mx1<-max(x[,1,])
mx2<-max(x[,2,])
mx3<-max(x[,3,])
    xl <- -max(-mn1,mx1)
    yl <- -max(-mn2,mx2)
    zl <- -max(-mn3,mx3)
    width <- max(-2 * xl, -2 * yl, -2*zl)
    out$xl <- -width/2 * 1.2 + xm
    out$yl <- -width/2 * 1.2 + ym
   out$zl <- -width/2 * 1.2 + zm
     out$xu <- width/2 * 1.2 + xm
    out$yu <- width/2 * 1.2 + ym
    out$zu <- width/2 * 1.2 + zm
    out$width <- width * 1.2
    out
}



procOPA<-function(A,B,scale=TRUE,reflect=FALSE){
out<-list(R=0,s=0,Ahat=0,Bhat=0,OSS=0,rmsd=0)
if (is.complex(sum(A))==TRUE){
k<-length(A)
Areal<-matrix(0,k,2)
Areal[,1]<-Re(A)
Areal[,2]<-Im(A)
A<-Areal
}
if (is.complex(sum(B))==TRUE){
k<-length(B)
Breal<-matrix(0,k,2)
Breal[,1]<-Re(B)
Breal[,2]<-Im(B)
B<-Breal
}
k<-dim(A)[1]
if (reflect==FALSE){
R<-fort.ROTATION(A,B)} else 
{
R<-fort.ROTATEANDREFLECT(A,B)
}
s<-1
if (scale==TRUE){
s<-fos(A,B)

if (reflect==TRUE){
s<-fos.REFLECT(A,B)
}

}
Ahat<-fcnt(A)
Bhat<-fcnt(B)%*%R*s
resid<-Ahat-Bhat
OSS<-sum(diag(t(resid)%*%resid))
out$R<-R
out$s<-s
out$Ahat<-Ahat
out$Bhat<-Bhat
m<-dim(Ahat)[2]
out$OSS<-OSS
out$rmsd <- sqrt(OSS/(k))
out
}
defplotsize2<-function(Y,project=c(1,2))
{
    out <- list(xl = 0, yl = 0, xu = 0, yu = 0, width = 0)
    n <- dim(Y)[3]
        xm <- mean(Y[, project[1], ])
        ym <- mean(Y[, project[2], ])
    x <- Y
    x[, project[1], ] <- Y[, project[1], ] - xm
    x[, project[2], ] <- Y[, project[2], ] - ym
    out <- list(xl = 0, yl = 0, width = 0)
mn1<-min(x[,project[1],])
mn2<-min(x[,project[2],])
mx1<-max(x[,project[1],])
mx2<-max(x[,project[2],])
    xl <- -max(-mn1, mx1)
    yl <- -max(-mn2,mx2)
    width <- max(-2 * xl, -2 * yl)
    out$xl <- -width/2 * 1.2 + xm
    out$yl <- -width/2 * 1.2 + ym
    out$xu <- width/2 * 1.2 + xm
    out$yu <- width/2 * 1.2 + ym
    out$width <- width * 1.2
    out
}


plotshapes<-function(A,B=0,joinline=c(1,1),orthproj=c(1,2),color=1,symbol=1){

CHECKOK<-TRUE
if (is.array(A)==FALSE){
if (is.matrix(A)==FALSE){
cat("Error !! argument should be an array or matrix \n")
CHECKOK<-FALSE
}
}

if (CHECKOK){
k<-dim(A)[1]
m<-dim(A)[2]
kk<-k
if (k >= 15){kk<-1}
par(pty="s")
#if (length(c(B))==1){
#par(mfrow=c(1,1))
#}
if (length(c(B))!=1){
par(mfrow=c(1,2))
}
if (length(dim(A))==3){
A<-A[,orthproj,]
}
if (is.matrix(A)==TRUE){
a<-array(0,c(k,2,1))
a[,,1]<-A[,orthproj]
A<-a
}
out<-defplotsize2(A)
width<-out$width
if (length(c(B))!=1){
if (length(dim(B))==3){ 
B<-B[,orthproj,]
}
if (is.matrix(B)==TRUE){
a<-array(0,c(k,2,1))
a[,,1]<-B[,orthproj]
B<-a
}
ans<-defplotsize2(B)
width<-max(out$width,ans$width)
} 
n<-dim(A)[3]

lc<-length(color)
lt<- k*m*n/lc
color<-rep(color,times=lt)
lc<-length(symbol)
lt<- k*m*n/lc
symbol<-rep(symbol,times=lt)



plot(A[,,1],xlim=c(out$xl,out$xl+width),ylim=c(out$yl,
out$yl+width),type="n",xlab=" ",ylab=" ")
for (i in 1:n){
select<- ((i-1)*k*m +1):(i*k*m)
points(A[,,i],pch=symbol[select],col=color[select])
lines(A[joinline,,i])
}
if (length(c(B))!=1){
A<-B
if (is.matrix(A)==TRUE){
a<-array(0,c(k,2,1))
a[,,1]<-A
A<-a
}
out<-defplotsize2(A)
n<-dim(A)[3]
plot(A[,,1],xlim=c(ans$xl,ans$xl+width),ylim=c(ans$yl,
ans$yl+width),type="n",xlab=" ",ylab=" ")
for (i in 1:n){
points(A[,,i],pch=symbol[select],col=color[select])
lines(A[joinline,,i])
}
}
}
}






#
#
#
#
BoxM<-function(A, B, npc)
{
#carries out Box's M test
#(see Mardia, Kent, Bibby 1979, p140)
#in: data arrays A, B
#out: z$M   M statistic
#     z$df degrees of freedom for approx distn of chi-squared statistic
#    z$pval  p-value
	z <- list(M = 0, df = 0, pval = 0)
	n1 <- dim(A)[3]
	n2 <- dim(B)[3]
	k <- dim(A)[1]
	m <- dim(A)[2]
if (m > 2){
print("Only works for 2D data at the moment!")
}
if (m == 2){
	C <- array(0, c(k, m, n1 + n2))
	C[,  , 1:n1] <- A
	C[,  , (n1 + 1):(n1 + n2)] <- B
	Cpr <- procrustes2d(C, 1, 2)
	p <- npc
	ng <- 2
	n <- n1 + n2
	S1 <- var(t(Cpr$tan[1:npc, 1:n1]))
	S2 <- var(t(Cpr$tan[1:npc, (n1 + 1):(n1 + n2)]))
	Su <- ((n1 - 1) * S1 + (n2 - 1) * S2)/(n1 + n2 - 2)
	S1inv <- eigen(S1)$vectors %*% diag(1/eigen(S1)$values) %*% t(eigen(S1)$
		vectors)
	S2inv <- eigen(S2)$vectors %*% diag(1/eigen(S2)$values) %*% t(eigen(S2)$
		vectors)
	logdet1 <- sum(log(eigen(S1inv %*% Su)$values))
	logdet2 <- sum(log(eigen(S2inv %*% Su)$values))
	gam <- 1 - ((2 * p^2 + 3 * p - 1)/(6 * (p + 1) * (ng - 1))) * (1/(n1 - 
		1) + 1/(n2 - 1) - 1/(n - ng))
	M <- gam * ((n1 - 1) * logdet1 + (n2 - 1) * logdet2)
	df <- (p * (p + 1) * (ng - 1))/2
	pval <- 1 - pchisq(M, df)
	z$M <- M
	z$df <- df
	z$pval <- pval 
}
	return(z)
}
Goodall2D<-function(A, B)
{
#Calculates Goodall's two sample F test for 2d data only
#in: data arrays A, B k x 2 x n data arrays
#out: z$F  F statistic
#     z$df1, z$df2  degrees of freedom
#     z$pval: p-value
        z <- list(F = 0, pval = 0, df1 = 0, df2 = 0)
        n1 <- dim(A)[3]
        n2 <- dim(B)[3]
        k <- dim(A)[1]
        m <- dim(A)[2]
        if(m != 2) {
                print("Data not two dimensional")
                return(z)
        }
        p <- 2 * k - 4
        Apr <- procrustes2d(A, 1, 2)
        Bpr <- procrustes2d(B, 1, 2)
        top <- sin(riemdist(Apr$mshape, Bpr$mshape))^2
        bot <- Apr$rmsd1^2 * n1 + Bpr$rmsd1^2 * n2
        Fstat <- ((n1 + n2 - 2)/(1/n1 + 1/n2) * top)/bot
        pval <- 1 - pf(Fstat, p, (n1 + n2 - 2) * p)
        z$F <- Fstat
        z$pval <- pval
        z$df1 <- p
        z$df2 <- (n1 + n2 - 2) * p
        return(z)
}

Goodalltest<-function(A, B,tol1=1e-07,tol2=tol1)
{
#Calculates Goodall's two sample F test
#in: data arrays A, B: 
#out: z$F  F statistic
#     z$df1, z$df2  degrees of freedom
#     z$pval: p-value
        z <- list(F = 0, pval = 0, df1 = 0, df2 = 0)
        n1 <- dim(A)[3]
        n2 <- dim(B)[3]
        k <- dim(A)[1]
        m <- dim(A)[2]
        p <- min(k * m - (m * (m - 1))/2 - 1 - m, n1 + n2 - 2)
        Apr <- procrustesGPA(A,tol1,tol2)
        Bpr <- procrustesGPA(B,tol1,tol2)
        top <- sin(riemdist(Apr$mshape, Bpr$mshape))^2
        bot <- Apr$rmsd1^2 * n1 + Bpr$rmsd1^2 * n2
        Fstat <- ((n1 + n2 - 2)/(1/n1 + 1/n2) * top)/bot
        pval <- 1 - pf(Fstat, p, (n1 + n2 - 2) * p)
        z$F <- Fstat
        z$pval <- pval
        z$df1 <- p
        z$df2 <- (n1 + n2 - 2) * p
        return(z)
}

Hotelling2D<-function (A, B) 
{
    z <- list(Tsq.partition = 0, Tsq = 0, F.partition = 0, F = 0, 
        pval = 0, df1 = 0, df2 = 0, T.df1 = 0, T.df2 = 0)
    n1 <- dim(A)[3]
    n2 <- dim(B)[3]
    n <- n1 + n2
    k <- dim(A)[1]
    m <- dim(B)[2]
    if (m != 2) {
        print("Data not two dimensional")
        return(z)
    }
    else {
        pool <- array(0, c(k, m, n))
        pool[, , 1:n1] <- A
        pool[, , (n1 + 1):n] <- B
        poolpr <- procrustes2d(pool, 1, 2)
        S1 <- var(t(poolpr$tan[, 1:n1]))
        S2 <- var(t(poolpr$tan[, (n1 + 1):(n1 + n2)]))
        gamma <- realtocomplex(preshape(poolpr$mshape))
        Sw <- ((n1 - 1) * S1 + (n2 - 1) * S2)/(n1 + n2 - 2)
        p <- 2 * k - 4
#        pcar <- eigen(Sw,EISPACK=TRUE)$vectors[, 1:p]
        pcar <- eigen(Sw)$vectors[, 1:p]
        pcasd <- sqrt(abs(eigen(Sw)$values[1:p]))
####### add small offset if defecient in rank
if (pcasd[p] < 0.000001)
{offset<- 0.000001
cat("*")
pcasd <- sqrt(pcasd**2 + offset**2)
}
#######################################
        pcax <- t(poolpr$tan) %*% pcar
        h <- defh(k - 1)
        zero <- matrix(0, k - 1, k)
        H <- cbind(h, zero)
        H1 <- cbind(zero, h)
        H <- rbind(H, H1)
        meanxy <- t(H) %*% V(gamma)
        realrot <- t(H) %*% pcar
        one1 <- matrix(1/n1, n1, 1)
        one2 <- matrix(1/n2, n2, 1)
        oneone <- rbind(one1, -one2)
        vbar <- poolpr$tan %*% oneone
        scores1 <- matrix(vbar, 1, (2 * k - 2)) %*% pcar
        scores <- scores1/pcasd
        F.partition <- ((scores[1:p]^2) * (n1 * n2 * (n1 + n2 - 
            p - 1)))/((n1 + n2) * (n1 + n2 - 2) * p)
        FF <- sum(F.partition)
        pval <- 1 - pf(FF, p, (n1 + n2 - p - 1))
        z$F.partition <- F.partition
        z$F <- FF
        z$pval <- pval
        z$df1 <- p
        z$T.df1 <- p
        z$df2 <- (n1 + n2 - p - 1)
        mm <- n - 2
        z$T.df2 <- mm
        z$Tsq <- FF * (n1+n2)*(n1+n2-2)*p/(n1*n2)/(n1+n2-p-1)
    z$Tsq.partition <- F.partition * (n1+n2)*(n1+n2-2)*p/(n1*n2)/(n1+n2-p-1)
    return(z)
    }
}



Hotellingtest<-function (A, B, tol1 = 1e-07, tol2 = 1e-07) 
{
    z <- list(Tsq.partition = 0, Tsq = 0, F.partition = 0, F = 0, 
        pval = 0, df1 = 0, df2 = 0, T.df1 = 0, T.df2 = 0)
    n1 <- dim(A)[3]
    n2 <- dim(B)[3]
    n <- n1 + n2
    k <- dim(A)[1]
    m <- dim(B)[2]
    pool <- array(0, c(k, m, n))
    pool[, , 1:n1] <- A
    pool[, , (n1 + 1):n] <- B
    poolpr <- procrustesGPA(pool, tol1, tol2,approxtangent=FALSE)
    S1 <- var(t(poolpr$tan[, 1:n1]))
    S2 <- var(t(poolpr$tan[, (n1 + 1):(n1 + n2)]))
    Sw <- ((n1 - 1) * S1 + (n2 - 1) * S2)/(n1 + n2 - 2)
    p <- min(k * m - (m * (m - 1))/2 - 1 - m, n1 + n2 - 2)
    eva<-eigen(Sw,symmetric=TRUE)
    pcar <- eva$vectors[, 1:p]
    pcasd <- sqrt(abs(eva$values[1:p]))

####### add small offset if defecient in rank
if (pcasd[p] < 0.000001)
{
offset<- 0.000001
cat("*")
pcasd <- sqrt(pcasd**2 + offset)
}
#######################################

    lam<-rep(0,times=(k*m-m))
    lam[1:p]<-1/pcasd**2
    Suinv<-eva$vectors%*%diag(lam)%*%t(eva$vectors)
#    check <- p
#    for (i in 1:p) {
#        if (pcasd[p + 1 - i] < 1e-04) {
#            check <- p + 1 - i - 1
#        }
#    }
#    p <- check
    pcax <- t(poolpr$tan) %*% pcar
    one1 <- matrix(1/n1, n1, 1)
    one2 <- matrix(1/n2, n2, 1)
    oneone <- rbind(one1, -one2)
    vbar <- poolpr$tan %*% oneone
    scores1 <- matrix(vbar, 1, m * k-m ) %*% pcar
    scores <- scores1/pcasd
#    tem<-c(t(vbar)%*%Suinv%*%vbar)  #(=Dsq)#
    F.partition <- ((scores[1:p]^2) * (n1 * n2 * (n1 + n2 - p - 
        1)))/((n1 + n2) * (n1 + n2 - 2) * p)
    FF <- sum(F.partition)
    pval <- 1 - pf(FF, p, (n1 + n2 - p - 1))
    z$F.partition <- F.partition
    z$F <- FF
    z$pval <- pval
    z$df1 <- p
    z$T.df1 <- p
    z$df2 <- (n1 + n2 - p - 1)
    mm <- n - 2
    z$T.df2 <- mm
    z$Tsq <- FF * (n1+n2)*(n1+n2-2)*p/(n1*n2)/(n1+n2-p-1)
    z$Tsq.partition <- F.partition * (n1+n2)*(n1+n2-2)*p/(n1*n2)/(n1+n2-p-1)
    return(z)
}



# Hotellingtest<-function(A, B, tol1=1e05,tol2=1e05)
# OLD VERSION using $tan rather than $tanpartial
#{
#Calculates two sample Hotelling Tsq test for testing whether 
#mean shapes are equal (m - Dimensions where m >= 2)
#in: A, B the k x m x n arrays of data for each group
#out: z$F : F-statistic
#     z$df1, z$df2 : dgrees of freedom
#     z$pval: pvalue
#        z <- list(Tsq.partition = 0, Tsq = 0, F.partition = 0, F = 0, pval = 0, 
#                df1 = 0, df2 = 0, T.df1 = 0, T.df2 = 0)
#        n1 <- dim(A)[3]
#        n2 <- dim(B)[3]
#        n <- n1 + n2
#        k <- dim(A)[1]
#        m <- dim(B)[2]
#        pool <- array(0, c(k, m, n))
#        pool[,  , 1:n1] <- A
#        pool[,  , (n1 + 1):n] <- B
#        poolpr <- procrustesGPA(pool,tol1,tol2)
#        S1 <- var(t(poolpr$tan[, 1:n1]))
#        S2 <- var(t(poolpr$tan[, (n1 + 1):(n1 + n2)]))
#        Sw <- ((n1 - 1) * S1 + (n2 - 1) * S2)/(n1 + n2 - 2)
#        p <- min(k * m - (m * (m - 1))/2 - 1 - m, n1 + n2 - 2)
#        pcar <- eigen(Sw)$vectors[, 1:p]
#        pcasd <- sqrt(eigen(Sw)$values[1:p])
#        check<-p
## checks to see if rank is reasonable
#        for (i in 1:p){
#        if (pcasd[p+1-i] < 0.0001){
#        check<-p+1-i-1
#        }
#        }
#        p<-check
#        pcax <- t(poolpr$tan) %*% pcar
#        one1 <- matrix(1/n1, n1, 1)
#        one2 <- matrix(1/n2, n2, 1)
#        oneone <- rbind(one1,  - one2)
#        vbar <- poolpr$tan %*% oneone
#        scores1 <- matrix(vbar, 1, m*k) %*% pcar
#        scores <- scores1/pcasd
#        F.partition <- ((scores[1:p]^2) * (n1 * n2 * (n1 + n2 - p - 1)))/((n1 + 
#                n2) * (n1 + n2 - 2) * p)
#        FF <- sum(F.partition)
#        pval <- 1 - pf(FF, p, (n1 + n2 - p - 1))
#        z$F.partition <- F.partition
#        z$F <- FF
#        z$pval <- pval
#        z$df1 <- p
#        z$T.df1 <- p
#        z$df2 <- (n1 + n2 - p - 1)
#        mm <- n - 2
#        z$T.df2 <- mm
#        z$Tsq <- (FF * (mm * p))/(mm - p + 1)
#        z$Tsq.partition <- (F.partition * (mm * p))/(mm - p + 1)
#        return(z)
#}


I2mat<-function(Be)
{
	zero <- rep(0, times = dim(Be)[1]^2)
	zero <- matrix(zero, dim(Be)[1], dim(Be)[2])
	temp <- cbind(Be, zero)
	temp1 <- cbind(zero, Be)
	tem <- rbind(temp, temp1)
	tem
}

tpsgrid.old<-function (TT, YY, xbegin = -999, ybegin = -999, xwidth= -999, opt = 2, ext = 0.1,
    ngrid = 22,cex=1,pch=20,col=2)
{
    k <- nrow(TT)
if (xwidth== -999){
bb<-array(TT,c(dim(TT),1))
aa<-defplotsize2(bb)
xwidth<-aa$width
}
if (xbegin== -999){
bb<-array(TT,c(dim(TT),1))
aa<-defplotsize2(bb)
xbegin<-aa$xl
}
if (ybegin== -999){
bb<-array(TT,c(dim(TT),1))
aa<-defplotsize2(bb)
ybegin<-aa$yl
}
    xstart <- xbegin
    ystart <- ybegin
    ngrid <- trunc(ngrid/2) * 2
    kx <- ngrid
    ky <- ngrid - 1
    l <- kx * ky
    step <- xwidth/(kx - 1)
    r <- 0
    X <- rep(0, times = kx)
    Y2 <- rep(0, times = ky)
    for (p in 1:kx) {
        ystart <- ybegin
        xstart <- xstart + step
        for (q in 1:ky) {
            ystart <- ystart + step
            r <- r + 1
            X[r] <- xstart
            Y2[r] <- ystart
        }
    }
    refc <- matrix(c(X, Y2), kx * ky, 2)
    TPS <- bendingenergy(TT)
    gamma11 <- TPS$gamma11
    gamma21 <- TPS$gamma21
    gamma31 <- TPS$gamma31
    W <- gamma11 %*% YY
    ta <- t(gamma21 %*% YY)
    B <- gamma31 %*% YY
    WtY <- t(W) %*% YY
    trace <- c(0)
    for (i in 1:2) {
        trace <- trace + WtY[i, i]
    }
    benergy <- (1/(8 * pi)) * trace
    l <- kx * ky
    phi <- matrix(0, l, 2)
    s <- matrix(0, k, 1)
    for (i in 1:l) {
        s <- matrix(0, k, 1)
        for (m in 1:k) {
            s[m, ] <- sigma(refc[i, ] - TT[m, ])
        }
        phi[i, ] <- ta + t(B) %*% refc[i, ] + t(W) %*% s
    }
    par(pty = "s")
    if (opt == 2) {
        par(mfrow=c(1,2))
        order <- linegrid(refc, kx, ky)
        plot(order[1:l, 1], order[1:l, 2], type = "l", xlim = c(xbegin -
            xwidth * ext, xbegin + xwidth * (1 + ext)), ylim = c(ybegin -
            (xwidth * ky)/kx * ext, ybegin + (xwidth * ky)/kx *
            (1 + ext)), xlab = " ", ylab = " ")
        lines(order[(l + 1):(2 * l), 1], order[(l + 1):(2 * l),
            2], type = "l")
        points(TT, cex = cex,pch=pch,col=col)
    }
    order <- linegrid(phi, kx, ky)
    plot(order[1:l, 1], order[1:l, 2], type = "l", xlim = c(xbegin -
        xwidth * ext, xbegin + xwidth * (1 + ext)), ylim = c(ybegin -
        (xwidth * ext * ky)/kx, ybegin + (xwidth * (1 + ext) *
        ky)/kx), xlab = " ", ylab = " ")
    lines(order[(l + 1):(2 * l), 1], order[(l + 1):(2 * l), 2],
        type = "l")
    points(YY, cex = cex,pch=pch,col=col)
}
#
#
V<-function(z)
{
#input complex k -vector
#ouput vectorized 2k vector of real stacked on imaginary components
	x <- c(Re(z), Im(z))
	x
}
Vinv<-function(x)
{
#input vectorized 2k vector of x1 stacked on x2 components
#input complex k -vector of the form x1 + 1i*x2
	nel <- length(x)/2
	zx <- x[1:nel]
	zy <- x[(nel + 1):(2 * nel)]
	z <- zx + (1i) * zy
	z
}
Vmat<-function(z)
{
#as Vinv but input is a k x n complex matrix
# output 2k x n matrix of stacked real then complex components
	x <- rbind(Re(z), Im(z))
	x
}

bendingenergy<-function (TT)
{
    z <- list(gamma11 = 0, gamma21 = 0, gamma31 = 0, prinwarps = 0,
        prinwarpeval = 0, Un = 0)


    k <- nrow(TT)
    m <- dim(TT)[2]
    S <- matrix(0, k, k)
    for (i in 1:k) {
        for (j in 1:k) {
            S[i, j] <- sigma(TT[i, ] - TT[j, ])
        }
    }

    one <- matrix(1, k, 1)
    zero <- matrix(0, m+1, m+1)
#    P <- cbind(S, one, TT)
    P <- rbind(S, t(one))
    Q <- rbind(P, t(TT))
    O <- cbind(one, TT)
    U <- rbind(O, zero)
    star <- cbind(Q, U)
    star <- matrix(star, k + m+1, k + m+1)
    A <- eigen(star,symmetric=TRUE)
    deltainv <- diag(1/A$values)
    gamma <- A$vectors
    starinv <- gamma %*% deltainv %*% t(gamma)
    gamma11 <- matrix(0, k, k)
    for (i in 1:k) {
        for (j in 1:k) {
            gamma11[i, j] <- starinv[i, j]
        }
    }
    gamma21 <- matrix(0, 1, k)
    for (i in 1:1) {
        for (j in 1:k) {
            gamma21[i, j] <- starinv[k + 1, j]
        }
    }
    gamma31 <- matrix(0, m, k)
    for (i in 1: (m) ) {
        for (j in 1:k) {
            gamma31[i, j] <- starinv[i + k + 1, j]
        }
    }
    prinwarp <- eigen(gamma11, symmetric = TRUE)
    prinwarps <- prinwarp$vectors
    prinwarpeval <- prinwarp$values

####need to rotate to compute affine components
Rot <- prcomp(TT)$rotation
TT<- TT  %*% Rot
    if (m == 2){
    meanxy <- c(TT[, 1], TT[, 2])
    alpha <- sum(meanxy[1:k]^2)
    beta <- sum(meanxy[(k + 1):(2 * k)]^2)
    u1 <- c(alpha * meanxy[(k + 1):(2 * k)], beta * meanxy[1:k])
    u2 <- c(-beta * meanxy[1:k], alpha * meanxy[(k + 1):(2 *
        k)])
    u1 <- u1/sqrt(alpha * beta)/sqrt(alpha + beta)
    u2 <- u2/sqrt(alpha * beta)/sqrt(alpha + beta)


    Un <- matrix(0, 2 * k, 2)
    Un[, 1] <- u1
    Un[, 2] <- u2
    Vn<-Un
    Vn[,1] <-  cbind(Un[1:k,1],Un[(k+1):(2*k),1]) %*%t(Rot)
    Vn[,2] <-  cbind(Un[1:k,2],Un[(k+1):(2*k),2]) %*%t(Rot)
    Un <- Vn
 }

    if (m == 3){
    meanxy <- c(TT[, 1], TT[, 2], TT[,3])
    alpha <- sum(meanxy[1:k]^2)
    beta <- sum(meanxy[(k + 1):(2 * k)]^2)
    gamma <- sum(meanxy[(2*k + 1):(3 * k)]^2)

mu<- meanxy[1:k]
nu <- meanxy[(k+1):(2*k)]
omega<- meanxy[(2*k+1):(3*k)]
ze <- rep(0,times=k)
    u1 <- c( ze , alpha*beta* omega , alpha*gamma* nu) /
                   sqrt(alpha^2*beta^2*gamma + alpha^2*gamma^2*beta)
    u2 <- c( alpha*beta* omega , ze,  beta*gamma* mu ) /
                   sqrt(beta^2*alpha^2*gamma + beta^2*gamma^2*alpha)
    u3 <- c( alpha*gamma*nu , beta*gamma*mu, ze  )/
                   sqrt(alpha^2*gamma^2*beta  + beta^2*gamma^2*alpha)

    u4 <- c( ze , ze , omega) /
                   sqrt( gamma)

    u5 <- c( -beta*gamma*mu , alpha*gamma*nu, ze  )/
                   sqrt(alpha*gamma^2*beta^2  + beta*gamma^2*alpha^2)
    
    tem <- c( -gamma*beta* mu , ze,  beta*alpha* omega ) /
                   sqrt(beta^2*alpha*gamma^2 + beta^2*gamma*alpha^2) 
    tem2 <- tem - u5*sum(u5*tem)
    u4 <- tem2/Enorm(tem2) 

     

    Un <- matrix(0, 3 * k, 5)
    Un[, 1] <- u1
    Un[, 2] <- u2
    Un[,3] <- u3
    Un[,4] <- u4
    Un[,5] <- u5
    Vn<-Un
    Vn[,1] <-  cbind(Un[1:k,1],Un[(k+1):(2*k),1],Un[(2*k+1):(3*k),1]) %*%t(Rot)
    Vn[,2] <-  cbind(Un[1:k,2],Un[(k+1):(2*k),2],Un[(2*k+1):(3*k),2]) %*%t(Rot)
    Vn[,3] <-  cbind(Un[1:k,3],Un[(k+1):(2*k),3],Un[(2*k+1):(3*k),3]) %*%t(Rot)
    Vn[,4] <-  cbind(Un[1:k,4],Un[(k+1):(2*k),4],Un[(2*k+1):(3*k),4]) %*%t(Rot)
    Vn[,5] <-  cbind(Un[1:k,5],Un[(k+1):(2*k),5],Un[(2*k+1):(3*k),5]) %*%t(Rot)
    Un <- Vn
 }


    z$gamma11 <- gamma11
    z$gamma21 <- gamma21
    z$gamma31 <- gamma31
    z$prinwarps <- prinwarps
    z$prinwarpeval <- prinwarpeval
    z$Un <- Un 
    return(z)
}

shaperw <- function( proc ,alpha=1,affine=FALSE){


rw<-proc

if ((alpha != 0)||(affine == TRUE)) {


k<-dim(proc$mshape)[1]
m<-dim(proc$mshape)[2]
n<-dim(proc$mshape)[3]

if (dim(proc$tan)[1]== (k*m-m) ){

if (m == 2){
He <- t( defh(k-1) )
Ze<-He*0
HH <- rbind( cbind( He, Ze ) , cbind( Ze, He ) )
proc$tan <- HH %*% proc$tan
}
if (m == 3){
He <- t( defh(k-1) )
Ze <- He*0
HH <- rbind( cbind( defh(k-1), Ze, Ze ) , cbind(Ze, defh(k-1) , Ze ) , cbind( Ze, Ze, defh(k-1) ) )
proc$tan <- HH %*% proc$tan
}
}




nconstr <- m + m*(m-1)/2 + 1
M <- k*m - nconstr

if (m == 2){
bb<- bendingenergy( proc$mshape )
Gamma11<- bb$gamma11
Be <- rbind( cbind( Gamma11, Gamma11*0 ) , cbind(Gamma11*0, Gamma11) )
Un <- bb$Un
Bedim<- 2
}
if (m == 3){
bb<- bendingenergy( proc$mshape )
Gamma11<- bb$gamma11
Ze<-Gamma11*0
Be <- rbind( cbind( Gamma11, Ze,Ze) , cbind(Ze, Gamma11,Ze) , cbind(Ze,Ze,Gamma11)  )
Un <- bb$Un
Bedim <- 5
}

ev <- eigen(Be,symmetric=TRUE)

Beminusalpha <- ev$vectors %*% diag( c(ev$values[ 1:(M-Bedim)]**(-alpha/2), rep(0, times = nconstr + Bedim)) )%*%t( ev$vectors)
Bealpha <- ev$vectors %*% diag( c(ev$values[ 1:(M-Bedim) ]**(alpha/2),  rep(0, times = nconstr + Bedim)  ))%*%
             t( ev$vectors)

evbe<-ev

SS <- Beminusalpha %*% var(t(proc$tan)) %*% Beminusalpha

ev <- eigen(SS)


relw.vec<- ev$vectors
relw.sd <- sqrt(abs(ev$values))


# ratio of eigenvalues of warps (quoted in book)
rw$percent<- relw.sd**2/sum(relw.sd**2)*100




sgnchange <- sample(c(-1,1),size= m*k ,replace=TRUE)

rw$pcar<- Bealpha %*% relw.vec%*%diag(sgnchange)
rw$pcasd<- relw.sd

rw$rawscores <-  t(  t( relw.vec ) %*% Beminusalpha %*%proc$tan  )

sd<-sqrt(abs(diag(var((rw$rawscores)))))

rw$scores <- (rw$rawscores) %*% diag(1/sd) 
rw$stdscores<-rw$scores

rw$scores<- rw$rawscores 



## partial warp scores



n<- proc$n
evbend<-eigen(Gamma11,symmetric=TRUE)

partialwarpscores<- array( 0 , c( n , m , k) ) 

for (i in 1:m){
partialwarpscores[,i,] <- t( t(evbend$vectors)%*%proc$rotated[,i,] )
}

rw$principalwarps<- evbe$vectors[,(k-m-1):1]
rw$principalwarps.eigenvalues <- evbe$values[(k-m-1):1]
rw$partialwarpscores <- partialwarpscores[,,(k-m-1):1] 

sumvar <- rep(0, times= (k-m-1) )

for (i in 1:(k-m-1)){
sumvar[i] <- sum(diag( var( partialwarpscores[,,k-m-i] ) ))
}

rw$partialwarps.percent<- sumvar/sum(proc$pcasd**2)*100

}

if (affine==TRUE){
dimun<- dim(Un)[2]
rw$pcar<- Un%*%diag(sgnchange[1:(dimun)])
pcno<-c(1: dimun)
rw$rawscores <- t(Un)%*%proc$tan


sd<-sqrt(abs(diag(var(t(rw$rawscores)))))
rw$pcasd <- sd

rw$percent <- sd**2 / sum( proc$pcasd**2 )*100

rw$scores <- t(rw$rawscores) %*% diag(1/sd) 
rw$rawscores<- t(rw$rawscores)
#######

tem<-prcomp1( (rw$rawscores) )

npc<-0
rw$stdscores<-tem$x
    for (i in 1:length(tem$sdev)) {
        if (tem$sdev[i] > 1e-07) {
            npc <- npc + 1
        }
    }
    for (i in 1:npc) {
        rw$stdscores[, i] <- tem$x[, i]/tem$sdev[i]
    }
rw$pcasd <- tem$sdev
rw$percent <- tem$sdev**2 / sum( proc$pcasd**2 )*100
rw$pcar <- Un %*% tem$rotation
rw$rawscores <- tem$x
rw$scores <- rw$rawscores


}

rw

}


bookstein2d<-function(A,l1=1,l2=2){
#input:  A: k x 2 x n array of 2D data, or  k x n complex matrix 
#l1,l2: baseline choice for sending to (-0.5,0),(0.5,0)
#output: z$bshpv - Bookstein shape variables array (including baseline)
# z$mshape - Bookstein mean shape (including baseline points) 
z <- list(k = 0, n = 0, mshape=0,bshpv=0)
if (is.complex(sum(A))==TRUE){
n<-dim(A)[2]
k<-dim(A)[1]
B<-array(0,c(k,2,n))
B[,1,]<-Re(A)
B[,2,]<-Im(A)
A<-B
}

if (is.matrix(A)==TRUE){
bb<-array(A,c(dim(A),1))
A<-bb
}
k<-dim(A)[1]
m<-2
n<-dim(A)[3]
reorder<-c(l1,l2,c(1:k)[-c(l1,l2)])
A[,,]<-A[reorder,,1:n]
bshpv<-array(0,c(k,m,n))
for (i in 1:n)
{
bshpv[,,i]<-bookstein.shpv(A[,,i])
}
bookmean<-matrix(0,k,m)
for (i in 1:n)
{
bookmean<-bookmean+bshpv[,,i]
}
bookmean<-bookmean/n
bookmean[reorder,] <- bookmean
bshpv[reorder, , ] <- bshpv
glim <- max( - min(bshpv), max(bshpv))
#par(pty="s")
#par(mfrow=c(1,1))
#plot(bshpv[,,1],xlim=c(-glim,glim),ylim=c(-glim,glim),type="n",xlab="u",ylab="v")
#for (i in 1:n)
#{
#for (j in 1:k){
#text(bshpv[j,1,i],bshpv[j,2,i],as.character(j))
#}
#}
z$mshape<-bookmean
z$bshpv<-bshpv
z$k<-k
z$n<-n
return(z)
}
bookstein.shpv<-function(x)
{
#input x:  k x 2 matrix or complex k-vector
#output u: k x 2 matrix of Bookstein shape variables 
#           with baseline sent to (-0.5,0) (0.5,0)
	if(is.complex(x)) {
		x <- complextoreal(x)
	}
	nj <- dim(x)[1]
	j <- rep(1, times = nj)
	w <- (x[, 1] + (1i) * x[, 2] - (j * (x[1, 1] + (1i) * x[1, 2])))/(x[2, 
		1] + (1i) * x[2, 2] - x[1, 1] - (1i) * x[1, 2]) - 0.5
	w <- w[1:nj]
	y <- (Re(w))
	z <- (Im(w))
	u <- cbind(y, z)
	u <- matrix(u, nj, 2)
	u
}
bookstein.shpv.complex<-function(z)
{
#input z:  complex k vector
#output u: k-2  complex vector of Bookstein shape variables 
#           with baseline sent to (-0.5) (0.5)
	nj <- length(z)
	x <- matrix(cbind(Re(z), Im(z)), nj, 2)
	j <- rep(1, times = nj)
	w <- (x[, 1] + (1i) * x[, 2] - (j * (x[1, 1] + (1i) * x[1, 2])))/(x[2, 
		1] + (1i) * x[2, 2] - x[1, 1] - (1i) * x[1, 2]) - 0.5
	u <- w[3:nj]
	u
}
cbevec<-function(z)
{
	t1 <- reassqpr(z)
#	t2 <- eigen(t1,symmetric=TRUE,EISPACK=TRUE)
t2 <- eigen(t1,symmetric=TRUE)
	reagamma <- t2$vectors[, 1]
#	print(t2$values/sum(t2$values))

	gamma <- Vinv(reagamma)
	gamma
}
cbevectors<-function(z,j)
{
	t1 <- reassqpr(z)
#	t2 <- eigen(t1,symmetric=TRUE,EISPACK=TRUE)
	t2 <- eigen(t1,symmetric=TRUE)
	reagamma <- t2$vectors[, j]
	gamma <- Vinv(reagamma)
	gamma
}
centroid.size<-function(x)
{
#returns the centroid size of a configuration (or configurations)
#input: k x m matrix/or a complex k-vector
# or input a real k x m x n array to get a vector of sizes for a sample
        if ((is.vector(x)==FALSE)&&is.complex(x)){
k <- nrow(x)
n<-ncol(x)
tem<-array(0,c(k,2,n))
tem[,1,]<-Re(x)
tem[,2,]<-Im(x)
x<-tem
        }
{
        if(length(dim(x))==3){
        n<-dim(x)[3]
        sz<-rep(0,times=n)
	k <- dim(x)[1]
	h <- defh(k - 1)
        for (i in 1:n){
  	xh <- h %*% x[,,i]
	sz[i] <- sqrt(sum(diag(t(xh) %*% xh)))      
        }
        sz
        }
        else
        {
	if(is.vector(x)&&is.complex(x)) {
		x <- cbind(Re(x), Im(x))
	}
	k <- nrow(x)
	h <- defh(k - 1)
	xh <- h %*% x
	size <- sqrt(sum(diag(t(xh) %*% xh)))
	size
        }
}        
}
centroid.size.complex<-function(zstar)
{
#returns the centroid size of a complex vector zstar
	h <- defh(nrow(as.matrix(zstar)) - 1)
	ztem <- h %*% zstar
	size <- sqrt(diag(Re(st(ztem) %*% ztem)))
	size
}
centroid.size.mD<-function(x)
{
#returns the centroid size of a  k x m matrix
	if(is.complex(x)) {
		x <- cbind(Re(x), Im(x))
	}
	k <- nrow(x)
	h <- defh(k - 1)
	xh <- h %*% x
	size <- sqrt(sum(diag(t(xh) %*% xh)))
	size
}
complextoreal<-function(z)
{
#input complex k-vector  - return k x 2 matrix 
	nj <- length(z)
	x <- matrix(cbind(Re(z), Im(z)), nj, 2)
	x
}
defh<-function(nrow)
{
#Defines and returns an nrow x (nrow+1) Helmert sub-matrix
	k <- nrow
	h <- matrix(0, k, k + 1)
	j <- 1
	while(j <= k) {
		jj <- 1
		while(jj <= j) {
			h[j, jj] <- -1/sqrt(j * (j + 1))
			jj <- jj + 1
		}
		h[j, j + 1] <- j/sqrt(j * (j + 1))
		j <- j + 1
	}
	h
}
full.procdist<-function(x, y)
{
#input  k x 2 matrices x, y
#output full Procrustes distance rho between x,y
	sin(riemdist(x, y))
}
genpower<-function(Be, alpha)
{
	k <- dim(Be)[1]
	if(alpha == 0) {
		gen <- diag(rep(1, times = k))
	}
	else {
		l <- k - 3
#		eb <- eigen(Be, symmetric = TRUE,EISPACK=TRUE)
		eb <- eigen(Be, symmetric = TRUE)
		ev <- c(eb$values[1:l]^( - alpha/2), 0, 0, 0)
		gen <- eb$vectors %*% diag(ev) %*% t(eb$vectors)
		gen
	}
}
isotropy.test<-function(sd, p, n)
{
#LR test for isotropy with Bartlett adjustment 
#in: sd - square roots of eigenvalues of covariance matrix
#     p - the number of larger eigenvalues to consider
#     n - sample size
#out: z$bartlett  - test statistic (e.g. see Mardia, Kent, Bibby, 1979, p235)
#     z$pval - p-value
	z <- list(bartlett = 0, pval = 0)
	tem <- sd^2
	bartlett <- (log(mean(tem[1:p])) - mean(log(tem[1:p]))) * p * (n - (2 * 
		p + 11)/6)
	pval <- 1 - pchisq(bartlett, ((p + 2) * (p - 1))/2)
	z$bartlett <- bartlett
	z$pval <- pval
	return(z)
}
linegrid<-function(ref, kx, ky)
{
	n <- ky
	m <- kx
	w <- n * m
	newgrid1 <- matrix(0, w, 2)
	v <- m * 0.5
	k <- 0
	for(l in 1:v) {
		k <- k + 1
		a <- (n + m - 1) * (k - 1) + 1
		b <- n * ((2 * k) - 1)
		d <- 2 * n * k
		for(j in a:b) {
			newgrid1[j,  ] <- ref[j,  ]
		}
		for(u in 1:n) {
			down <- d - u + 1
			up <- b + u
			newgrid1[up,  ] <- ref[down,  ]
		}
	}
	newgrid2 <- matrix(0, w, 2)
	for(i in 1:v) {
		z <- (2 * i) - 1
		for(x in 1:m) {
			r1 <- m * (z - 1) + x
			e <- n * (x - 1) + z
			newgrid2[r1,  ] <- ref[e,  ]
		}
	}
	y <- v - 1
	for(p in 1:y) {
		f <- 2 * p
		for(q in 1:m) {
			r2 <- m * (f - 1) + q
			s <- n * (m - 1) + f - n * (q - 1)
			newgrid2[r2,  ] <- ref[s,  ]
		}
	}
	order <- rbind(newgrid1, newgrid2)
	order
}
mahpreshapedist<-function(z, m, pcar, pcasdev)
{
if (is.double(z)==TRUE) z<-realtocomplex(z)
if (is.double(m)==TRUE) m<-realtocomplex(m)
	w <- preshape(z)
	y <- preshape(m)
	zp <- project(w, y)
	k <- length(pcasdev)/2
	if(pcasdev[2 * k - 1] < 1e-07)
		pcasdev[2 * k - 1] <- 1e+22
	if(pcasdev[2 * k] < 1e-07)
		pcasdev[2 * k] <- 1e+22
	Sinv <- (pcar) %*% diag(1/pcasdev^2) %*% t(pcar)
	Z <- V(zp)
	d2 <- t(Z) %*% Sinv %*% (Z)
	dist <- sqrt(d2)
	dist
}
makearray<-function(x,k,m,n){
#makes a k x m x n array from a dataset read in as a table
tem<-c(t(x))
tem <- array(tem, c(m, k, n))
tem <- aperm(tem, c(2, 1, 3))
tem
}
 
movie<-function(mean, pc, sd, xl, xu, yl, yu, lineorder,movielength=20)
{
	k <- length(mean)/2
	for(i in 1:movielength) {
		plotPDMnoaxis(mean, pc * (-1)^i, sd, xl, xu, yl, yu, lineorder)
	}
	plot(mean[c(1:k)], mean[c((k + 1):(2 * k))], xlim = c(xl, xu), ylim = c(
		yl, yu), xlab = " ", ylab = " ", axes = FALSE)
}
Enorm<-function(X)
{
#finds Euclidean norm of real matrix X
	if(is.complex(X)) {
		n <- sqrt(Re(c(st(X) %*% X)))
	}
	else {
		n <- sqrt(sum(diag(t(X) %*% X)))
	}
	n
}
partial.procdist<-function(x, y)
{
#input  k x 2 matrices x, y
#output partial Procrustes distance rho between x,y
	sqrt(2) * sqrt(1 - cos(riemdist(x, y)))
}
partialwarpgrids<-function(TT, YY, xbegin, ybegin, xwidth, nr, nc, mag)
{
#
#affine grid and partial warp grids for the TPS deformation of TT to YY
#displayed as an nr x nc array of plots
#mag = magnification effect
	k <- nrow(TT)
	YY <- TT + (YY - TT) * mag
	xstart <- xbegin
	ystart <- ybegin
	kx <- 22
	ky <- 21
	l <- kx * ky
	step <- xwidth/(kx - 1)
	r <- 0
	X <- rep(0, times = 220)
	Y2 <- rep(0, times = 220)
	for(p in 1:kx) {
		ystart <- ybegin
		xstart <- xstart + step
		for(q in 1:ky) {
			ystart <- ystart + step
			r <- r + 1
			X[r] <- xstart
			Y2[r] <- ystart
		}
	}
	refc <- matrix(c(X, Y2), kx * ky, 2)
	TPS <- bendingenergy(TT)
	gamma11 <- TPS$gamma11
	gamma21 <- TPS$gamma21
	gamma31 <- TPS$gamma31
	W <- gamma11 %*% YY
	ta <- t(gamma21 %*% YY)
	B <- gamma31 %*% YY
	WtY <- t(W) %*% YY
	R <- matrix(0, k, 2)
	par(mfrow = c(nr, nc))
	par(pty = "s")  #AFFINEPART
	phi <- matrix(0, l, 2)
	s <- matrix(0, k, 1)
	for(i in 1:l) {
		s <- matrix(0, k, 1)
		for(m in 1:k) {
			s[m,  ] <- sigma(refc[i,  ] - TT[m,  ])
		}
		phi[i,  ] <- ta + t(B) %*% refc[i,  ]
	}
	newpt <- matrix(0, k, 2)
	for(i in 1:k) {
		s <- matrix(0, k, 1)
		for(m in 1:k) {
			s[m,  ] <- sigma(TT[i,  ] - TT[m,  ])
		}
		newpt[i,  ] <- ta + t(B) %*% TT[i,  ]
	}
	order <- linegrid(phi, kx, ky)
	plot(order[1:l, 1], order[1:l, 2], type = "l", xlim = c(xbegin - xwidth/
		10, xbegin + (xwidth * 11)/10), ylim = c(ybegin - (xwidth/10 * 
		ky)/kx, ybegin + ((xwidth * 11)/10 * ky)/kx), xlab = " ", ylab
		 = " ")
	lines(order[(l + 1):(2 * l), 1], order[(l + 1):(2 * l), 2], type = "l")
	points(newpt, cex = 2)
	for(jnw in 1:(k - 3)) {
		nw <- k - 2 - jnw
		phi <- matrix(0, l, 2)
		s <- matrix(0, k, 1)
		for(i in 1:l) {
			s <- matrix(0, k, 1)
			for(m in 1:k) {
				s[m,  ] <- sigma(refc[i,  ] - TT[m,  ])
			}
			phi[i,  ] <- refc[i,  ] + TPS$prinwarpeval[nw] * t(YY) %*% 
				TPS$prinwarps[, nw] %*% t(TPS$prinwarps[, nw]) %*% 
				s
		}
		newpt <- matrix(0, k, 2)
		for(i in 1:k) {
			s <- matrix(0, k, 1)
			for(m in 1:k) {
				s[m,  ] <- sigma(TT[i,  ] - TT[m,  ])
			}
			newpt[i,  ] <- TT[i,  ] + TPS$prinwarpeval[nw] * t(YY) %*% 
				TPS$prinwarps[, nw] %*% t(TPS$prinwarps[, nw]) %*% 
				s
		}
		R <- newpt - TT + R
		order <- linegrid(phi, kx, ky)
		plot(order[1:l, 1], order[1:l, 2], type = "l", xlim = c(xbegin - 
			xwidth/10, xbegin + (xwidth * 11)/10), ylim = c(ybegin - (
			xwidth/10 * ky)/kx, ybegin + ((xwidth * 11)/10 * ky)/kx
			), xlab = " ", ylab = " ")
		lines(order[(l + 1):(2 * l), 1], order[(l + 1):(2 * l), 2], 
			type = "l")
		points(newpt, cex = 2)
	}
#percentage  (need to normalize)
	d2 <- sin(riemdist(YY, TT))^2
	d3 <- sin(riemdist(R + TT, TT))^2
	percentaff <- (d2 - d3)/d2 * 100
	print("percent affine")
	print(percentaff)
}
partialwarps<-function(mshape, rotated)
{
#obtain the affine and partial warp scores for a dataset 
#where the reference configuration is mshape and the full procrustes
#rotated figures are given in the array rotated
#output: y$pwpwercent percentage of variability (squared Procrustes distance)
#  in the direction of each of the affine and principal warps
#        y$pwscores: the affine and partial warps scores
#
	y <- list(pwpercent = 0, pwscores = 0, unpercent = 0)
	k <- nrow(mshape)
	n <- dim(rotated)[3]
	msh <- mshape
	rot <- rotated
	TPS <- bendingenergy(msh)
	FX <- rot[, 1,  ]
	FY <- rot[, 2,  ]
	U <- TPS$prinwarps[, 1:(k - 3)]
	partialX <- t(U) %*% FX
	partialY <- t(U) %*% FY
	Un <- TPS$Un
	UnXY <- t(Un) %*% rbind(FX, FY)
	scores <- matrix(0, 2 * (k - 3), n)
	for(i in 1:(k - 3)) {
		r <- 2 * i - 1
		scores[r,  ] <- partialX[k - 2 - i,  ]
		scores[r + 1,  ] <- partialY[k - 2 - i,  ]
	}
	scores <- rbind(UnXY, scores)
	percwarp <- rep(0, times = (k - 2))
	sumev <- sum(eigen(var(t(scores)))$values)
	for(i in 1:(k - 2)) {
		sum1 <- sum(eigen(var(t(scores[(2 * i - 1):(2 * i),  ])))$
			values)
		percwarp[i] <- sum1/sumev
	}
	unpercent <- c(0, 0)
	unpercent[1] <- var(scores[1,  ])/sumev
	unpercent[2] <- var(scores[2,  ])/sumev
	y$unpercent <- unpercent
	y$pwpercent <- percwarp
	y$pwscores <- t(scores)
	return(y)
}
plot2rwscores<-function(rwscores, rw1, rw2, ng1, ng2)
{
	par(pch = "x")
	glim <- max( - min(rwscores), max(rwscores))
	plot(rwscores[1:ng1, rw1], rwscores[1:ng1, rw2], xlim = c( - glim, glim
		), ylim = c( - glim, glim), xlab = " ", ylab = " ")
	par(pch = "+")
	points(rwscores[(ng1 + 1):(ng1 + ng2), rw1], rwscores[(ng1 + 1):(ng1 + 
		ng2), rw2])
}
plotPDM<-function(mean, pc, sd, xl, xu, yl, yu, lineorder)
{
	for(i in c(-3,0,3)) {
		fig <- mean + i * pc * sd
		k <- length(mean)/2
		figx <- fig[1:k]
		figy <- fig[(k + 1):(2 * k)]
		plot(figx, figy, axes = TRUE, xlab = " ", ylab = " ", ylim = c(yl, 
			yu), xlim = c(xl, xu))  #               par(lty = i + 1)
		lines(figx[lineorder], figy[lineorder])
            if (i == -3) title(sub="mean - c sd")
            if (i == 0) title(sub="mean")
            if (i == 3) title(sub="mean + c sd")
		par(lty = 1)
	}
}
plotPDM2<-function(mean, pc, sd, xl, xu, yl, yu, lineorder)
{
	par(lty = 1)
	k <- length(mean)/2
	plot(mean[1:k], mean[(k + 1):(2 * k)], axes = TRUE, xlab = " ", ylab = 
		"  ", ylim = c(yl, yu), xlim = c(xl, xu))
	for(i in c(-3:3)) {
		fig <- mean + i * pc * sd
		figx <- fig[1:k]
		figy <- fig[(k + 1):(2 * k)]    #       
		if(i < 0) {
			par(lty = 1)
			par(pch = "*")
		}
		if(i == 0) {
			par(lty = 4)
			par(pch = 1)
		}
		if(i > 0) {
			par(lty = 2)
			par(pch = "+")
		}
		points(figx, figy)
		lines(figx[lineorder], figy[lineorder])
	}
}
plotPDM3<-function(mean, pc, sd, xl, xu, yl, yu, lineorder)
{
	par(lty = 1)
	k <- length(mean)/2
	figx <- matrix(0, 2 * k, 7)
	figy <- figx
	plot(mean[1:k], mean[(k + 1):(2 * k)], axes = TRUE, xlab = " ", ylab = 
		"  ", ylim = c(yl, yu), xlim = c(xl, xu))
	for(i in c(-3:3)) {
		fig <- mean + i * pc * sd
		figx[, i + 4] <- fig[1:k]
		figy[, i + 4] <- fig[(k + 1):(2 * k)]
	}
	for(i in 1:k) {
#               par(lty = 2)
#               lines(figx[i, 1:4], figy[i, 1:4])
		par(lty = 1)
		lines(figx[i, 4:7], figy[i, 4:7])
	}
}
plotPDMbook<-function(mean, pc, sd, xl, xu, yl, yu, lineorder)
{
	par(lty = 1)
	k <- length(mean)/2
	figx <- matrix(0, 2 * k, 7)
	figy <- figx
	plot(bookstein.shpv(cbind(mean[1:k], mean[(k + 1):(2 * k)])), axes = TRUE, 
		xlab = " ", ylab = "  ", ylim = c(yl, yu), xlim = c(xl, xu))
	for(i in c(-3:3)) {
		fig <- mean + i * pc * sd
		figx[, i + 4] <- fig[1:k]
		figy[, i + 4] <- fig[(k + 1):(2 * k)]
		u <- bookstein.shpv(cbind(figx[, i + 4], figy[, i + 4]))
		figx[, i + 4] <- u[, 1]
		figy[, i + 4] <- u[, 2]
	}
	for(i in 1:k) {
#               par(lty = 2)
#               lines(figx[i, 1:4], figy[i, 1:4])
		par(lty = 1)
		lines(figx[i, 4:7], figy[i, 4:7])
	}
}
plotPDMnoaxis<-function(mean, pc, sd, xl, xu, yl, yu, lineorder)
{
	for(i in c(-3:3)) {
		fig <- mean + i * pc * sd
		k <- length(mean)/2
		figx <- fig[1:k]
		figy <- fig[(k + 1):(2 * k)]
		plot(figx, figy, axes = FALSE, xlab = " ", ylab = " ", ylim = c(yl, 
			yu), xlim = c(xl, xu))
		lines(figx[lineorder], figy[lineorder])
            for (ii in 1:1000) { aa<-1}
	}
}
pointsPDMnoaxis3<-function(mean, pc, sd, xl, xu, yl, yu, lineorder,i)
{
		fig <- mean + i * pc * sd
		k <- length(mean)/2
		figx <- fig[1:k]
		figy <- fig[(k + 1):(2 * k)]
 		points(figx, figy)
                text(figx, figy,1:k)
		lines(figx[lineorder], figy[lineorder])
}
plotpairscores<-function(scores, nr, nc, ng1, ng2, ch1, ch2)
{
#plots pairs of scores score 2 vs score 1, score 4 vs score 3 etc
#in an nr x nc grid of plots
	par(pty = "s")
	par(cex = 2)
	par(mfrow = c(nr, nc))
	k <- ncol(scores)/2 + 2
	glim <- max( - min(scores), max(scores))
	for(i in 1:(k - 2)) {
		plot(scores[1:ng1, (2 * i - 1)], scores[1:ng1, (2 * i)], pch =
			ch1, xlim = c( - glim, glim), ylim = c( - glim, glim),
			xlab = " ", ylab = " ")
		points(scores[(ng1 + 1):(ng1 + ng2), (2 * i - 1)], scores[(ng1 +
			1):(ng1 + ng2), (2 * i)], pch = ch2)
	}
}
#################################
plotpca<-function (proc, pcno, type, mag, xl, yl, width, joinline = c(1,
    1),project=c(1,2))
{
    k <- proc$k
    zero <- matrix(0, k - 1, k)
    h <- defh(k - 1)
    H <- cbind(h, zero)
    H1 <- cbind(zero, h)
    H <- rbind(H, H1)
    if (project[1]==1){
    select1<-1:k
    }
    if (project[1]==2){
    select1<-(k+1):(2*k)
    }
    if (project[1]==3){
    select1<-(2*k+1):(3*k)
    }
   if (project[2]==1){
    select2<-1:k
    }
    if (project[2]==2){
    select2<-(k+1):(2*k)
    }
    if (project[2]==3){
    select2<-(2*k+1):(3*k)
    }
    select<-c(select1,select2)
    meanxy <- c(proc$mshape[, project[1]], proc$mshape[,project[2]])

    if (dim(proc$pcar)[1] == (2 * (k - 1))) {
        pcarot <- (t(H) %*% proc$pcar)[select,]
    }
    if (dim(proc$pcar)[1] != (2 * (k-1))) {
        pcarot <- proc$pcar[select,]
    }
    par(pty = "s")
    par(lty = 1)
    np <- length(pcno)
    nr <- trunc((length(pcno) + 1)/2)
    if (type == "g") {
        par(mfrow = c(nr, 2))
        if (np == 1) {
            par(mfrow = c(1, 1))
        }
        for (i in 1:np) {
            j <- pcno[i]
            fig <- meanxy + pcarot[, j] * 3 * mag * proc$pcasd[j]
            figx <- fig[1:k]
            figy <- fig[(k + 1):(2 * k)]
            YY <- cbind(figx, figy)
            tpsgrid(cbind(proc$mshape[,project[1]],proc$mshape[,project[2]])
	               , YY, xl, yl, width, 1, 0.1, 22)
        }
    }
    else {
        if (type == "r") {
            par(mfrow = c(np, 3))
            for (i in 1:np) {
                j <- pcno[i]
                plotPDM(meanxy, pcarot[, j], mag * proc$pcasd[j],
                  xl, xl + width, yl, yl + width, joinline)
                title(as.character(paste("PC ", as.character(pcno[i]),
                  ": ", as.character(round(proc$percent[i], 1)),
                  "%")))
            }
        }
        else {
            if (type == "v") {
                par(mfrow = c(nr, 2))
                if (np == 1) {
                  par(mfrow = c(1, 1))
                }
                for (i in 1:np) {
                  j <- pcno[i]
                  plotPDM3(meanxy, pcarot[, j], mag * proc$pcasd[j],
                    xl, xl + width, yl, yl + width, joinline)
                  title(as.character(paste("PC ", as.character(pcno[i]),
                    ": ", as.character(round(proc$percent[i],
                      1)), "%")))
                }
            }
            else {
                if (type == "b") {
                  par(mfrow = c(nr, 2))
                  if (np == 1) {
                    par(mfrow = c(1, 1))
                  }
                  for (i in 1:np) {
                    j <- pcno[i]
                    plotPDMbook(meanxy, pcarot[, j], mag * proc$pcasd[j],
                      -0.6, 0.6, -0.6, 0.6, joinline)
                    title(as.character(paste("PC ", as.character(pcno[i]),
                      ": ", as.character(round(proc$percent[i],
                        1)), "%")))
                  }
                }
                else {
                  if (type == "s") {
                    par(mfrow = c(nr, 2))
                    if (np == 1) {
                      par(mfrow = c(1, 1))
                    }
                    for (i in 1:np) {
                      j <- pcno[i]
                      plotPDM2(meanxy, pcarot[, j], mag * proc$pcasd[j],
                        xl, xl + width, yl, yl + width, joinline)
                      title(as.character(paste("PC ", as.character(pcno[i]),
                        ": ", as.character(round(proc$percent[i],
                          1)), "%")))
                    }
                  }
                  else {
                    if (type == "m") {
                      par(mfrow = c(1, 1))
                      for (i in 1:np) {
                        j <- pcno[i]
                        cat(paste("PC ", pcno[i], " \n"))
                        movie(meanxy, pcarot[, j], mag * proc$pcasd[j],
                          xl, xl + width, yl, yl + width, joinline,
                          20)
                      }
                    }
                  }
                }
            }
        }
    }
    par(mfrow = c(1, 1))
}

##############################################





plotprinwarp<-function(TT, xbegin, ybegin, xwidth, nr, nc)
{
#
#plots the principal warps of TT as perspective plots
#the plots are displayed  in an nr x nc array of plots
	kx <- 21
	k <- nrow(TT)
	l <- kx^2
	xstart0 <- xbegin
	ystart0 <- ybegin
	xstart <- xstart0
	ystart <- ystart0
	step <- xwidth/kx
	r <- 0
	X <- rep(0, times = l)
	Y2 <- rep(0, times = l)
	for(p in 1:kx) {
		ystart <- ystart0
		xstart <- xstart + step
		for(q in 1:kx) {
			ystart <- ystart + step
			r <- r + 1
			X[r] <- xstart
			Y2[r] <- ystart
		}
	}
	refperp <- matrix(c(X, Y2), l, 2)
	xstart <- xstart0
	xgrid <- rep(0, times = kx)
	for(i in 1:kx) {
		xstart <- xstart + step
		xgrid[i] <- xstart
	}
	ystart <- ystart0
	ygrid <- rep(0, times = kx)
	for(i in 1:kx) {
		ystart <- ystart + step
		ygrid[i] <- ystart
	}
	TPS <- bendingenergy(TT)
	prinwarp <- TPS$prinwarps
	phi <- matrix(0, l, k - 3)
	s <- matrix(0, k, 1)
	for(i in 1:l) {
		s <- matrix(0, k, 1)
		for(m in 1:k) {
			s[m,  ] <- sigma(refperp[i,  ] - TT[m,  ])
		}
		phi[i,  ] <- diag(sqrt(TPS$prinwarpeval[1:(k - 3)])) %*% t(
			prinwarp[, 1:(k - 3)]) %*% s
	}
	phiTT <- matrix(0, k, k - 3)
	for(i in 1:k) {
		s <- matrix(0, k, 1)
		for(m in 1:k) {
			s[m,  ] <- sigma(TT[i,  ] - TT[m,  ])
		}
		phiTT[i,  ] <- diag(sqrt(TPS$prinwarpeval[1:(k - 3)])) %*% t(
			prinwarp[, 1:(k - 3)]) %*% s
	}
	par(mfrow = c(nr, nc))
	for(nw in 1:(k - 3)) {
		zgrid <- matrix(0, kx, kx)
		m <- 0
		for(i in 1:kx) {
			for(j in 1:kx) {
				m <- m + 1
				zgrid[i, j] <- phi[m, k - 2 - nw]
			}
		}
		zpersp <- persp(xgrid, ygrid, zgrid, axes = TRUE)
#  NB the following is an S-Plus function : use trans3d() in R
#		points(perspp(TT[, 1], TT[, 2], phiTT[, k - 2 - nw], zpersp), 
#			cex = 2)
	}
}
plotproc<-function(proc, xl, yl, width, joinline=c(1,1))
{
#provides plots of the full Procrustes rotated objects in proc
#proc is an S object of the type output from the function procrustes2d
#xl, yl lower xlimit and ylimit in plot
#width = width (and height) of the square plotting region
	par(pty = "s")
	plot(proc$rotated[,  , 1], xlim = c(xl, xl + width), ylim = c(yl, yl + 
		width),type="n", xlab = "", ylab = "")
	for(i in 1:proc$n) {
		points(proc$rotated[,  , i])
		lines(proc$rotated[joinline,  , i])
	}
}
plotrelwarp<-function(mshape, rotsd, pcno, type, mag, xl, yl, width, joinline)
{
#provides PC plots: similar to plotpca but different argument
#here rotsd is the rotation x s.d. , and can be from the usual 
# PCA or from using relative warps 
#pcno is a vector of the numbers (index) of PCs to be plotted
#e.g. pcno<-c(1,2,4,7) will plot the four PCs no. 1,2,4,7
#type = type of display
#  "r" : rows along PCs evaluated at c = -3,-2,-1,0,1,2,3 sd's along PC
#  "v" : vectors drawn from mean to +/- 3 sd's along PC
#  "b" : vectors drawn as in `v' but using Bookstein shape variables
#  "s" : plots along c= -3, -2, -1, 0, 1, 2, 3 superimposed
#  "m" : movie backward and forwards from -3 to +3 sd's along PC
#
#mag = magnification of effect (1 = use s.d.'s from the data)
#xl, yl lower xlimit and ylimit in plot
#width = width (and height) of the square plotting region
#joinline = vector of landmark numbers which are joined up in the plot by 
#straight lines: joinline = c(1,1) will give no lines
# 
	k <- nrow(mshape)
	pcarot <- rotsd
	par(pty = "s")
	par(lty = 1)
	meanxy <- c(mshape[, 1], mshape[, 2])
	np <- length(pcno)
	if(type == "g") {
		par(mfrow = c(1, np))
		for(i in 1:np) {
			j <- pcno[i]
			fig <- meanxy + pcarot[, j] * mag * 3
			figx <- fig[1:k]
			figy <- fig[(k + 1):(2 * k)]
			YY <- cbind(figx, figy)
			tpsgrid(mshape, YY, xl, yl, width, 1, 0.1, 22)
		}
	}
	else {
		if(type == "r") {
			par(mfrow = c(np, 7))
			for(i in 1:np) {
				j <- pcno[i]
				plotPDM(meanxy, pcarot[, j], mag, xl, xl + 
				  width, yl, yl + width, joinline)
			}
		}
		else {
			if(type == "v") {
				par(mfrow = c(1, np))
				for(i in 1:np) {
				  j <- pcno[i]
				  plotPDM3(meanxy, pcarot[, j], mag, xl, xl + 
				    width, yl, yl + width, joinline)
				}
			}
			else {
				if(type == "b") {
				  par(mfrow = c(1, np))
				  for(i in 1:np) {
				    j <- pcno[i]
				    plotPDMbook(meanxy, pcarot[, j], mag, xl, 
				      xl + width, yl, yl + width, joinline)
				  }
				}
				else {
				  if(type == "s") {
				    par(mfrow = c(1, np))
				    for(i in 1:np) {
				      j <- pcno[i]
				      plotPDM2(meanxy, pcarot[, j], mag, xl, xl +
					width, yl, yl + width, joinline)
				    }
				  }
				  else {
				    if(type == "m") {
				      par(mfrow = c(1, 1))
				      for(i in 1:np) {
					j <- pcno[i]
					movie(meanxy, pcarot[, j], mag, xl, xl + 
					  width, yl, yl + width, joinline,20)
				      }
				    }
				  }
				}
			}
		}
	}
	par(mfrow = c(1, 1))
}
preshape<-function(x)
{
#input k x m matrix / complex k-vector
#output k-1 x m matrix / k-1 x 1 complex matrix
	if(is.complex(x)) {
		k <- nrow(as.matrix(x))
		h <- defh(k - 1)
		zstar <- x
		ztem <- h %*% zstar
		size <- sqrt(diag(Re(st(ztem) %*% ztem)))
		if(is.vector(zstar))
			z <- ztem/size
		if(is.matrix(zstar))
			z <- ztem %*% diag(1/size)
	}
	else {
		if(length(dim(x)) == 3) {
			k <- dim(x)[1]
			h <- defh(k - 1)
			n <- dim(x)[3]
			m <- dim(x)[2]
			z <- array(0, c(k - 1, m, n))
			for(i in 1:n) {
				z[,  , i] <- h %*% x[,  , i]
				size <- centroid.size(x[,  , i])
				z[,  , i] <- z[,  , i]/size
			}
		}
		else {
			k <- nrow(as.matrix(x))
			h <- defh(k - 1)
			ztem <- h %*% x
			size <- centroid.size(x)
			z <- ztem/size
		}
	}
	z
}
preshape.mD<-function(x)
{
#input k x m matrix
#output k-1 x 1 matrix 
	h <- defh(nrow(x) - 1)
	ztem <- h %*% x
	size <- centroid.size.mD(x)
	z <- ztem/size
	z
}
preshape.mat<-function(zstar)
{
	h <- defh(nrow(as.matrix(zstar)) - 1)
	ztem <- h %*% zstar
	size <- sqrt(diag(Re(st(ztem) %*% ztem)))
	if(is.vector(zstar))
		z <- ztem/size
	if(is.matrix(zstar))
		z <- ztem %*% diag(1/size)
	z
}
preshapetoicon<-function(z)
{
#convert a preshape (real or complex) to an icon in configuration space
	h <- defh(nrow(z))
	t(h) %*% z
}
#
#
#
#prcomp1<-function(x, retx = TRUE)
#{
#	s <- svd(scale(x, scale = FALSE), nu = 0)	# remove column means
#	rank <- sum(s$d > 0)
#	if(rank < ncol(x))
#		s$v <- s$v[, 1:rank]
#	s$d <- s$d/sqrt(max(1, nrow(x) - 1))
#	if(retx)
#		list(sdev = s$d, rotation = s$v, x = x %*% s$v)
#	else list(sdev = s$d, rotation = s$v)
#}
prinwscoregrids<-function(TT, TPS, score, xbegin, ybegin, xwidth, nr, nc)
{
#grids displaying the effect of each principal warp at `score' 
#along each warp. Grids displayed in an nr x nc array
	par(pty = "s")
	par(mfrow = c(nr, nc))
	k <- nrow(TT)
	xstart <- xbegin
	ystart <- ybegin
	kx <- 22
	ky <- 21
	l <- kx * ky
	step <- xwidth/(kx - 1)
	r <- 0
	X <- rep(0, times = 220)
	Y2 <- rep(0, times = 220)
	for(p in 1:kx) {
		ystart <- ybegin
		xstart <- xstart + step
		for(q in 1:ky) {
			ystart <- ystart + step
			r <- r + 1
			X[r] <- xstart
			Y2[r] <- ystart
		}
	}
	refc <- matrix(c(X, Y2), kx * ky, 2)    #       TPS <- bendingenergy(TT)
	for(jnw in 1:(k - 3)) {
		nw <- k - 2 - jnw
		phi <- matrix(0, l, 2)
		s <- matrix(0, k, 1)
		for(i in 1:l) {
			s <- matrix(0, k, 1)
			for(m in 1:k) {
				s[m,  ] <- sigma(refc[i,  ] - TT[m,  ])
			}
			phi[i,  ] <- refc[i,  ] + sqrt(TPS$prinwarpeval[nw]) * 
				score * t(TPS$prinwarps[, nw]) %*% s
		}
		newpt <- matrix(0, k, 2)
		for(i in 1:k) {
			s <- matrix(0, k, 1)
			for(m in 1:k) {
				s[m,  ] <- sigma(TT[i,  ] - TT[m,  ])
			}
			newpt[i,  ] <- TT[i,  ] + sqrt(TPS$prinwarpeval[nw]) * 
				score * t(TPS$prinwarps[, nw]) %*% s
		}
		order <- linegrid(phi, kx, ky)
		plot(order[1:l, 1], order[1:l, 2], type = "l", xlim = c(xbegin - 
			xwidth/10, xbegin + (xwidth * 11)/10), ylim = c(ybegin - (
			xwidth/10 * ky)/kx, ybegin + ((xwidth * 11)/10 * ky)/kx
			), xlab = " ", ylab = " ")
		lines(order[(l + 1):(2 * l), 1], order[(l + 1):(2 * l), 2], 
			type = "l")
		points(newpt, cex = 2)
	}
}
procdistreflect<-function(x, y)
{
#input  k x m matrices x, y
#output reflection shape distance (rho*) between them
#if x, y are not too far apart then (rho*)=rho (Riemannian dist)
   if (sum((x - y)^2) == 0) {
        riem <- 0
    }
    if (sum((x - y)^2) != 0) {
	m <- ncol(x)
	z <- preshape(x)
	w <- preshape(y)
	Q <- t(z) %*% w %*% t(w) %*% z
	ev <- sqrt(eigen(Q, symmetric = TRUE)$values)
#	riem <- acos(sum(ev))
	riem<-acos(min(sum(ev),1))
	} 
	riem
}
procrustes2d<-function(x, l1=1, l2=2, approxtangent=FALSE, expomap=FALSE)
{
#input k x 2 x n real array, or k x n complex matrix
#mean shape will have landmarks l1, l2 horizontal (l1 left, l2 right)
#
#output:
# z$k : no of landmarks
# z$m : no of dimensions (=2 here)
# z$n : sample size
# z$tan : the real 2k-2 x n matrix of partial Procrustes tangent coordinates
#          with pole given by the preshape of the full Procrustes mean
# z$rotated : the k x m x n array of real full Procrustes rotated data 
# z$pcar : the columns are eigenvectors (PCs) of the sample covariance Sv of z$tan
# z$pcasd : the square roots of eigenvalues of Sv (s.d.'s of PCs)
# z$percent : the % of variability explained by the PCs
# z$scores : PC scores normalised to have unit variance
# z$rawscores : PC scores (unnormalised)
# z$size : the centroid sizes of the configurations
# z$rho : Kendall's Procrustean (Riemannian) distance rho to the mean shape
# z$rmsrho : r.m.s. of rho 
# z$rmsd1 : r.m.s. of full Procrustes distances to the mean shape d1
#
	z <- list(k = 0, m = 0, n = 0, 
	rotated = 0, 
	tan = 0, 
	pcar = 0, 
	scores= 0, 
	rawscores=0, 
	pcasd = 0,  
	percent = 0, 
		 size = 0, rho = 0, rmsrho = 0, 
		rmsd1 = 0, mshape = 0)
	if(is.complex(x) == FALSE) {
		x <- x[, 1,  ] + (1i) * x[, 2,  ]
	}
#	cat("Procrustes 2D eigenanalysis \n")
	k <- nrow(x)
	n <- ncol(x)
	h <- defh(k - 1)
	zp <- preshape(x)
	gamma <- cbevec(zp)
	cbmean <- t(h) %*% gamma
	theta <- Arg(cbmean[l2] - cbmean[l1])
	cbmeanrot <- exp((-0-1i) * theta) * cbmean
	gamma <- h %*% cbmeanrot
	tan <- project(zp, gamma)
	icon <- array(0, c(k, 2, n))
	tanapprox<-matrix(0,2*k,n)
	size <- rep(0, times = n)
	rho <- rep(0, times = n)
	mu<-complextoreal(cbmeanrot)
	sum<-0
	for(i in 1:n) {
		tem <- tantofigurefull(tan[, i], gamma)
		icon[, 1, i] <- Re(tem)
		icon[, 2, i] <- Im(tem)
		sum<-sum+icon[,,i]
		size[i] <- centroid.size(x[, i])
		rho[i] <- riemdist(x[, i], c(cbmeanrot))
	}
	xbar<-sum/n

	rv <- Vmat(tan)
	if (approxtangent==TRUE){
	for (i in 1:n){
	tanapprox[,i]<-as.vector(icon[,,i])-as.vector(xbar)
	}
	tanapprox<-tanapprox/centroid.size(xbar)
	pca <- prcomp1(t(tanapprox))
	z$tan<-tanapprox
	} 

if (expomap==TRUE){
temp <- rv
for (i in 1:(n)){
temp[,i] <- rv[,i]/Enorm(rv[,i])*rho[i]
}
rv <- temp
}

	if (approxtangent==FALSE){
	pca <- prcomp1(t(rv))
	z$tan<-rv
	}
	z$pcar <- pca$rotation
	z$pcasd <- pca$sdev
	z$percent <- z$pcasd^2/sum(z$pcasd^2) * 100
	z$rotated <- icon
	npc <- 0
	for(i in 1:length(pca$sdev)) {
		if(pca$sdev[i] > 1e-07) {
			npc <- npc + 1
		}
	}
     z$scores <- pca$x
     z$rawscores <- pca$x
	for (i in 1:npc) {
        z$scores[, i] <- pca$x[, i]/pca$sdev[i]
    }
	z$rho <- rho
	z$size <- size
	z$mshape <- mu
	z$k <- k
	z$m <- 2
	z$n <- n
	z$rmsrho <- sqrt(mean(rho^2))
	z$rmsd1 <- sqrt(mean(sin(rho)^2))
	return(z)
}

testmeanshapes.old<-function(A,B,Hotelling=TRUE,tol1=1e05,tol2=1e05){
	if(is.complex(A)) {
		tem <- array(0, c(nrow(A), 2, ncol(A)))
		tem[, 1,  ] <- Re(A)
		tem[, 2,  ] <- Im(A)
		A <- tem
      }
	if(is.complex(B)) {
		tem <- array(0, c(nrow(B), 2, ncol(B)))
		tem[, 1,  ] <- Re(B)
		tem[, 2,  ] <- Im(B)
		B <- tem
	}
	m <- dim(A)[2]
if (Hotelling==TRUE)  {
if (m==2) {test<-Hotelling2D(A,B)}
if (m>2) {test<-Hotellingtest(A,B,tol1=tol1,tol2=tol2)}
cat("Hotelling's T^2 test: ",c("Test statistic = ",round(test$F,2)),
c("\n p-value = ",round(test$pval,4)),c("Degrees of freedom = ",
test$df1,test$df2),"\n")
}

if (Hotelling==FALSE)  {
if (m==2) {test<-Goodall2D(A,B)}
if (m>2) {test<-Goodalltest(A,B,tol1=tol1,tol2=tol2)}
cat("Goodall's F test: ",c("Test statistic = ",round(test$F,2)),
c("\n p-value = ",round(test$pval,4)),c("Degrees of freedom = ",
test$df1,test$df2),"\n")
}
test
}






procGPA<-function(x,scale=TRUE,reflect=FALSE,eigen2d=FALSE,
tol1=1e-05,tol2=tol1,tangentcoords="residual",proc.output=FALSE,distances=TRUE,pcaoutput=TRUE,
alpha = 0, affine = FALSE)
{
#
#
n<-dim(x)[length(dim(x))]
if ((n>100)&(distances==TRUE)){
print("To speed up use option distances=FALSE")
}
if ((n>100)&(pcaoutput==TRUE)){
print("To speed up use option pcaoutput=FALSE")
}



if (scale==TRUE){
if (tangentcoords=="residual"){
tangentresiduals<-TRUE
expomap<-FALSE
}
if (tangentcoords=="partial"){
tangentresiduals<-FALSE
expomap<-FALSE
}
if (tangentcoords=="expomap"){
tangentresiduals<-FALSE
expomap<-TRUE
}
}
if (scale==FALSE){
#all three options are equivalent
if (tangentcoords=="residual"){
tangentresiduals<-TRUE
expomap<-FALSE
}
if (tangentcoords=="partial"){
tangentresiduals<-TRUE
expomap<-FALSE
}
if (tangentcoords=="expomap"){
tangentresiduals<-TRUE
expomap<-FALSE
}
}

approxtangent<-tangentresiduals
	if(is.complex(x)) {
		tem <- array(0, c(nrow(x), 2, ncol(x)))
		tem[, 1,  ] <- Re(x)
		tem[, 2,  ] <- Im(x)
		x <- tem
	}
	m <- dim(x)[2]
n<-dim(x)[3]
if (reflect==FALSE){
if ((m == 2)&&(scale==TRUE)){
if (eigen2d==TRUE){
out<-procrustes2d(x,approxtangent=approxtangent,expomap=expomap)
}
else
{
out<-procrustesGPA(x,tol1,tol2,approxtangent=approxtangent,proc.output=proc.output,
distances=distances,pcaoutput=pcaoutput,reflect=reflect,expomap=expomap)
}                 
}
if ((m > 2)&&(scale==TRUE)){
out<-procrustesGPA(x,tol1,tol2,approxtangent=approxtangent,proc.output=proc.output
,distances=distances,pcaoutput=pcaoutput,reflect=reflect,expomap=expomap)
}
if (scale==FALSE){
out<-procrustesGPA.rot(x,tol1,tol2,approxtangent=approxtangent,
proc.output=proc.output,distances=distances,pcaoutput=pcaoutput,reflect=reflect,expomap=expomap)
}
}
if (reflect==TRUE){
if (scale==TRUE){
out<-procrustesGPA(x,tol1,tol2,approxtangent=approxtangent,
proc.output=proc.output,distances=distances,pcaoutput=pcaoutput,reflect=reflect,expomap=expomap)
}
if (scale==FALSE){
out<-procrustesGPA.rot(x,tol1,tol2,approxtangent=approxtangent,
proc.output=proc.output,distances=distances,pcaoutput=pcaoutput,reflect=reflect,expomap=expomap)
}
}
out$stdscores<-out$scores
out$scores<-out$rawscores

if (approxtangent==FALSE){

out$mshape <- out$mshape/centroid.size( out$mshape)
for (i in 1:n){
out$rotated[,,i]<- out$rotated[,,i]/centroid.size(out$rotated[,,i])
}

}

rw<-out
rw <- shaperw(out, alpha=alpha , affine=affine )
rw$GSS <- sum((n-1)*rw$pcasd**2)
rw
}


procrustesGPA<-function (x, tol1 = 1e-05, tol2 = 1e-05,distances=TRUE,pcaoutput=TRUE,
approxtangent=TRUE,proc.output=FALSE,reflect=FALSE,expomap=FALSE) 
{
	z <- list(k = 0, m = 0, n = 0, 
	rotated = 0, 
	tan = 0, 
	pcar = 0, 
	scores= 0, 
	rawscores=0, 
	pcasd = 0,  
	percent = 0, 
		 size = 0, rho = 0, rmsrho = 0, 
		rmsd1 = 0, mshape = 0)
    if (is.complex(x)) {
        tem <- array(0, c(nrow(x), 2, ncol(x)))
        tem[, 1, ] <- Re(x)
        tem[, 2, ] <- Im(x)
        x <- tem
    }
    k <- dim(x)[1]
    m <- dim(x)[2]
    n <- dim(x)[3]
    x<-cnt3(x)
    zgpa <- fgpa(x, tol1, tol2,proc.output=proc.output,reflect=reflect)

    if (distances == TRUE) {
    if (proc.output){cat("Shape distances and sizes calculation ...\n")}
    size <- rep(0, times = n)
    rho <- rep(0, times = n)
    size <- apply(x, 3, centroid.size)
    rho <- apply(x, 3, y <- function(x) {
        riemdist(x, zgpa$mshape)
    })
    }



    tanpartial <- matrix(0, k * m-m , n)
    ident <- diag(rep(1, times = (m * k-m )))
    gamma <- as.vector(preshape(zgpa$mshape))
    for (i in 1:n) {
        tanpartial[, i] <- (ident - gamma %*% t(gamma)) %*% 
        as.vector(preshape(zgpa$r.s.r[, , i]))
    }

if (expomap==TRUE){
temp <- tanpartial
for (i in 1:(n)){
temp[,i] <- tanpartial[,i]/Enorm(tanpartial[,i])*rho[i]
}
tanpartial<-temp
}

    tan <- zgpa$r.s.r[, 1, ] - zgpa$mshape[, 1]
    for (i in 2:m) {
        tan <- rbind(tan, zgpa$r.s.r[, i, ] - zgpa$mshape[, i])
    }

    if (pcaoutput==TRUE){
    if (proc.output){cat("PCA calculation ...\n")}
    if (approxtangent==FALSE){
    pca <- prcomp1(t(tanpartial))
    }
    if (approxtangent==TRUE){
    pca<-prcomp1(t(tan))
    }
    npc <- 0
    for (i in 1:length(pca$sdev)) {
        if (pca$sdev[i] > 1e-07) {
            npc <- npc + 1
        }
    }
    z$scores <- pca$x
    z$rawscores <- pca$x
    for (i in 1:npc) {
        z$scores[, i] <- pca$x[, i]/pca$sdev[i]
    }
    z$pcar <- pca$rotation
    z$pcasd <- pca$sdev
    z$percent <- z$pcasd^2/sum(z$pcasd^2) * 100
    }


   if (approxtangent==FALSE){
    z$tan <- tanpartial
    }
    if (approxtangent==TRUE){
    z$tan <- tan
    }


   if (distances == TRUE) {
    z$rho <- rho
    z$size <- size
    z$rmsrho <- sqrt(mean(rho^2))
    z$rmsd1 <- sqrt(mean(sin(rho)^2))
    }
 
    z$rotated <- zgpa$r.s.r    
    z$mshape <- zgpa$mshape
    z$k <- k
    z$m <- m
    z$n <- n

              if (proc.output){cat("Finished.\n")}

    return(z)
}


procrustesGPA.rot<-function (x, tol1 = 1e-05, tol2 = 1e-05, distances = TRUE, 
pcaoutput=TRUE, approxtangent=TRUE,proc.output=FALSE,reflect=FALSE, expomap=FALSE) 
{
	z <- list(k = 0, m = 0, n = 0, 
	rotated = 0, 
tan = 0, 
pcar = 0, 
scores= 0, 
rawscores=0, 
pcasd = 0,  
percent = 0, 
		 size = 0, rho = 0, rmsrho = 0, 
		rmsd1 = 0, mshape = 0)
    if (is.complex(x)) {
        tem <- array(0, c(nrow(x), 2, ncol(x)))
        tem[, 1, ] <- Re(x)
        tem[, 2, ] <- Im(x)
        x <- tem
    }
    k <- dim(x)[1]
    m <- dim(x)[2]
    n <- dim(x)[3]



#    print("GPA (rotation only)")
    x<-cnt3(x)
    zgpa <- fgpa.rot(x, tol1, tol2,proc.output=proc.output,reflect=reflect)

    if (distances == TRUE) {
       if (proc.output){cat("Shape distances and sizes calculation ...\n")}
        size <- rep(0, times = n)
        rho <- rep(0, times = n)
        size <- apply(x, 3, centroid.size)
        rho <- apply(x, 3, y <- function(x) {
            riemdist(x, zgpa$mshape)
        })
}


    
    tanpartial <- matrix(0, k * m-m, n)
    ident <- diag(rep(1, times = (m * k-m)))
    gamma <- as.vector(preshape(zgpa$mshape))
    for (i in 1:n) {
        tanpartial[, i] <- (ident - gamma %*% t(gamma)) %*% 
        as.vector(preshape(zgpa$r.s.r[, , i]))
    }

if (expomap==TRUE){
temp <- tanpartial
for (i in 1:(n)){
temp[,i] <- tanpartial[,i]/Enorm(tanpartial[,i])*rho[i]
}
tanpartial<-temp
}



    tan <- zgpa$r.s.r[, 1, ] - zgpa$mshape[, 1]
    for (i in 2:m) {
        tan <- rbind(tan, zgpa$r.s.r[, i, ] - zgpa$mshape[, i])
    }

        if (approxtangent==FALSE){
    z$tan <- tanpartial
    }
    if (approxtangent==TRUE){
    z$tan <- tan
    }


  if (pcaoutput==TRUE){
      if (proc.output){cat("PCA calculation ...\n")}

        if (approxtangent==FALSE){
    pca <- prcomp1(t(tanpartial))
    }
    if (approxtangent==TRUE){
    pca<-prcomp1(t(tan))
    }



    npc <- 0
    for (i in 1:length(pca$sdev)) {
        if (pca$sdev[i] > 1e-07) {
            npc <- npc + 1
        }
    }
    
    z$scores <- pca$x
    z$rawscores <- pca$x
    for (i in 1:npc) {
        z$scores[, i] <- pca$x[, i]/pca$sdev[i]
    }

    z$pcar <- pca$rotation
    z$pcasd <- pca$sdev
    z$percent <- z$pcasd^2/sum(z$pcasd^2) * 100
    }

    if (distances == TRUE) {
    z$rho <- rho
    z$size <- size
    z$rmsrho <- sqrt(mean(rho^2))
    z$rmsd1 <- sqrt(mean(sin(rho)^2))
    } 
    z$rotated <- zgpa$r.s.r
    z$mshape <- zgpa$mshape
    z$k <- k
    z$m <- m
    z$n <- n
           if (proc.output){cat("Finished.\n")}
    return(z)
}




project<-function(z, gamma)
{
#input z: preshape, gamma: preshape (k-1 x 1 matrices)
#output Kent's tangent plane coordinates 
#of z at the pole gamma (k-1 complex vector)
	nr <- nrow(z)
	nc <- ncol(z)
	g <- matrix(gamma, nr, 1)
	ident <- diag(nr)
	theta <- diag(c(exp((-0-1i) * Arg(st(g) %*% z))), nc, nc)
	v <- (ident - g %*% st(g)) %*% z %*% theta
	v
}
read.array<-function(name, k, m, n)
{
#input name : filename, k: no of points, m: no of dimensions, n: sample size
#output x: k x m x n array of data
#e.g. for 2D data assume file format x1 y1 x2 y2 .. xn yn for each object
	tem <- scan(name)
	tem <- array(tem, c(m, k, n))
	tem <- aperm(tem, c(2, 1, 3))
	x <- tem
	x
}
read.in<-function(name, k, m)
{
#input name : filename, k: no of points, m: no of dimensions
#output x: k x m x n array of data ( n: sample size)
#e.g. for m=2-D data assume file format x1 y1 x2 y2 ... xk yk for each object
#for m=3-D data: x1 y1 z1 x2 y2 z2 ... xk yk zk
	tem <- scan(name)
	n <- length(tem)/(k * m)
	tem <- array(tem, c(m, k, n))
	tem <- aperm(tem, c(2, 1, 3))
	x <- tem
	x
}
realtocomplex<-function(x)
{
#input k x 2 matrix - return complex k-vector 
	k <- nrow(x)
	zstar <- x[, 1] + (1i) * x[, 2]
	zstar
}
reassqpr<-function(z)
{
	j <- 1
	nc <- ncol(z)
	nr <- nrow(z)
	stemp <- matrix(0, 2 * nr, 2 * nr)
	repeat {
		t1 <- matrix(z[, j], nr, 1)
		vz <- rbind(Re(t1), Im(t1))
		viz <- rbind(Re((1i) * t1), Im((1i) * t1))
		stemp <- stemp + vz %*% t(vz) + viz %*% t(viz)
		if(j == nc)
			break
		j <- j + 1
	}
	stemp
}
relwarps<-function(mshape, rotated, alpha)
{
#find the relative warps for a dataset with mshape as the reference  
#and `rotated' as the array of Procrustes rotated figures
#alpha is the power of the bending energy
#  alpha=+1 : emphasizes large scale 
#  alpha=-1 : emphasizes small scale
#output: 
#     z$rwarps  : the relative warps
#     z$rwscores : the relative warp scores
#     z$rwpercent : the percentage of total variability explained by each            #relative warp
	z <- list(rwarps = 0, rwscores = 0, rwpercent = 0, ev = 0, unif = 0, 
		unscores = 0, lengths = 0)
	k <- nrow(mshape)
	TPS <- bendingenergy(mshape)
	Be <- TPS$gamma11
	stackxy <- rbind(rotated[, 1,  ], rotated[, 2,  ])
	n <- dim(rotated)[3]
	msum <- rep(0, times = 2 * k)
	for(i in 1:n) {
		msum <- msum + stackxy[, i]
	}
	msum <- msum/n
	meanxy <- msum
	cstackxy <- matrix(0, 2 * k, n)
	for(i in 1:n) {
		cstackxy[, i] <- stackxy[, i] - meanxy
	}
	Bpow <- genpower(Be, alpha)
	Bpowinv <- genpower(Be,  - alpha)
	IBpow <- I2mat(Bpow)
	IBpowinv <- I2mat(Bpowinv)
	if(alpha == 0) {
		IBpow <- diag(rep(1, times = (2 * k)))
		IBpowinv <- diag(rep(1, times = (2 * k)))
	}
	stacknew <- IBpow %*% cstackxy
	gamma <- matrix(0, 2 * k, 2 * k)
	pcarotation <- eigen(stacknew %*% t(stacknew)/n, symmetric = TRUE)$vectors
	pcaev <- eigen(stacknew %*% t(stacknew)/n, symmetric = TRUE)$values
	pcasdev <- rep(0, times = 2 * k)
	for(i in 1:(2 * k)) {
		pcasdev[i] <- sqrt(abs(pcaev[i]))
	}
	scores <- t(IBpow %*% pcarotation) %*% cstackxy
	percent <- rep(0, times = 2 * k)
	for(i in 1:(2 * k)) {
		percent[i] <- pcasdev[i]^2
	}
	Un <- TPS$Un
	UnXY <- t(Un) %*% cstackxy
	z$unif <- Un %*% t(matrix(c(sqrt(var(UnXY[1,  ])), 0, 0, sqrt(var(UnXY[
		2,  ]))), 2, 2))
	z$unscores <- t(UnXY)
	z$lengths <- sqrt(abs(percent))
	z$rwarps <- IBpowinv %*% pcarotation %*% diag(pcasdev)
	z$rwscores <- t(scores)
	z$ev <- pcaev
	percentrw <- percent/sum(percent) * 100
	z$rwpercent <- percentrw
	return(z)
}

ssriemdist<-function(x,y,reflect=FALSE){
sx <- centroid.size(x)
sy <- centroid.size(y)
sd <- sx**2+sy**2-2*sx*sy*cos(riemdist(x,y,reflect=reflect))
sqrt(abs(sd))
}

riemdist<-function(x, y, reflect=FALSE)
{
#input two k x m matrices x, y or complex k-vectors
#output Riemannian distance rho between them 

if (sum((x-y)**2)==0){
riem <- 0
}
if (sum((x-y)**2)!=0){

if (reflect==FALSE) { 
	if(ncol(as.matrix(x)) < 3) {
          if (is.complex(x)==FALSE){x<-realtocomplex(x)}
          if (is.complex(y)==FALSE){y<-realtocomplex(y)} 
#riem <- c(acos(Mod(st(preshape(x)) %*% preshape(y))))
riem<-c(acos(min(1,(Mod(st(preshape(x)) %*% preshape(y))))))
	}
	else {
		m <- ncol(x)
		z <- preshape(x)
		w <- preshape(y)
		Q <- t(z) %*% w %*% t(w) %*% z
		ev <- eigen(t(z) %*% w)$values
		check <- 1
		for(i in 1:m) {
			check <- check * ev[i]
		}
		ev <- sqrt(abs(eigen(Q, symmetric = TRUE)$values))
		if(Re(check) < 0)
			ev[m] <-  - ev[m]
		riem <- acos(min(sum(ev),1))
	}
}
if (reflect==TRUE){

		m <- ncol(x)
		z <- preshape(x)
		w <- preshape(y)
		Q <- t(z) %*% w %*% t(w) %*% z

		ev <- sqrt(abs(eigen(Q, symmetric = TRUE)$values))

		riem <- acos(min(sum(ev),1))
}

}
riem
}

riemdist.complex<-function(z, w)
{
#input complex k-vectors z, w
#output Riemannian distance rho between them
	c(acos(min(Mod(st(preshape(z)) %*% preshape(w)),1)))
}

riemdist.mD<-function(x, y)
{
#input  k x m matrices x, y
#output Riemannian distance rho between them
	m <- ncol(x)
	z <- preshape.mD(x)
	w <- preshape.mD(y)
	Q <- t(z) %*% w %*% t(w) %*% z
		ev <- eigen(t(z) %*% w)$values
		check <- 1
		for(i in 1:m) {
			check <- check * ev[i]
		}
	ev <- sqrt(eigen(Q, symmetric = TRUE)$values)
	if(check < 0)
		ev[m] <-  - ev[m]
	riem <- acos(min(sum(ev),1))
	riem
}
rotateaxes<-function(mshapein, rotatedin)
{
#Rotates a mean shape and the Procrustes rotated data to have 
#horizontal and vertical principal axes
#output: z$mshape rotated mean shape
#        z$rotated   rotated procrustes registered data
#        z$R  the rotation matrix
#
	z <- list(mshape = 0, rotated = 0, R = 0)
	n <- dim(rotatedin)[3]
	S <- var(mshapein)
	R <- eigen(S)$vectors
	msh <- mshapein %*% R
	ico <- rotatedin
	for(i in 1:n) {
		ico[,  , i] <- rotatedin[,  , i] %*% R
	}
	z$mshape <- msh
	z$rotated <- ico
	z$R <- R
	return(z)
}

#sigma<-function(x)
#{
#	length <- sqrt(x[1]^2 + x[2]^2)
#	if(length == 0)
#		sig <- 0
#	else sig <- length^2 * log(length^2)
#	sig
#}

sigma<-function( x )
{
#  other radial basis functions/covariance functions are possible of course

	hh <- Enorm(x)
	if( hh == 0)
		sig <- 0
	else
{
if ( length(x) == 2  ){
sig <- hh^2 * log( hh^2)   # null space includes affine terms (2D data)
}
if ( length(x)==3 ){
	sig <-  -hh   # null space includes affine terms (3D data)
}
}
sig
}


st<-function(zstar)
{
#input complex matrix
#output transpose of the complex conjugate 
	st <- t(Conj(zstar))
	st
}
tanfigure<-function(vv, gamma)
{
#inverse projection from complex tangent plane coordinates vv, using pole gamma
#output centred icon
	k <- nrow(gamma) + 1
	h <- defh(k - 1)
	zvv <- tanpreshape(vv, gamma)
	zstvv <- t(h) %*% zvv
	zstvv
}
tanfigurefull<-function(vv, gamma)
{
#inverse projection from complex tangent plane coordinates vv, using pole gamma
#using Procrustes to with scaling to the pole gamma 
#output centred icon 
	k <- nrow(gamma) + 1
	f1 <- tanfigure(vv, gamma)
	h <- defh(k - 1)
	f2 <- t(h) %*% gamma
	beta <- Mod(st(f1) %*% f2)
	f1 <- f1 * c(beta)
	f1
}
tanpreshape<-function(vv, gamma)
{
#inverse projection from tangent plane coordinates vv, using pole gamma
#output preshape
	z <- c((1 - st(vv) %*% vv)^0.5) * gamma + vv
	z
}
tantofigure<-function(vv, gamma)
{
#inverse projection from complex tangent plane coordinates vv, using pole gamma
#output centred icon
	k <- nrow(gamma) + 1
	h <- defh(k - 1)
	zvv <- tantopreshape(vv, gamma)
	zstvv <- t(h) %*% zvv
	zstvv
}
tantofigurefull<-function(vv, gamma)
{
	k <- nrow(gamma) + 1
	h <- defh(k - 1)
	f1 <- tantofigure(vv, gamma)
	f2 <- t(h) %*% gamma
	beta <- Mod(st(f1) %*% f2)
	f1 <- f1 * c(beta)
	f1
}
tantopreshape<-function(vv, gamma)
{
#inverse projection from complex tangent plane coordinates vv, using pole gamma
#output preshape
	z <- c((1 - st(vv) %*% vv)^0.5) * gamma + vv
	z
}
plot3Ddata<-function(dna.data,land=1:k,objects=1:n,joinline=c(1,1)){
dna<-procGPA(dna.data[,,1:2])
w1<-defplotsize2(dna.data[,1:2,])
w2<-defplotsize2(dna.data[,c(1,3),])
w3<-defplotsize2(dna.data[,c(2,3),])
width<-max(c(w1$width,w2$width,w3$width))
xl<-min(c(w1$xl,w2$xl,w3$xl))
xu<-xl+width
yl<-min(c(w1$yl,w2$yl,w3$yl))
yu<-yl+width
n<-dim(dna.data)[3]
k<-dim(dna.data)[1]
m<-dim(dna.data)[2]
par(mfrow=c(1,1))
par(pty="s")
view1<-1
view2<-2
view3<-3
lineorder<-joinline
for (j in 1:1){
for (ii in objects){
par(mfrow=c(2,2))
mag<-0
pcno<-1
plotPDMnoaxis3(c(dna.data[land,view2,ii],dna.data[land,view3,ii]),
c(dna$pcar[((view1-1)*k+(land)),pcno],dna$pcar[((view2-1)*k+(land)),pcno]),
mag*dna$pcasd[pcno],xl,xu,yl,yu,lineorder,1)
mag<-0
pcno<-1
plotPDMnoaxis3(c(dna.data[land,view1,ii],dna.data[land,view3,ii]),
c(dna$pcar[((view1-1)*k+(land)),pcno],dna$pcar[((view2-1)*k+(land)),pcno]),
mag*dna$pcasd[pcno],xl,xu,yl,yu,lineorder,1)
mag<-0
pcno<-1
plotPDMnoaxis3(c(dna.data[land,view1,ii],dna.data[land,view2,ii]),
c(dna$pcar[((view1-1)*k+(land)),pcno],dna$pcar[((view2-1)*k+(land)),pcno]),
mag*dna$pcasd[pcno],xl,xu,yl,yu,lineorder,1)
plot(c(0,0),c(50,50),xlim=c(0,0),ylim=c(0,0),type="n",xlab=" ",ylab=" ",axes=FALSE)
title(as.character(ii))
}
}
}

plot3Ddata.static<-function(dna.data,land=1:k,objects=1:n,joinline=c(1,1)){
dna<-procGPA(dna.data[,,1:2])
w1<-defplotsize2(dna.data[,1:2,])
w2<-defplotsize2(dna.data[,c(1,3),])
w3<-defplotsize2(dna.data[,c(2,3),])
width<-max(c(w1$width,w2$width,w3$width))
xl<-min(c(w1$xl,w2$xl,w3$xl))
xu<-xl+width
yl<-min(c(w1$yl,w2$yl,w3$yl))
yu<-yl+width
n<-dim(dna.data)[3]
k<-dim(dna.data)[1]
m<-dim(dna.data)[2]
par(mfrow=c(1,1))
par(pty="s")

lineorder<-joinline

par(mfrow=c(2,2))
mag<-0
pcno<-1
ii<-1
view1<-1
view2<-2
view3<-3
plotPDMnoaxis3(c(dna.data[land,view1,ii],dna.data[land,view2,ii]),
c(dna$pcar[((view1-1)*k+(land)),pcno],dna$pcar[((view2-1)*k+(land)),pcno]),
mag*dna$pcasd[pcno],xl,xu,yl,yu,lineorder,1)
for (ii in objects){
pointsPDMnoaxis3(c(dna.data[land,view1,ii],dna.data[land,view2,ii]),
c(dna$pcar[((view1-1)*k+(land)),pcno],dna$pcar[((view2-1)*k+(land)),pcno]),
mag*dna$pcasd[pcno],xl,xu,yl,yu,lineorder,1)
}
view1<-1
view2<-3
view3<-2
plotPDMnoaxis3(c(dna.data[land,view1,ii],dna.data[land,view2,ii]),
c(dna$pcar[((view1-1)*k+(land)),pcno],dna$pcar[((view2-1)*k+(land)),pcno]),
mag*dna$pcasd[pcno],xl,xu,yl,yu,lineorder,1)
for (ii in objects){
pointsPDMnoaxis3(c(dna.data[land,view1,ii],dna.data[land,view2,ii]),
c(dna$pcar[((view1-1)*k+(land)),pcno],dna$pcar[((view2-1)*k+(land)),pcno]),
mag*dna$pcasd[pcno],xl,xu,yl,yu,lineorder,1)
}
view1<-2
view2<-3
view3<-1
plotPDMnoaxis3(c(dna.data[land,view1,ii],dna.data[land,view2,ii]),
c(dna$pcar[((view1-1)*k+(land)),pcno],dna$pcar[((view2-1)*k+(land)),pcno]),
mag*dna$pcasd[pcno],xl,xu,yl,yu,lineorder,1)
for (ii in objects){
pointsPDMnoaxis3(c(dna.data[land,view1,ii],dna.data[land,view2,ii]),
c(dna$pcar[((view1-1)*k+(land)),pcno],dna$pcar[((view2-1)*k+(land)),pcno]),
mag*dna$pcasd[pcno],xl,xu,yl,yu,lineorder,1)
}


}

plot3Dmean<-function(dna){
land<-1:dim(dna$mshape)[1]
w1<-defplotsize2(dna$rotated[,1:2,])
w2<-defplotsize2(dna$rotated[,c(1,3),])
w3<-defplotsize2(dna$rotated[,c(2,3),])
width<-max(c(w1$width,w2$width,w3$width))
xl<-min(c(w1$xl,w2$xl,w3$xl))
xu<-xl+width
yl<-min(c(w1$yl,w2$yl,w3$yl))
yu<-yl+width
par(mfrow=c(2,2))
par(pty="s")
plot(dna$mshape[land,1],dna$mshape[land,2],xlim=c(xl,xu),ylim=c(yl,yu),xlab=" ",ylab=" ")
text(dna$mshape[land,1],dna$mshape[land,2],land)
lines(dna$mshape[land,1],dna$mshape[land,2])
plot(dna$mshape[land,1],dna$mshape[land,3],xlim=c(xl,xu),ylim=c(yl,yu),xlab=" ",ylab=" ")
text(dna$mshape[land,1],dna$mshape[land,3],land)
lines(dna$mshape[land,1],dna$mshape[land,3])
plot(dna$mshape[land,2],dna$mshape[land,3],xlim=c(xl,xu),ylim=c(yl,yu),xlab=" ",ylab=" ")
text(dna$mshape[land,2],dna$mshape[land,3],land)
lines(dna$mshape[land,2],dna$mshape[land,3])
title("Procrustes mean shape estimate")
}

plot3Dpca<-function(dna,pcno,joinline=c(1,1)){
#choose subset
w1<-defplotsize2(dna$rotated[,1:2,])
w2<-defplotsize2(dna$rotated[,c(1,3),])
w3<-defplotsize2(dna$rotated[,c(2,3),])
width<-max(c(w1$width,w2$width,w3$width))
xl<-min(c(w1$xl,w2$xl,w3$xl))-width/4
xu<-xl+width*1.5
yl<-min(c(w1$yl,w2$yl,w3$yl))-width/4
yu<-yl+width*1.5
k<-dim(dna$mshape)[1]
lineorder<-joinline
par(mfrow=c(1,1))
cat("X-Y view \n")
view1<-1
view2<-2
view3<-3
land<-c(1:k)
for (j in 1:10){
for (ii in -12:12){
mag<-ii/4
plotPDMnoaxis3(c(dna$mshape[land,view1],dna$mshape[land,view2]),
c(dna$pcar[((view1-1)*k+(land)),pcno],dna$pcar[((view2-1)*k+(land)),pcno]),
mag*dna$pcasd[pcno],xl,xu,yl,yu,lineorder,1)
}
for (ii in -11:11){
mag<--ii/4
plotPDMnoaxis3(c(dna$mshape[land,view1],dna$mshape[land,view2]),
c(dna$pcar[((view1-1)*k+(land)),pcno],dna$pcar[((view2-1)*k+(land)),pcno]),
mag*dna$pcasd[pcno],xl,xu,yl,yu,lineorder,1)
}
}
#choose subset
par(mfrow=c(1,1))
cat("X-Z view \n")
view1<-1
view2<-3
view3<-2
land<-c(1:k)
for (j in 1:10){
for (ii in -12:12){
mag<-ii/4
plotPDMnoaxis3(c(dna$mshape[land,view1],dna$mshape[land,view2]),
c(dna$pcar[((view1-1)*k+(land)),pcno],dna$pcar[((view2-1)*k+(land)),pcno]),
mag*dna$pcasd[pcno],xl,xu,yl,yu,lineorder,1)
}
for (ii in -11:11){
mag<--ii/4
plotPDMnoaxis3(c(dna$mshape[land,view1],dna$mshape[land,view2]),
c(dna$pcar[((view1-1)*k+(land)),pcno],dna$pcar[((view2-1)*k+(land)),pcno]),
mag*dna$pcasd[pcno],xl,xu,yl,yu,lineorder,1)
}
}
#choose subset
par(mfrow=c(1,1))
cat("Y-Z view \n")
view1<-2
view2<-3
view3<-1
land<-c(1:k)
for (j in 1:10){
for (ii in -12:12){
mag<-ii/4
plotPDMnoaxis3(c(dna$mshape[land,view1],dna$mshape[land,view2]),
c(dna$pcar[((view1-1)*k+(land)),pcno],dna$pcar[((view2-1)*k+(land)),pcno]),
mag*dna$pcasd[pcno],xl,xu,yl,yu,lineorder,1)
}
for (ii in -11:11){
mag<--ii/4
plotPDMnoaxis3(c(dna$mshape[land,view1],dna$mshape[land,view2]),
c(dna$pcar[((view1-1)*k+(land)),pcno],dna$pcar[((view2-1)*k+(land)),pcno]),
mag*dna$pcasd[pcno],xl,xu,yl,yu,lineorder,1)
}
}
}






banner1<-function(char)
{
par(mfrow=c(1,1))
plot(c(0,0),c(1,1),axes=FALSE,type="n",xlab=" ",ylab=" ")
a1<-char
if (length(a1)==2) a1<-paste(a1[1],a1[2])
if (length(a1)==3) a1<-paste(a1[1],a1[2],a1[3])
if (is.character(a1)==FALSE) char<-as.character(a1)
title(a1)
}

banner4<-function(a1,a2,a3,a4)
{
par(mfrow=c(2,2))
plot(c(0,0),c(1,1),axes=FALSE,type="n",xlab=" ",ylab=" ")
if (length(a1)==2) a1<-paste(a1[1],a1[2])
if (length(a1)==3) a1<-paste(a1[1],a1[2],a1[3])
if (is.character(a1)==FALSE) a1<-as.character(a1)
title(a1)
plot(c(0,0),c(1,1),axes=FALSE,type="n",xlab=" ",ylab=" ")
if (length(a2)==2) a2<-paste(a2[1],a2[2])
if (length(a2)==3) a2<-paste(a2[1],a2[2],a2[3])
if (is.character(a2)==FALSE) a2<-as.character(a2)
title(a2)
plot(c(0,0),c(1,1),axes=FALSE,type="n",xlab=" ",ylab=" ")
if (length(a3)==2) a3<-paste(a3[1],a3[2])
if (length(a3)==3) a3<-paste(a3[1],a3[2],a3[3])
if (is.character(a3)==FALSE) a3<-as.character(a3)
title(a3)
plot(c(0,0),c(1,1),axes=FALSE,type="n",xlab=" ",ylab=" ")
if (length(a4)==2) a4<-paste(a4[1],a4[2])
if (length(a4)==3) a4<-paste(a4[1],a4[2],a4[3])
if (is.character(a4)==FALSE) a4<-as.character(a4)
title(a4)
}

#data

gorf.dat<-array(c(5,193,53,-27,0,0,0,33,-2,105,18,176,72,114,92,38
,51,191,55,-31,0,0,0,33,25,106,56,171,98,105,99,15
,36,187,59,-31,0,0,0,36,12,102,38,171,91,103,100,19
,23,202,48,-30,0,0,0,39,3,103,33,180,84,112,94,28
,30,185,62,-25,0,0,0,32,11,101,37,168,85,106,96,21
,4,195,65,-21,0,0,0,34,-4,100,15,180,69,120,102,34
,37,195,62,-32,0,0,0,35,20,101,50,173,102,105,105,22
,41,191,58,-34,0,0,0,34,15,100,47,175,93,105,99,18
,40,190,52,-33,0,0,0,38,13,107,44,176,88,113,102,31
,-4,179,62,-21,0,0,0,29,1,89,9,164,70,111,100,36
,41,206,53,-25,0,0,0,39,11,104,47,177,95,111,95,26
,33,197,55,-30,0,0,0,35,7,106,39,175,89,111,95,24
,-12,205,52,-15,0,0,0,38,-10,111,4,187,66,129,80,44
,13,186,56,-32,0,0,0,34,8,101,25,166,80,105,97,26
,20,186,45,-31,0,0,0,34,10,96,31,165,84,104,90,19
,29,183,55,-31,0,0,0,32,10,98,39,163,82,106,95,17
,11,203,57,-28,0,0,0,39,-2,106,23,182,77,122,100,36
,37,187,54,-27,0,0,0,34,11,100,43,171,84,106,93,28
,49,191,53,-31,0,0,0,35,21,102,54,172,94,98,99,18
,-8,191,57,-34,0,0,0,32,-7,93,6,173,71,119,101,30
,43,184,49,-32,0,0,0,33,14,100,49,165,91,99,98,20
,57,185,62,-37,0,0,0,35,22,103,61,169,96,100,104,24
,-10,196,55,-20,0,0,0,38,-10,107,5,181,73,123,88,46
,20,195,60,-28,0,0,0,32,6,101,33,173,84,114,100,30
,35,202,59,-27,0,0,0,34,6,108,41,182,83,117,99,31
,1,188,60,-19,0,0,0,35,-2,99,12,170,70,119,93,45
,24,194,52,-24,0,0,0,39,8,105,34,174,80,115,95,32
,25,204,55,-27,0,0,0,34,7,108,35,185,83,118,92,32
,36,198,47,-30,0,0,0,39,14,110,43,177,92,105,98,25
,8,198,53,-35,0,0,0,34,4,101,22,175,82,111,100,24),c(2,8,30))
gorf.dat<-aperm(gorf.dat,c(2,1,3))

gorm.dat<-array(c(53,220,46,-35,0,0,0,37,12,122,58,204,93,117,103,28
,57,219,50,-43,0,0,0,37,13,119,61,198,102,110,104,20
,89,227,52,-47,0,0,0,32,35,120,93,201,131,92,104,4
,46,222,51,-45,0,0,0,30,11,113,54,196,101,117,101,16
,85,220,48,-38,0,0,0,39,28,125,87,203,121,106,103,7
,64,208,43,-39,0,0,0,36,22,111,67,191,104,102,101,18
,67,216,51,-37,0,0,0,35,17,119,68,191,108,109,94,15
,35,236,61,-42,0,0,0,33,2,119,43,211,90,126,104,30
,116,218,40,-38,0,0,0,41,41,124,116,201,133,94,103,12
,56,234,60,-34,0,0,0,34,12,121,58,215,109,119,112,28
,40,223,58,-36,0,0,0,34,9,113,46,202,97,120,112,24
,94,223,49,-57,0,0,0,33,31,122,94,206,136,99,113,-1
,68,222,59,-41,0,0,0,30,18,119,68,204,104,114,98,11
,65,224,56,-33,0,0,0,35,15,130,67,205,108,115,95,20
,67,214,52,-47,0,0,0,36,26,115,74,192,114,105,104,11
,110,213,52,-46,0,0,0,37,42,121,109,190,133,97,108,-8
,46,219,50,-42,0,0,0,36,11,121,56,199,104,108,102,21
,79,209,66,-43,0,0,0,35,24,114,84,193,108,115,109,14
,58,244,74,-22,0,0,0,37,7,131,64,219,98,128,100,29
,43,236,64,-43,0,0,0,33,12,124,52,215,110,121,105,7
,70,226,54,-37,0,0,0,39,28,121,74,204,122,105,107,7
,68,224,55,-37,0,0,0,35,18,121,71,207,109,108,98,13
,34,247,63,-35,0,0,0,35,4,124,45,225,104,135,110,29
,49,236,59,-40,0,0,0,38,19,127,59,219,105,121,109,25
,98,195,44,-44,0,0,0,36,30,116,98,177,121,89,105,10
,109,208,49,-40,0,0,0,36,29,125,105,189,120,102,102,12
,61,224,51,-35,0,0,0,41,15,122,67,206,107,121,103,25
,43,213,49,-57,0,0,0,28,20,111,58,194,112,106,108,6
,26,249,67,-14,0,0,0,38,-11,130,33,225,87,148,97,53),c(2,8,29))
gorm.dat<-aperm(gorm.dat,c(2,1,3))
qset2.dat<-array(c(117.98,219.62,114.52,41.93,166.15,113.59,206.54,121.79,165.11,142.92,62.07,136.58
,105.52,235.08,109.96,57.31,165.44,126.42,223.05,140.88,169.58,156.83,59.5,138.02
,142.83,222,132.08,40.2,188.8,115.39,236.9,125.08,193.18,142.22,83.37,141.47
,126.35,204.27,113.2,36.04,171.77,102.43,228.93,111.75,173.55,131.8,66.26,126.65
,99.11,231.52,119.58,44.52,169.28,129.51,219.93,141.98,167.65,154.41,56.83,141.64
,134.14,228.48,129.09,56.98,175.85,125.39,217.57,138.01,179.29,152.67,73.84,153.08
,119.8,219.48,122.36,48.52,171.44,115.58,216.57,131.13,170.86,148.36,65.8,138.8
,105.62,222.74,96.27,51.3,152.13,119.07,205.02,131.95,157.36,141.61,52.35,142.73
,122.15,202.41,128.75,32.09,176.88,102.84,220.81,119.41,177.78,131.04,76.65,120.61
,127.78,223.48,117.65,53.94,183.45,123.87,236.96,133.84,184.24,149.04,77.96,149.95
,123.4,200.6,116.57,28.43,172.67,97.33,221.34,106.24,175.29,122.05,64.69,115.53
,132.4,227.96,127.46,50.4,171.57,118.67,203.72,131.91,175.25,152.38,71.98,138.5
,123.39,216.23,113.96,29.9,167.68,105.58,202.15,114.18,172.01,133.8,65.78,125.86
,136.83,207.84,117.37,33.67,178.1,100.65,233.17,111.44,184.05,124.85,69.8,129.69
,143.31,219.71,122.1,48.08,177.23,113.52,221.38,122.32,180.86,140.87,80.11,151.72
,105.96,218.06,100.46,48.73,155.23,119.51,217.76,130.01,159.11,145.49,51.29,139.66
,115.73,234.14,115.74,56.11,168.42,128.73,220.32,135.6,169.4,154.65,66.02,152.32
,101.73,215.74,102.67,35.73,151.69,110.19,199.12,120.24,151.6,136.36,48.82,132.42
,124.93,222.42,130.09,46.89,163.86,118.81,197.76,137.59,165.11,148.92,67.45,139.67
,104.68,231.95,88.17,53.4,150.87,128.76,203.38,143.96,151.59,152.11,41.85,151.15
,123.93,242.1,121.99,64.08,175.47,143.77,211.56,155.4,173.07,164.91,68.98,158.26
,137.21,207.27,130.76,34.56,188.92,105.07,242.85,109.74,187.67,127.43,80.29,127.14
,107.35,212.21,89.5,38.88,151.63,105.73,196.56,113.82,152.69,133.73,46.53,135.99),c(2,6,23))
qset2.dat<-aperm(qset2.dat,c(2,1,3))
qcet2.dat<-array(c(168.59,35.66,159.99,215.17,104.93,141.07,49.01,115.13,101.69,110.87,227.11,120.51
,165,38.35,163.89,214.32,108.26,144.79,50.17,124.02,103.5,113.13,220.43,122.65
,166.97,39.44,163.63,221.62,108.18,147.07,54.11,137.43,109.15,118.11,227.89,128.32
,164.02,44.22,171.76,237.29,95.93,161.44,32,131.85,95.91,126.54,222.54,133.51
,153.46,40.67,156.62,225.12,103.81,152.44,43.65,139.01,99.95,120.7,216.62,132.9
,148.3,52.66,147.61,239.16,90.37,164.39,32.39,145.13,88.1,130.85,209.22,145.88
,141.33,32.7,128.49,215.09,71.39,135.31,18.32,115.62,77.48,102.92,186.28,125.49
,130.21,18.71,136.38,201.17,68.85,125.19,11.43,94.64,64.22,90.83,184.64,107
,130,26.8,134.24,217.41,77.07,141.93,14.66,122.68,74.56,105.14,189.9,109.31
,134.99,22.44,116.08,205.87,74.69,130.05,22.18,96.09,68.08,95.91,185.97,115.99
,146.87,22.6,111.5,201.74,77.67,124.1,23.04,97.03,80.91,85.93,192.2,118.4
,146.26,23.38,119.94,209.46,75.94,125.31,19.72,100.92,82.97,87.47,192.22,109.63
,138.32,25.44,119.25,208.88,71.18,133.17,16.68,103.45,69.21,98.08,178.43,123.71
,142.57,19.25,99.1,197.59,62.7,111.06,3.94,86.78,69.17,79.06,180.48,117.38
,144.57,21.44,129,204.76,80.93,121.11,18.96,103.58,82.35,90.59,191.23,111.98
,137.6,19.71,120.16,214.71,77.7,128.06,16.08,102.22,73.54,90.19,193.21,118.14
,123.81,22.67,118.93,207.44,73.63,129.11,4.94,104.21,71.72,100.35,191.62,119.61
,131.1,28.64,94.82,211.86,59.44,129.47,3.89,102.71,70.28,95.52,178.71,131.69
,155.84,22.41,112.53,207.45,68.64,118.65,8.37,97.65,78.78,85.28,184.6,114.92
,174.04,27.44,108.13,202.49,80.74,111.51,20.1,79.49,96.21,80.84,205.25,127.23
,127.85,20.44,119.59,208.01,61.98,125.83,4.2,112.66,80.06,90.68,182.45,108.83
,146.48,28.36,113.82,214.92,75.14,127.12,14.06,98.64,81.57,96.68,198.33,125.63
,139.46,33.99,99.48,220.24,73.24,135.58,8.99,105,75.01,103.38,191.51,135.63
,152.27,30.86,138.67,226.92,91.78,137.46,24.93,121.51,92.28,107.03,216.37,131.07
,139.15,28.78,133.63,210.56,88.84,131.16,26.42,112.92,87.97,99.71,201.01,119.24
,153.32,20.49,109.35,201.5,76.58,109.87,10.91,91.13,74.45,85.03,183.14,119.64
,166.02,23.04,140.31,215.68,91.82,128.43,26.01,109.78,94.61,99.14,212.09,119.48
,127.79,22.26,136.79,202.76,89.87,131.73,26.63,113.72,88.85,90.81,198.88,105.03
,169.22,26.47,158.71,210.32,102.37,134.39,30.64,108.87,100.56,97.55,220.44,125.16
,146.88,23.83,116.78,218.55,84.03,132.22,25.08,111.84,89.31,94.31,194.7,129.03), 
c(2,6,30))
qcet2.dat<-aperm(qcet2.dat,c(2,1,3))
qlet2.dat<-array(c(134.26,224.28,122.34,36.86,174.35,111.24,235.79,127.9,174.49,142.43,51.5,141.89
,139.29,231.38,80.82,47.82,159.37,105.97,236.93,120.91,162.59,139.69,39.92,163.16
,151.57,219.95,104.42,49.26,162.95,105.68,241.12,123.23,181.45,133.49,60.47,151.83
,146.16,231.16,95.78,46.87,170.85,111.77,234.72,109.48,178.8,140.05,56.79,157.8
,150.16,222.81,81.59,53.08,161.66,102.76,238.7,94,173.74,131.99,42.42,166.01
,134.04,218.32,84.43,40.37,155.91,104.39,227.07,112.75,163.69,135.45,45.04,154
,141,221.6,98.86,46.69,160.43,112.77,218.15,114.61,162.83,136.17,37.06,153.86
,123.63,231.23,66.42,46.17,140.02,108.31,217.82,118,147.98,137.78,25.35,159.24
,137.62,226.64,96.5,39.64,164.59,105.39,235.04,102.87,169.37,134.5,46.83,150.8
,173.79,206.17,90.39,27.9,161.84,84.83,234.85,91.6,180.59,110.28,55.57,160.27
,117.39,233.75,114.47,42.63,167.1,124.53,229.45,135.47,167.66,150.42,40.92,142.56
,131.55,225.48,102.05,26.81,161.73,103.06,244.48,115.44,170.49,133.86,38.65,144.13
,134.78,226.28,110.56,41.28,170.17,112.5,231.45,114.1,176.5,139.26,51.28,145.88
,115.95,227.77,85.99,46.33,156.48,111.17,220.59,117.43,161.6,140.85,43.06,154.56
,133.09,226.95,99.38,40.42,160.22,111.1,236.36,117.74,168.23,138.63,44.79,149.37
,125.36,216.11,93.92,36.95,161.42,103.33,220.04,104.98,165.19,129.44,43.48,133.74
,123.1,202.86,99.27,42.82,157.22,98.93,206.46,111.18,162.05,129.3,64.39,123.39
,121.37,217.11,108.09,34.09,155.47,105.03,214.66,121.61,160.24,135.04,48.64,133.8
,120.54,232.1,85.37,53.23,157.81,109.98,213.66,117.86,169.01,146.01,38.88,151.63
,126.42,222.64,82.96,46.57,157.34,100.52,218.67,101.9,165.27,131.51,42.51,147.56
,126.91,220.11,105.38,35.76,158.15,109.6,207.12,117.23,160.11,138.72,40.17,142.75
,117.87,227.17,113.96,49.9,163.19,119.44,228.16,124.54,163.3,152.93,42.17,141.99
,125.23,224.48,93.4,39.47,166.22,109.39,234.83,108.47,165.65,139.89,42.06,144.54), 
c(2,6,23))
qlet2.dat<-aperm(qlet2.dat,c(2,1,3))
digit3.dat<-c(9,27,12,31,17,36,26,39,34,37,36,33,38,27,35,19,30,15,21,14,21,8,16,6,8,5
,17,40,21,38,26,36,27,32,25,28,22,27,19,29,24,25,26,20,28,16,26,13,18,14,15,17
,19,38,24,38,29,33,30,29,27,24,21,25,17,26,27,24,30,22,31,19,31,16,27,15,24,15
,9,40,15,43,24,41,29,36,24,30,20,26,12,22,20,22,24,20,21,16,18,14,13,12,9,10
,14,41,21,42,29,42,35,37,32,33,26,30,16,26,25,26,29,24,33,20,30,16,23,11,16,12
,24,39,28,40,35,38,38,35,34,30,29,27,22,24,27,24,29,22,31,19,28,15,20,11,13,12
,9,39,15,39,21,40,25,36,23,31,21,27,19,25,21,25,23,24,25,22,22,19,15,17,8,17
,8,38,14,41,25,43,29,38,25,33,18,29,8,28,12,27,16,25,18,23,13,21,7,21,1,22
,4,34,12,39,22,42,31,36,27,30,23,28,11,25,20,25,22,24,22,22,19,19,13,18,8,18
,21,36,25,37,31,36,33,32,32,28,29,25,27,22,29,21,31,20,31,18,28,16,24,16,20,16
,14,40,20,39,25,37,27,31,26,28,20,29,16,31,21,28,25,23,28,16,25,13,17,15,13,18
,12,40,20,42,30,42,36,33,31,24,23,22,16,23,25,22,31,18,33,13,31,9,24,8,17,8
,9,35,17,36,26,34,30,31,26,27,20,25,13,27,19,25,23,21,26,15,22,12,12,12,7,13
,17,38,24,39,30,37,34,34,31,28,22,25,16,28,21,26,27,24,30,20,26,15,18,14,10,17
,21,35,27,36,36,35,39,28,38,22,34,18,28,19,31,18,33,17,31,15,26,15,20,17,14,20
,16,40,20,43,25,39,27,31,24,24,19,21,17,23,19,22,21,21,23,21,22,18,19,16,15,16
,15,41,21,45,34,44,40,39,36,35,26,30,16,29,24,25,28,20,31,16,28,14,21,14,12,12
,11,42,22,42,32,39,35,34,32,29,25,26,20,27,25,26,31,23,35,19,31,14,21,12,16,15
,5,44,15,43,24,41,29,36,22,28,13,28,5,29,14,28,24,26,29,22,26,19,17,17,10,20
,14,37,19,39,25,38,28,32,25,26,20,22,14,23,17,23,21,20,23,17,21,15,16,15,11,15
,16,35,22,38,30,36,32,29,29,23,23,20,17,20,20,19,24,17,26,14,21,11,16,12,12,15
,14,38,17,40,25,42,28,38,27,32,24,28,20,25,23,25,26,24,28,21,24,18,18,17,10,18
,7,40,13,43,22,45,31,42,27,38,21,34,13,32,18,31,24,30,27,27,23,23,15,22,6,22
,14,35,21,36,26,34,31,30,28,26,25,22,21,18,21,17,22,16,23,15,20,12,13,10,5,10
,10,46,17,47,27,43,29,36,26,30,22,29,16,28,20,27,21,25,23,21,21,19,15,20,9,20
,18,39,24,42,33,41,38,35,37,30,32,28,28,27,33,22,37,18,41,15,37,13,29,11,21,12
,18,38,22,42,30,42,34,36,33,32,29,30,22,28,25,26,28,24,28,20,27,19,22,18,18,18
,9,41,17,43,30,40,34,31,30,23,23,19,11,19,15,17,18,13,21,10,17,8,12,7,5,7
,8,36,12,42,20,43,25,38,24,35,23,33,21,32,20,31,20,30,20,27,16,25,9,24,2,25
,19,41,24,45,33,45,38,38,36,31,28,27,21,23,24,22,26,20,28,17,26,14,20,13,14,11)
digit3.dat<-array(digit3.dat,c(2,13,30))
digit3.dat<-aperm(digit3.dat,c(2,1,3))
digit3.dat[,2,]<- -digit3.dat[,2,]
dna.dat<-array(c(23.825,12.021,13.002,25.742,18.923,13.225,22.784,25.153,15.445,17.418
,28.811,17.309,12.468,29.084,21.876,8.439,26.152,26.442,5.606,21.087
,29.894,6.200,14.944,33.112,8.585,12.034,38.280,14.445,8.903,39.890
,20.638,10.881,42.477,21.716,27.837,44.502,24.662,23.765,40.767,25.795
,17.985,37.076,24.130,12.686,33.111,20.240,10.380,28.222,14.758,7.779
,24.003,9.126,10.368,21.474,7.514,14.143,16.620,5.636,20.061,13.796
,9.181,24.591,9.975,15.219,27.372,7.792

,24.426,12.320,12.815,25.129,18.383,14.600,22.555,25.151,16.215,17.277
,28.975,17.620,12.170,28.852,22.133,7.779,25.623,26.239,5.167,20.176
,30.154,6.302,14.201,32.961,8.863,11.842,38.282,15.019,9.028,40.102
,21.051,10.202,42.586,22.535,27.544,45.013,24.570,23.693,40.375,25.299
,18.777,35.736,23.716,12.495,33.049,19.815,9.500,29.441,13.961,7.614
,25.021,9.187,10.748,21.303,7.266,14.440,16.882,5.328,20.518,13.417
,9.450,24.724,9.887,16.128,27.037,7.757

,24.824,12.355,12.721,25.459,18.787,14.482,22.545,24.678,16.606,17.068
,28.216,18.062,11.696,28.643,22.449,7.256,25.495,26.539,4.443,20.008
,29.999,5.302,14.541,33.584,8.475,12.111,38.538,14.696,9.392,40.333
,21.118,10.517,42.136,23.203,27.877,44.420,25.594,23.483,39.833,25.459
,18.596,35.894,23.239,12.499,33.283,19.753,10.123,28.758,14.685,7.357
,24.571,9.730,10.141,20.369,7.686,14.172,15.883,5.535,20.148,13.033
,9.722,25.009,10.287,16.126,27.525,8.623

,24.562,12.271,12.907,25.641,18.820,14.322,22.811,24.763,17.048,17.248
,28.619,18.489,11.519,27.989,22.598,8.070,25.157,27.032,3.894,19.924
,29.213,5.220,15.143,33.490,8.348,12.542,38.506,14.119,9.535,40.897
,20.276,10.987,42.136,23.192,28.071,44.539,25.274,23.775,39.509,25.630
,18.350,35.879,23.652,12.505,32.860,19.781,9.775,28.642,14.454,6.712
,24.583,9.440,9.756,20.462,7.865,14.225,16.493,5.837,20.340,13.680
,9.425,25.287,10.594,15.637,27.792,8.320

,24.498,12.655,13.535,26.012,19.043,14.582,22.839,24.695,17.030,17.673
,28.305,18.898,12.275,28.466,22.641,7.960,25.669,26.958,4.000,20.256
,29.476,5.319,15.623,33.821,8.755,12.476,38.043,14.771,9.263,39.758
,20.522,10.340,41.694,22.720,27.780,45.260,25.492,23.558,40.577,25.904
,18.294,35.791,23.664,13.116,32.898,19.581,9.743,29.357,14.351,7.149
,25.092,10.135,9.427,20.128,7.067,13.562,15.964,5.322,20.148,13.918
,9.684,24.717,10.987,15.626,27.045,8.706

,24.253,12.608,13.053,25.931,19.065,13.470,23.190,24.682,16.697,17.355
,28.266,18.178,11.941,28.445,21.624,8.302,26.165,26.289,4.645,20.902
,29.472,5.834,15.445,33.346,9.014,13.039,38.182,14.821,8.706,39.192
,20.416,10.356,42.138,21.676,27.355,45.520,25.249,22.900,41.575,26.118
,18.221,36.773,23.771,12.865,32.905,19.117,9.404,29.920,14.366,7.244
,25.218,9.684,9.810,21.100,7.052,13.603,16.361,5.871,20.211,14.759
,9.332,24.692,10.714,14.984,27.381,8.114

,23.794,12.783,13.221,26.116,19.069,14.260,22.944,24.986,16.310,17.113
,28.409,17.881,11.562,28.432,21.504,7.479,25.690,25.928,4.418,20.599
,30.049,5.580,15.584,34.067,9.173,12.544,38.631,15.060,9.018,38.984
,20.102,10.329,42.574,22.345,27.701,44.920,25.689,22.778,41.294,26.502
,17.802,36.716,23.573,12.346,33.446,19.038,9.376,29.982,14.359,7.379
,25.213,9.768,10.069,20.982,7.294,13.751,16.439,5.559,20.199,13.972
,9.213,24.830,10.359,15.581,27.031,7.719

,23.617,12.350,13.813,25.639,19.130,14.870,22.756,25.155,16.183,16.854
,28.676,17.591,11.233,28.507,21.359,7.934,26.050,26.376,5.116,20.641
,29.944,5.718,15.029,32.991,8.692,12.395,37.840,14.732,8.728,39.201
,19.511,9.906,43.019,22.376,27.740,44.419,25.920,23.168,41.030,26.208
,18.416,36.239,23.729,12.432,33.529,19.435,9.740,29.434,14.798,7.666
,24.759,9.937,10.295,20.754,7.938,13.711,15.998,5.671,20.047,14.131
,9.114,24.716,10.494,14.954,26.951,7.747

,24.012,12.503,13.953,25.912,19.342,14.449,23.029,25.100,16.549,17.484
,28.586,18.903,11.667,28.481,21.397,7.785,26.260,26.607,4.911,20.648
,29.762,5.753,15.350,33.300,8.908,12.434,37.994,14.578,8.736,39.552
,19.979,10.470,43.119,22.662,27.821,43.727,25.443,22.420,41.221,26.247
,17.848,36.494,24.167,12.379,32.929,19.582,9.488,29.216,14.702,7.711
,24.595,9.384,9.829,20.604,7.131,13.683,16.123,5.060,20.085,13.915
,8.908,25.200,11.222,14.847,27.059,8.911

,23.851,12.509,13.804,25.990,19.002,14.469,22.931,24.771,16.634,17.508
,28.373,18.286,11.915,28.212,21.270,7.770,25.930,26.098,5.343,20.386
,29.828,5.846,15.077,33.966,9.031,12.428,38.933,15.099,8.924,39.810
,20.518,10.602,43.163,21.689,28.245,42.375,25.115,22.774,40.248,25.516
,17.567,36.126,23.642,11.533,33.571,19.658,9.291,29.426,14.590,7.776
,24.159,9.043,9.992,20.933,7.292,14.245,16.582,5.694,20.509,13.985
,8.699,25.401,10.715,15.055,27.404,8.281

,23.695,12.829,13.504,25.764,19.437,13.969,22.597,25.450,15.948,17.092
,28.423,18.022,11.790,28.676,21.803,7.059,25.934,26.430,5.091,20.793
,30.234,5.589,15.473,34.066,8.997,11.769,38.160,14.696,8.522,40.353
,19.819,10.878,44.196,22.554,28.318,43.074,25.597,23.197,40.337,25.806
,17.780,36.375,23.787,11.741,33.154,19.557,8.962,29.841,14.286,7.574
,25.083,9.196,9.944,20.748,7.313,14.236,16.045,5.646,20.541,13.471
,9.889,25.682,11.896,15.361,27.634,8.380

,24.570,13.140,12.786,26.102,19.946,13.576,22.686,25.498,15.318,17.702
,28.165,18.723,12.523,28.503,22.326,7.233,26.220,26.025,4.907,21.266
,30.023,5.678,15.793,33.796,9.233,12.043,37.348,14.858,8.771,40.203
,19.898,10.903,44.171,21.444,28.378,42.556,25.142,23.227,39.449,25.636
,17.412,36.446,23.332,11.752,33.264,18.848,9.020,29.643,13.468,7.620
,24.944,9.220,9.848,21.160,7.691,13.422,15.760,6.713,19.640,13.470
,9.683,25.976,12.855,14.939,27.417,8.986

,24.627,13.046,12.380,25.895,19.776,12.833,23.067,25.155,15.595,18.103
,27.715,18.877,13.546,27.639,23.414,7.305,26.190,26.006,4.724,20.975
,29.586,5.140,15.664,33.846,9.195,11.925,37.119,15.034,8.789,39.473
,20.075,10.553,43.273,21.398,28.628,42.349,25.145,23.406,40.183,25.822
,18.112,36.110,23.461,12.655,32.655,18.873,9.201,29.785,13.440,7.458
,25.414,8.322,9.906,21.597,7.603,13.366,16.335,6.232,19.783,14.282
,9.705,25.919,12.853,14.983,27.636,9.035

,24.787,13.151,13.411,26.120,19.950,13.492,22.909,25.361,15.449,17.724
,28.036,18.727,13.426,27.879,23.848,7.001,25.988,26.517,5.155,20.872
,30.584,5.420,15.382,34.616,9.694,11.348,37.645,15.821,8.699,39.613
,20.186,10.669,43.251,21.809,28.473,42.674,25.263,23.078,39.879,25.499
,17.890,35.848,23.637,12.841,32.064,18.802,8.835,28.782,13.400,7.644
,24.509,8.059,10.687,21.628,7.038,13.879,16.268,6.444,19.486,12.936
,9.772,25.582,11.871,15.490,27.474,8.755

,24.879,12.862,13.363,25.855,19.829,14.067,22.454,25.488,15.988,16.932
,28.288,18.371,13.232,28.098,23.747,7.640,26.351,27.493,5.631,20.657
,30.620,5.432,14.891,34.055,10.077,11.331,37.364,15.912,8.673,39.810
,20.365,10.935,43.227,22.006,28.874,42.468,25.377,23.089,40.458,25.434
,18.097,36.265,23.154,12.725,32.439,18.303,8.809,29.316,12.224,7.154
,25.247,8.256,10.648,21.434,6.867,13.846,16.261,6.495,19.229,12.785
,9.444,25.660,11.445,16.015,27.106,9.083

,24.569,13.132,13.356,26.118,20.059,13.598,22.943,25.493,16.332,17.156
,28.220,18.641,12.640,27.993,23.201,6.940,26.026,26.789,5.822,20.299
,30.624,5.577,14.426,34.013,9.260,11.671,38.441,15.405,9.037,40.045
,20.732,11.190,43.389,22.255,28.629,41.898,25.166,23.030,39.331,25.211
,17.886,34.989,22.927,11.813,32.321,18.034,7.986,30.112,12.507,7.039
,25.429,7.760,10.850,22.086,6.699,14.098,16.608,6.578,19.622,12.735
,9.516,25.603,10.904,15.494,27.452,8.388

,24.374,13.368,13.547,25.968,20.247,13.246,23.008,25.503,15.953,17.295
,28.819,18.144,12.941,28.253,23.238,7.486,25.753,26.647,5.568,20.694
,30.577,5.546,14.753,34.301,9.336,11.829,38.770,14.971,8.889,41.320
,20.733,10.856,43.703,22.427,28.706,42.112,25.228,23.369,39.109,25.342
,17.710,35.477,23.173,11.708,32.903,18.173,7.920,30.467,12.415,7.036
,25.282,8.608,10.579,22.073,6.742,14.416,16.828,6.456,19.579,12.318
,9.521,25.792,11.092,15.795,27.217,8.544

,24.523,13.575,13.730,26.283,20.256,14.113,22.742,25.391,16.837,16.910
,28.361,18.931,12.169,28.452,22.750,7.157,25.331,26.758,5.413,20.606
,30.489,5.381,15.019,34.458,9.238,12.654,38.882,15.055,9.116,40.835
,20.612,11.073,44.104,22.835,28.368,41.865,25.681,23.191,38.900,25.820
,17.319,35.189,23.205,11.268,33.056,18.536,8.055,29.949,12.768,7.279
,25.385,7.923,10.343,21.699,7.260,13.681,16.484,6.793,20.112,13.511
,9.567,26.174,11.279,15.334,27.686,8.372

,24.006,13.405,12.592,26.309,20.061,13.619,23.040,25.554,16.139,17.822
,28.432,19.405,12.417,28.345,23.022,7.234,25.262,26.351,5.410,20.248
,30.270,5.563,14.632,34.272,9.276,12.834,39.126,14.824,9.049,40.659
,19.558,11.043,44.656,22.847,28.845,41.482,25.470,23.095,38.692,25.220
,16.932,36.035,22.990,11.096,33.874,18.029,7.900,30.415,12.609,6.544
,25.900,8.550,9.951,21.714,7.448,13.457,16.513,6.870,19.847,13.544
,8.872,25.764,11.318,14.641,27.537,7.839

,23.936,13.085,12.376,25.865,19.427,13.664,23.213,25.058,16.283,18.112
,28.607,19.261,12.743,28.401,22.584,7.267,25.067,26.234,6.259,20.058
,30.295,5.223,14.693,34.091,8.846,12.606,38.941,14.581,8.891,40.367
,19.148,10.977,44.551,22.080,28.672,41.808,25.625,23.181,39.061,25.517
,17.658,35.650,23.555,11.141,33.784,18.422,7.879,30.872,13.176,6.986
,25.093,9.595,10.402,21.787,7.283,14.052,17.155,6.272,19.831,13.550
,8.912,25.658,10.834,15.066,27.919,8.719

,24.228,13.221,11.812,25.933,19.620,13.289,23.382,25.319,15.563,18.068
,28.246,18.967,13.613,28.286,22.944,8.360,25.533,27.070,6.114,20.845
,30.299,5.328,14.990,34.186,8.986,11.480,38.107,14.578,9.012,40.742
,18.999,11.541,45.355,21.695,28.682,41.301,25.220,23.429,38.858,25.281
,17.872,35.522,23.677,11.484,33.448,19.107,8.155,29.984,12.684,7.155
,25.575,10.248,10.120,21.324,7.475,14.113,16.755,5.670,20.340,13.921
,9.459,25.472,10.983,15.771,28.089,9.318

,24.316,12.458,11.843,25.508,18.965,13.747,23.278,24.916,16.394,17.774
,28.402,19.187,12.378,28.542,22.605,7.252,25.365,26.306,6.069,20.368
,29.967,5.496,14.233,33.573,8.681,11.885,38.394,14.356,9.525,41.153
,19.296,11.200,44.990,21.956,28.698,41.790,25.477,23.507,38.898,25.134
,17.916,35.541,23.580,11.891,33.423,18.965,7.818,30.491,13.199,7.030
,25.775,9.665,9.849,21.557,7.725,13.950,16.492,6.011,20.248,14.478
,9.515,25.732,11.493,16.071,28.073,8.995

,24.069,12.514,11.576,26.193,19.281,13.560,23.049,25.025,15.906,17.225
,28.166,19.122,12.145,27.688,22.386,7.635,25.248,26.821,5.438,20.677
,30.447,5.051,15.003,34.500,9.197,12.486,38.468,14.699,9.325,40.616
,19.694,10.810,44.394,21.784,28.486,41.204,25.422,23.303,39.170,25.769
,17.606,36.130,23.125,11.743,33.599,18.675,8.250,30.057,12.399,6.446
,25.964,9.671,9.713,21.051,8.200,13.572,16.156,6.096,20.125,14.059
,9.421,25.580,11.507,15.633,28.168,9.245

,24.382,13.415,12.005,26.227,19.807,13.778,23.018,25.170,16.477,17.422
,28.183,18.648,12.599,27.936,22.634,7.930,25.547,26.937,5.015,20.995
,30.642,5.322,15.062,34.391,9.387,12.104,38.000,14.870,9.037,40.766
,20.116,11.114,44.270,22.390,28.993,40.814,25.331,22.929,39.238,25.429
,17.057,36.390,22.793,10.993,34.070,18.588,7.965,30.069,11.850,6.919
,26.147,9.683,9.913,20.736,7.799,13.653,15.823,5.725,20.032,14.331
,9.347,25.125,11.310,15.548,27.697,8.797

,25.071,13.359,12.219,25.922,20.228,13.700,22.969,25.646,16.479,17.218
,28.753,18.370,12.285,28.254,22.454,8.012,26.006,26.558,5.684,20.677
,30.618,5.143,14.471,34.316,8.841,12.183,39.055,14.522,8.799,41.346
,19.438,10.811,44.824,22.870,28.797,40.740,25.663,23.231,38.089,25.721
,17.236,35.440,22.782,11.278,33.869,18.262,8.303,30.135,11.771,6.803
,25.902,9.355,9.961,21.065,8.367,13.499,15.573,5.697,19.805,14.183
,8.787,25.293,11.196,14.928,27.345,8.810

,24.555,13.167,12.945,26.062,19.730,14.230,23.077,25.395,16.617,17.485
,28.522,18.766,12.427,28.226,22.678,7.346,25.153,26.424,5.141,20.778
,31.057,5.181,14.573,34.397,9.210,12.140,38.451,14.821,8.758,40.470
,20.120,10.521,44.182,22.452,28.993,40.984,25.474,23.382,38.705,25.796
,17.536,35.590,23.012,12.041,32.760,18.130,9.197,29.415,11.158,6.757
,26.335,9.712,9.912,20.939,7.884,13.871,15.872,5.667,20.538,14.465
,9.330,25.550,11.486,15.120,27.100,8.018

,25.137,12.870,12.482,26.450,19.304,14.192,22.839,25.036,16.498,17.098
,27.988,18.477,12.516,27.594,23.342,6.652,25.339,26.480,5.038,20.611
,30.998,5.596,14.238,33.760,9.367,12.383,38.383,14.898,8.935,40.587
,20.080,11.064,43.809,22.269,28.515,41.425,25.532,22.994,39.101,25.689
,17.577,35.290,23.238,12.103,32.263,18.131,8.696,29.280,11.572,6.737
,25.943,9.747,10.059,20.420,7.506,14.235,15.951,5.520,20.764,14.883
,9.477,25.594,12.102,14.816,27.279,8.341

,24.811,13.177,12.248,26.359,19.647,13.791,22.819,25.076,15.935,17.259
,28.137,18.816,12.031,27.817,22.524,7.139,25.285,26.476,5.233,20.844
,30.918,5.023,14.625,34.304,8.951,12.113,38.641,14.903,9.183,40.263
,19.489,10.837,44.214,22.943,28.849,41.934,25.785,23.491,38.818,25.458
,17.501,35.809,23.084,12.499,32.507,18.711,8.415,29.596,11.995,6.648
,25.674,9.782,9.834,20.682,7.330,14.037,16.282,5.679,20.535,14.148
,9.529,25.973,12.347,14.857,27.412,8.891

,25.461,13.475,12.459,26.240,19.803,14.240,22.439,25.184,16.246,17.018
,28.482,18.729,12.138,27.810,22.628,6.935,25.500,26.392,5.181,21.303
,30.774,5.269,15.405,34.596,9.057,12.253,38.675,14.733,9.098,40.950
,18.928,10.930,44.884,22.702,28.490,41.516,25.756,22.649,39.048,26.150
,16.724,36.124,23.328,11.649,32.925,18.644,8.006,30.218,12.586,6.889
,25.838,9.876,10.168,20.938,7.544,13.875,16.115,5.352,20.441,14.098
,9.147,25.602,11.691,14.869,27.495,8.845

,24.646,13.201,12.498,26.020,19.910,14.081,22.700,25.203,16.348,16.921
,28.147,18.532,12.219,27.992,22.139,7.297,25.197,26.330,5.477,21.167
,31.099,5.289,15.822,35.342,8.801,12.200,39.234,14.381,9.267,41.504
,19.133,11.253,44.763,22.661,28.604,41.983,25.373,23.002,39.097,25.724
,16.979,36.597,23.274,11.205,33.696,18.583,7.895,30.512,12.581,6.908
,26.393,10.461,9.886,20.595,7.139,14.317,16.469,5.308,21.069,14.978
,8.870,25.450,10.652,14.988,27.488,7.901)
,c(3,22,30))
dna.dat<-aperm(dna.dat,c(2,1,3))
#
#
#
macf.dat<-c(54.33203,24.10905,69.5
,141.80250,21.59643,69.5
,132.23880,62.78124,69.5
,88.22106,52.28123,69.5
,147.10890,26.61518,90.0
,107.42800,27.77578,101.0
,99.74427,46.85715,97.0
,58.35540,29.57223,67.0
,134.87930,26.38988,67.0
,120.57150,66.43083,67.0
,79.14922,52.19386,67.0
,138.17850,37.67144,86.0
,101.76780,34.47324,97.5
,90.83484,54.19082,92.0
,50.04349,15.22191,71.5
,139.76850,21.33820,71.5
,124.26710,63.79000,71.5
,80.94720,51.47673,71.5
,147.78980,32.26446,94.0
,102.32870,25.23133,105.5
,90.64163,47.04449,98.5
,41.93115,24.83244,70.5
,138.72930,22.35828,70.5
,122.68840,65.04043,70.5
,74.09913,53.69709,70.5
,142.63240,35.36625,92.0
,97.04822,33.37946,111.5
,89.51760,56.40900,96.5
,48.44877,35.75250,68.0
,134.68970,33.55962,68.0
,120.88830,73.77755,68.0
,77.94841,65.53058,68.0
,135.42950,42.75692,91.0
,101.92880,40.64505,99.5
,87.04530,59.04715,92.5
,44.05272,36.38397,70.0
,133.47320,45.07240,70.0
,114.88450,83.74655,70.0
,70.13158,71.71946,70.0
,139.94810,54.24338,87.0
,97.87708,47.19924,105.5
,80.04594,65.13947,96.5
,53.94042,32.50219,69.0
,136.19530,33.46588,69.0
,120.75410,74.12678,69.0
,78.24210,60.88469,69.0
,142.38210,44.35960,89.5
,100.75480,38.09180,101.0
,87.92361,55.80280,94.5
,45.11740,11.62884,68.5
,132.08950,17.75877,68.5
,112.66980,57.22899,68.5
,71.76939,47.64127,68.5
,136.93780,26.69818,88.5
,96.43125,21.78416,101.0
,84.39475,41.89694,95.0
,42.09966,16.11359,69.0
,131.92440,19.70687,69.0
,121.73400,57.71608,69.0
,72.94730,49.99466,69.0
,136.93340,26.92302,89.0
,96.84825,26.58288,103.0
,84.30907,48.96717,94.5)
macf.dat<-array(macf.dat,c(3,7,9))
macf.dat<-aperm(macf.dat,c(2,1,3))
macm.dat<-c(34.82811,16.50834,77.5
,138.91980,15.13858,77.5
,125.15760,58.60464,77.5
,72.28854,49.79207,77.5
,146.19080,22.68885,100.0
,99.30268,24.86908,117.0
,91.79910,46.49960,107.0
,40.40179,3.73932,73.0
,132.23560,7.56574,73.0
,114.63210,53.28955,73.0
,70.66502,33.57051,73.0
,139.58480,21.61227,90.5
,97.93692,8.49867,108.5
,79.90506,28.91153,100.5
,40.54510,9.51130,75.0
,136.61260,15.82863,75.0
,106.90960,63.82611,75.0
,76.19816,46.63517,75.0
,145.02210,30.40421,94.5
,101.72660,19.45746,113.5
,86.97967,43.98130,105.5
,21.11454,16.57673,75.0
,131.52700,23.12809,75.0
,109.44810,63.03707,75.0
,61.73774,53.69610,75.0
,135.91480,34.78890,101.5
,90.65395,25.30813,117.5
,75.66082,49.38123,105.5
,30.79976,19.21503,73.5
,134.92160,32.11148,73.5
,115.81510,69.88405,73.5
,67.15240,57.06633,73.5
,139.56950,44.82271,97.5
,95.38217,25.95223,112.0
,78.97741,47.89584,107.0
,18.88770,10.47136,74.5
,130.35790,12.40497,74.5
,114.85390,57.63774,74.5
,63.47649,48.13175,74.5
,138.25830,25.62929,97.5
,89.01810,18.95535,117.0
,75.67622,43.31009,104.5
,40.28789,14.90687,69.0
,134.29020,14.66977,69.0
,125.54870,56.83236,69.0
,75.68020,53.65364,69.0
,142.36350,24.22211,92.5
,99.44497,25.41932,106.0
,87.63929,45.45810,99.5
,25.38359,10.64805,72.5
,130.99770,9.63434,72.5
,118.65580,54.78021,72.5
,68.79280,48.67834,72.5
,139.77820,20.83856,97.0
,91.36346,16.29169,111.5
,75.55544,44.28398,101.0
,27.93545,5.21197,71
,130.98990,4.76235,71
,103.16230,51.99304,71
,70.59641,43.71388,71
,136.03820,13.90246,92
,91.75840,14.31955,109
,79.95213,35.70748,95)
macm.dat<-array(macm.dat,c(3,7,9))
macm.dat<-aperm(macm.dat,c(2,1,3))

sooty.dat<- c(-1426,-310.4167
,-1424,-160.4167
,-1117,320.5833
,-755,854.5833
,1238,1363.5833
,2330,471.5833
,1435,-748.4167
,771,-557.4167
,433,-395.4167
,-176,-299.4167
,-376,-290.4167
,-933,-248.4167
,-1000.20254,-1601.5969
,-1076.57007,-1266.4282
,-1124.65334,-635.6890
,-1193.94980,147.7853
,-61.16474,1895.7533
,1484.57069,2113.5422
,1649.32657,746.7048
,1212.33458,388.9087
,730.79486,166.1701
,156.62415,-383.9590
,-88.03479,-533.8656
,-689.07556,-1037.3256)
sooty.dat<-array(sooty.dat,c(2,12,2))
sooty.dat<-aperm(sooty.dat,c(2,1,3))



    schizo.dat<-c(  0.345632   , -0.0360314
    , -0.356301   ,  0.0234333
   ,  0.0119311     ,  0.17692
     ,  0.37789    ,  0.480402
    ,  0.719631     , -0.41189
    ,  0.397921    , -0.140558
    ,  0.351751    , -0.385748
    ,  0.333756    , -0.655051
    ,  0.032181    , -0.275235
    ,  0.112563    , -0.506533
    , -0.233126     , -0.28334
    , -0.667337   ,  0.0522613
    ,  0.188945    , -0.142714
    ,  0.237198    ,  0.048306
    , -0.340236   ,  0.0997385
   , -0.0161814    ,  0.201017
    ,  0.301584    ,  0.516546
    ,  0.510795    , -0.323537
    ,  0.269407   , -0.0562209
    ,  0.239301    , -0.253218
    ,  0.229339     , -0.48236
   ,  0.0201328    , -0.122625
    ,  0.048306    , -0.341874
    , -0.200997    , -0.122698
     , -0.62316    ,  0.120534
    ,  0.124688   , -0.0262477
    ,  0.341616    ,  0.048306
    , -0.408509    ,  0.119819
  , -0.00814926    ,  0.277322
     ,  0.39797     ,  0.62498
    ,  0.591116    , -0.299441
     ,  0.35776   , -0.0361405
     ,  0.27143    , -0.249202
    ,  0.273516    , -0.530553
    ,  0.032181    , -0.110577
   ,  0.0724024     , -0.34589
    , -0.188949    , -0.138762
    , -0.707498    ,  0.140615
    ,  0.124688   , -0.0262477
    ,  0.329567     , -0.10832
    , -0.359859   ,  0.0341931
    , -0.056342    ,  0.152824
    ,  0.377668    ,  0.474464
    ,  0.673363    , -0.462009
    ,  0.373825    , -0.228912
    ,  0.323639    , -0.450005
     ,  0.34179    , -0.735372
   ,  0.0924219    , -0.295315
    ,  0.100555    , -0.540522
    , -0.181195    , -0.260518
    , -0.651361   ,  0.0963372
    ,  0.174798    , -0.178741
    ,  0.193021   ,  0.0683864
    , -0.539516    ,  0.210918
     , -0.16074    ,  0.333555
    ,  0.293246    ,  0.639288
    ,  0.520753    , -0.301366
    ,  0.245311   , -0.0281084
    ,  0.142916     , -0.28133
    ,  0.128938    , -0.510473
    , -0.152558    , -0.118609
   , -0.0761516    , -0.379879
    , -0.414127    , -0.127988
    , -0.848201    ,  0.246974
  , -0.00592517  , -0.00203414
    ,  0.337599    , -0.076192
    , -0.356431    ,  0.146356
    ,  0.068156    ,  0.201017
    ,  0.520571    ,  0.464342
    ,  0.601074    , -0.518234
    ,  0.321616    , -0.208831
    ,  0.251349    , -0.409844
    ,  0.197211    , -0.695211
  , -0.00797966    , -0.251139
  , -0.00787855    , -0.500361
     , -0.26726    , -0.201009
    , -0.719602     ,  0.25526
    ,  0.110557    , -0.140611
    ,  0.223261    ,  0.228767
    , -0.431293    ,  0.230133
   , -0.0439952    ,  0.389752
     ,  0.23941    ,  0.819652
    ,  0.641234    , -0.088515
    ,  0.377841    ,  0.140566
    ,  0.339703     , -0.10864
    ,  0.393998    , -0.345813
   ,  0.0201328   , -0.0583678
    ,  0.172844    , -0.283494
    , -0.268058   , -0.0890055
    , -0.777652    ,  0.251188
    ,  0.154806     ,  0.12832
    ,  0.243341    , -0.052358
    , -0.451374    ,  0.145795
   , -0.0480112    ,  0.297382
    ,  0.408185    ,  0.574616
     ,  0.56291    , -0.538339
     ,  0.31561    , -0.190749
    ,  0.251349    , -0.437957
    ,  0.169098    , -0.715291
  , .00005         , -0.194914
 ,  0.000153631    , -0.492329
    , -0.278085    , -0.117118
    , -0.797732    ,  0.255204
    ,  0.122674   , -0.0905492
    ,  0.287518   ,  0.0600918
    , -0.503583    ,  0.121699
   , -0.0861578    ,  0.307424
    ,  0.355977    ,  0.614776
    ,  0.603071    , -0.277295
    ,  0.303562   , -0.0100258
    ,  0.259382    , -0.293379
    ,  0.241387    , -0.538584
   , -0.0521564    , -0.114593
   ,  0.0443303    , -0.343735
    , -0.274068   , -0.0930217
    , -0.805764    ,  0.138738
     ,  0.12669  ,  0.00182023
    ,  0.319646     ,  0.20467
    , -0.431294    ,  0.153827
   , -0.0660775    ,  0.391762
     ,  0.29172    ,  0.793493
    ,  0.771746  , -0.00420256
    ,  0.407979    ,  0.166681
    ,  0.387896     , -0.12872
    ,  0.393998    , -0.394006
   ,  0.0864079   , -0.0382915
    ,  0.168828    , -0.279477
    , -0.225876    , -0.109086
    , -0.779633    ,  0.126677
    ,  0.190947    ,  0.126318
    ,  0.303582    ,  0.208686
    , -0.395149    ,  0.302422
  , -0.00985258    ,  0.387746
    ,  0.388105     ,  0.72522
    ,  0.611103    , -0.289343
    ,  0.327658   ,  0.0743116
    ,  0.283478    , -0.229121
    ,  0.201226    , -0.494408
  ,  0.00608663   , -0.0141951
   ,  0.0362983    , -0.279477
     , -0.24194  , -0.00466831
    , -0.699312    ,  0.339529
     ,  0.12669   ,  0.0901736
    ,  0.287518   , -0.0362937
    , -0.428786    ,  0.133507
   , -0.0801519    ,  0.233129
    ,  0.374063    ,  0.614792
    ,  0.671344    , -0.478098
    ,  0.383883    , -0.134524
    ,  0.315606    , -0.417876
    ,  0.265483    , -0.675131
   ,  0.0342193     , -0.19491
   ,  0.0443223    , -0.466226
    , -0.266036    , -0.201455
    , -0.735453    ,  0.202987
    ,  0.176577   , -0.0906634
    ,  0.251373    ,  0.212702
     , -0.42477    ,  0.149571
    , -0.124329    ,  0.345579
    ,  0.145148    ,  0.767403
    ,  0.671344  , -0.00821852
    ,  0.424044    ,  0.138569
    ,  0.407976    , -0.116671
     ,  0.44219    , -0.377942
    ,  0.122573   , -0.0583636
    ,  0.202939    , -0.283647
     , -0.14957   , -0.0970378
    , -0.755533    ,  0.114633
    ,  0.180593   ,  0.0940755
     ,  0.29555   , -0.0362937
    , -0.348465 ,  0.000976592
   , -0.0480233    ,  0.148792
    ,  0.273662    ,  0.458166
     ,  0.65528    , -0.341552
    ,  0.371835    , -0.106411
    ,  0.327655    , -0.377715
    ,  0.333757     , -0.65505
   ,  0.0984763    , -0.247119
    ,  0.118601    , -0.500514
    , -0.181699    , -0.253664
    , -0.675212   ,  0.0182479
    ,  0.192642    , -0.122792
    ,  0.309487    ,  0.265173
      , -0.3438    ,  0.289185
   ,  0.0199632    ,  0.482141
    ,  0.355444    ,  0.735605
    ,  0.613122    , -0.156788
    ,  0.357761    ,  0.128518
    ,  0.315607    , -0.132736
    ,  0.305644     , -0.36991
   ,  0.0121007   ,  0.0380178
    ,  0.100555    , -0.223253
    , -0.182924   ,  0.0134128
    , -0.623245    ,  0.301185
    ,  0.162814    ,  0.142714
    ,  0.193021    ,  0.305334
    , -0.472766    ,  0.163996
    , -0.172808     ,  0.38174
   ,  0.0847168    ,  0.829799
    ,  0.538908    ,  0.114214
    ,  0.249327    ,  0.273096
    ,  0.199141   ,  0.0319229
    ,  0.273516    , -0.165091
   , -0.0360921   ,  0.0340017
   ,  0.0563381    , -0.141071
    , -0.257222   , -0.0263121
    , -0.791835   ,  0.0884059
  , -0.00382644    ,  0.178572
     ,  0.39784    , -0.011935
    , -0.440185    ,  0.132558
  , -0.00814926    ,  0.309451
    ,  0.470037    ,  0.542738
    ,  0.580994    , -0.385704
    ,  0.389889    , -0.212847
    ,  0.283478    , -0.401812
    ,  0.237372    , -0.675131
   ,  0.0201328    , -0.267203
   ,  0.0644106    , -0.508393
     , -0.23742    , -0.196261
    , -0.711598    ,  0.184719
   ,  0.0944764    , -0.130548
    ,  0.269326   ,  0.0844506
    , -0.327735    ,  0.048221
   , -0.0242135    ,  0.253226
    ,  0.257186    ,  0.643139
    ,  0.601074      , -0.1849
    ,  0.325632   ,  0.0160683
    ,  0.307574    , -0.217073
    ,  0.301629    , -0.442199
   ,  0.0803737    , -0.126641
    ,  0.112603    , -0.351767
    , -0.153083     , -0.13602
    , -0.635293   ,  0.0280923
    ,  0.158733  , -0.00605034
    ,  0.217118    ,  0.212965
    , -0.416089    ,  0.192799
    , -0.104535    ,  0.353627
    ,  0.229073    ,  0.727476
    ,  0.601074      , -0.1849
    ,  0.293503    ,  0.112454
    ,  0.267414    , -0.112655
    ,  0.285565    , -0.402039
   ,  0.0281649   , -0.0101749
   ,  0.0844909    , -0.259397
    , -0.221356   , -0.0275863
    , -0.743727    ,  0.188735
   ,  0.0944764    ,  0.118448
    ,  0.293423    , -0.076192
    , -0.327736   ,  0.0442049
   , -0.0121653    ,  0.160856
    ,  0.357587    ,  0.394143
    ,  0.572961    , -0.457993
    ,  0.281455    , -0.228912
    ,  0.263398    , -0.470085
    ,  0.189179    , -0.679147
  , .00005         , -0.259171
     ,  0.02425    , -0.524458
    , -0.193244    , -0.208309
    , -0.583084   ,  0.0883332
    ,  0.142669    , -0.178741
    ,  0.289407   ,  0.0201935
    , -0.351832    ,  0.164687
  ,  0.00389893    ,  0.265274
     ,  0.40578    ,  0.510609
     ,  0.49264    , -0.433896
      ,  0.3176    , -0.120478
    ,  0.235285    , -0.329523
    ,  0.201227    , -0.554649
  ,  0.00406852     , -0.17885
   ,  0.0081857    , -0.403976
    , -0.245452    , -0.107908
    , -0.647341    ,  0.228896
   ,  0.0583318   , -0.0622752
    ,  0.301455   ,  0.0603542
    , -0.391993    ,  0.224928
   , -0.0121653    ,  0.341579
    ,  0.441925    ,  0.594946
     ,  0.45248    , -0.389719
     ,  0.29752   , -0.0963815
     ,  0.17906     , -0.30141
    ,  0.140986    , -0.570714
     , -0.02806    , -0.106561
   , -0.0560713    , -0.391927
    , -0.245452   , -0.0556988
    , -0.683486    ,  0.293153
   ,  0.0583318  , -0.00605033
    ,  0.237198   , -0.0922563
    , -0.396005  ,  0.00807245
    , -0.116583     ,  0.11668
     ,  0.29333    ,  0.402175
     ,  0.47256    , -0.502169
     ,  0.22523     , -0.22088
     ,  0.17906    , -0.433941
    ,  0.157051    , -0.634971
   , -0.0320761    , -0.271219
   , -0.0400071    , -0.500361
    , -0.257501    , -0.260518
     , -0.69955   ,  0.0521888
   ,  0.0382515    , -0.158661
    ,  0.321535   , -0.0641438
    , -0.379941    ,  0.108474
   , -0.0201975    ,  0.192985
    ,  0.445941    ,  0.450368
     ,  0.49264    , -0.506185
    ,  0.325632    , -0.184735
    ,  0.203157    , -0.365667
    ,  0.124922    , -0.675131
   , -0.0119958    , -0.206962
   , -0.0279589    , -0.476265
    , -0.285613    , -0.192245
    , -0.651357    ,  0.176687
   ,  0.0944764    , -0.126532
    ,  0.317519   , -0.0761919
    , -0.219297   , -0.0300639
    ,  0.072172    ,  0.140776
     ,  0.40578    ,  0.430287
    ,  0.580994    , -0.457993
    ,  0.345712    , -0.216863
    ,  0.287494    , -0.417876
      ,  0.2695    , -0.663083
    ,  0.096438    , -0.259171
   ,  0.0764588    , -0.496345
    , -0.165131     , -0.24847
    , -0.518831  , -0.00808046
    ,  0.162749    , -0.178741
    ,  0.363819    , -0.132663
    , -0.295482    ,  0.140679
    ,  0.178874    ,  0.174894
    ,  0.568759     ,  0.29153
    ,  0.593042    , -0.550362
    ,  0.345712    , -0.297185
    ,  0.267414     , -0.50623
    ,  0.136969     , -0.74742
   , -0.0039636    , -0.275235
   ,  0.0724427    , -0.520442
    , -0.239946    , -0.153262
   , -0.524623    ,  0.229146
    ,  0.194967    , -0.196981
    ,  0.255389    ,  0.160493
    , -0.435309    ,  0.214068
    , -0.124316    ,  0.357623
    ,  0.267522     ,  0.65901
    ,  0.597058    , -0.136708
    ,  0.305552   ,  0.0763094
    ,  0.255365    , -0.140768
    ,  0.257451    , -0.402038
    ,  0.032181   , -0.0382875
    ,  0.084491    , -0.279478
    , -0.243962   , -0.0167163
    , -0.761587    ,  0.231107
   ,  0.0945651   ,  0.0881595
    ,  0.279486  ,  0.00788297
    , -0.330892    ,  0.141779
   ,  0.0222759    ,  0.267264
    ,  0.386119    ,  0.484248
    ,  0.574959    , -0.365648
     ,  0.33569    , -0.102395
    ,  0.263398    , -0.345587
    ,  0.233355    , -0.566697
   , -0.0180098    , -0.178854
   ,  0.0965392    , -0.387911
    , -0.233908     , -0.12515
    , -0.627023    ,  0.206999
   ,  0.0681437   , -0.0866473)


schizo.dat<-array(schizo.dat,c(2,13,28))
schizo.dat<-aperm(schizo.dat,c(2,1,3))
braincon.dat<-schizo.dat[,,1:14]
brainscz.dat<-schizo.dat[,,15:28]




panf.dat<-array(c(47,-23,0,0,0,32,12,87,21,156,31,133,66,92,83,20,
63,-22,0,0,0,29,6,89,14,157,23,136,62,101,95,24,
56,-11,0,0,0,31,2,89,4,159,17,141,54,107,86,34,
51,-27,0,0,0,29,12,89,30,156,39,135,71,94,86,18,
51,-19,0,0,0,35,8,97,23,169,36,143,67,101,85,25,
54,-23,0,0,0,34,14,84,36,155,44,133,76,87,89,19,
53,-21,0,0,0,30,5,90,24,162,31,139,63,99,86,21,
56,-23,0,0,0,33,16,92,30,156,40,131,75,95,95,15,
55,-20,0,0,0,33,12,89,13,157,31,136,63,97,91,22,
35,-23,0,0,0,26,9,81,11,153,26,131,64,92,84,18,
48,-25,0,0,0,30,23,89,44,160,49,139,81,95,95,18,
35,-34,0,0,0,30,10,84,23,153,36,128,68,95,94,13,
46,-23,0,0,0,28,4,92,14,163,27,137,60,101,86,27,
42,-23,0,0,0,30,19,88,33,163,45,130,77,95,93,20,
50,-19,0,0,0,32,7,90,14,157,25,139,68,101,87,23,
54,-19,0,0,0,29,6,92,9,163,20,140,64,104,92,26,
46,-31,0,0,0,25,3,89,2,167,24,136,63,95,85,19,
47,-19,0,0,0,29,7,84,6,148,22,126,58,89,80,24,
46,-27,0,0,0,31,12,86,23,156,35,130,72,93,89,18,
49,-23,0,0,0,29,14,88,25,159,36,133,76,90,90,16,
50,-23,0,0,0,32,9,92,14,167,31,141,71,97,87,28,
40,-23,0,0,0,30,13,92,30,167,40,140,74,99,91,21,
32,-25,0,0,0,36,7,96,12,170,22,147,61,108,81,35,
41,-30,0,0,0,29,18,87,31,160,43,135,74,87,89,10,
46,-25,0,0,0,29,15,88,23,163,37,134,70,91,85,22,
43,-22,0,0,0,33,-1,96,10,178,24,153,66,105,87,32),c(2,8,26))
panf.dat<-aperm(panf.dat,c(2,1,3))
select<-c(5,1,2,3,4,6,7,8)
panf.dat<-panf.dat[select,,]

panm.dat<-array(c(43,-21,0,0,0,34,14,101,25,179,40,150,75,104,90,31,
48,-23,0,0,0,31,5,92,11,166,24,144,63,99,82,24,
43,-23,0,0,0,29,13,92,21,161,33,138,68,100,84,21,
45,-32,0,0,0,30,8,100,14,163,28,143,74,102,97,21,
40,-27,0,0,0,29,7,93,12,166,23,147,69,102,90,25,
49,-19,0,0,0,33,3,94,11,165,23,144,64,108,89,30,
55,-17,0,0,0,31,6,97,17,168,29,144,69,101,85,23,
49,-27,0,0,0,26,11,92,16,178,30,152,78,100,89,23,
48,-23,0,0,0,29,10,96,32,166,42,139,69,96,87,21,
49,-19,0,0,0,29,3,100,14,172,25,152,63,109,87,32,
49,-26,0,0,0,34,8,91,22,168,34,140,74,93,93,17,
52,-26,0,0,0,32,7,92,10,172,28,143,66,100,93,23,
36,-26,0,0,0,28,11,92,25,165,34,140,71,98,90,18,
46,-26,0,0,0,33,12,94,20,174,35,145,71,106,93,24,
47,-22,0,0,0,31,10,88,29,156,34,138,68,93,91,20,
47,-25,0,0,0,30,0,97,2,172,19,147,62,102,82,26,
53,-19,0,0,0,29,6,80,2,142,13,123,54,100,91,31,
49,-24,0,0,0,29,5,90,14,168,31,144,77,96,89,22,
52,-25,0,0,0,30,6,93,15,167,25,147,74,99,92,24,
42,-27,0,0,0,32,11,87,26,159,37,137,76,95,87,21,
37,-25,0,0,0,29,-2,87,-2,175,16,152,63,107,91,33,
41,-27,0,0,0,25,6,87,13,167,29,141,70,99,89,28,
46,-24,0,0,0,35,12,98,15,175,31,152,76,106,92,27,
45,-22,0,0,0,29,8,90,12,176,24,151,70,104,88,24,
46,-27,0,0,0,29,2,88,6,166,23,138,67,100,86,24,
43,-20,0,0,0,33,7,94,17,165,30,143,70,93,81,24,
44,-20,0,0,0,29,1,87,0,154,19,133,64,98,80,16,
52,-28,0,0,0,28,2,84,18,155,25,135,63,97,93,21),c(2,8,28))
panm.dat<-aperm(panm.dat,c(2,1,3))
select<-c(5,1,2,3,4,6,7,8)
panm.dat<-panm.dat[select,,]


pongof.dat<-array(c(43,-31,0,0,0,29,13,90,45,150,48,126,74,80,91,19,
49,-31,0,0,0,33,28,93,70,152,72,130,86,78,85,-5,
51,-36,0,0,0,32,34,100,74,152,74,131,87,73,92,-2,
48,-26,0,0,0,32,10,89,35,154,40,136,65,88,86,14,
56,-29,0,0,0,24,11,87,44,155,48,133,72,82,91,12,
49,-30,0,0,0,29,25,96,72,159,69,137,85,78,93,9,
51,-28,0,0,0,26,11,94,50,162,49,132,75,87,89,13,
57,-24,0,0,0,35,10,98,55,164,53,138,77,88,93,21,
37,-29,0,0,0,38,30,93,63,161,68,129,86,68,87,6,
53,-30,0,0,0,29,5,88,39,147,38,131,64,84,90,14,
58,-24,0,0,0,30,0,88,14,151,24,133,59,91,88,21,
54,-26,0,0,0,36,11,106,41,173,48,146,78,94,99,24,
59,-25,0,0,0,29,4,102,28,177,34,156,63,100,90,25,
32,-36,0,0,0,25,26,90,71,144,69,124,82,72,91,3,
52,-30,0,0,0,35,21,99,55,164,59,143,79,90,95,22,
51,-27,0,0,0,35,11,92,37,152,42,132,69,88,95,26,
47,-24,0,0,0,35,7,98,36,163,38,138,66,89,87,21,
60,-23,0,0,0,23,-2,82,28,158,32,135,56,90,89,25,
46,-31,0,0,0,25,4,87,34,145,37,120,66,80,89,20,
46,-28,0,0,0,36,29,94,73,148,71,123,82,74,92,11,
32,-37,0,0,0,36,32,88,81,140,81,117,91,63,90,-5,
43,-27,0,0,0,25,2,90,32,159,37,131,62,91,87,22,
38,-27,0,0,0,30,4,93,30,160,36,136,60,92,86,27,
38,-27,0,0,0,34,14,92,53,155,48,129,71,82,84,24),c(2,8,24))
pongof.dat<-aperm(pongof.dat,c(2,1,3))
select<-c(5,1,2,3,4,6,7,8)
pongof.dat<-pongof.dat[select,,]


pongom.dat<-array(c(49,-45,0,0,0,26,25,106,68,190,68,151,84,80,100,14,
64,-28,0,0,0,31,10,106,46,185,53,156,77,95,99,14,
55,-31,0,0,0,33,23,113,72,186,72,160,92,95,102,21,
64,-25,0,0,0,36,5,109,35,188,42,165,69,102,101,33,
46,-39,0,0,0,31,36,106,97,155,89,133,95,74,92,1,
53,-28,0,0,0,33,17,111,54,185,59,159,87,88,98,16,
47,-36,0,0,0,36,35,114,72,183,74,150,94,81,105,15,
44,-35,0,0,0,37,15,110,69,183,74,153,84,89,94,11,
55,-43,0,0,0,26,24,105,71,172,75,146,91,81,104,1,
49,-33,0,0,0,35,11,113,58,188,56,159,79,92,93,24,
45,-32,0,0,0,29,20,107,67,184,67,155,85,82,93,11,
48,-34,0,0,0,32,33,116,89,192,91,160,100,81,99,5,
41,-51,0,0,0,29,50,100,127,139,115,119,105,53,96,-19,
60,-45,0,0,0,32,36,102,108,163,97,129,98,75,98,-6,
65,-35,0,0,0,30,24,112,65,188,72,158,88,87,103,18,
54,-28,0,0,0,33,7,114,33,206,42,172,74,98,94,28,
41,-39,0,0,0,34,24,115,79,188,79,154,93,79,99,2,
42,-40,0,0,0,39,41,114,112,187,102,147,112,70,97,-14,
65,-27,0,0,0,30,17,109,62,187,66,151,83,82,96,18,
54,-36,0,0,0,30,29,122,65,204,70,176,93,87,96,9,
55,-37,0,0,0,32,25,116,75,190,71,155,88,84,98,8,
50,-35,0,0,0,26,8,101,39,172,49,148,75,87,98,18,
78,-31,0,0,0,28,15,119,56,204,60,179,91,94,103,8,
42,-32,0,0,0,34,37,117,102,181,99,148,102,83,97,5,
52,-39,0,0,0,38,15,111,58,201,63,158,89,91,95,-8,
47,-37,0,0,0,23,37,98,85,160,83,136,89,75,94,-5,
49,-37,0,0,0,34,37,115,105,179,98,151,96,80,97,5,
48,-32,0,0,0,36,10,113,53,189,57,166,79,100,99,21,
41,-24,0,0,0,39,4,128,42,209,47,178,73,102,93,21,
50,-32,0,0,0,39,27,121,73,198,75,166,95,89,101,15),c(2,8,30))
pongom.dat<-aperm(pongom.dat,c(2,1,3))
select<-c(5,1,2,3,4,6,7,8)
pongom.dat<-pongom.dat[select,,]



#The Procrustes routines in the next part were initially
# written by Mohammad Faghihi (University of Leeds) 1993, although many improvements, corrections,
# and speed-ups have been done since then.



# add(a3) compute the summation of a3[,,i]'s

# bgpa(a3) compute the scaling coefficients (bi's)

# close1(a) adds one additional row to matrix a that is the same as the first row

# cnt3(a3) replace each a3[ , , i] by fcnt(a3[ , , i])

# del(po, w1) plots point of po and joins them by contiguity matrix w1.

# dif(a3) compute sum( tr (xi-xj)'(xi-xj) )/n^2 for i<j and ( xi=a3[ , , i] )
#NB RETURN$Gpa is  now GSS/n  (CHANGED to this in version 0.92!)

# dis(a, b, c) compute distances between three points a, b and c. Each 
#  point should contain two co-ordinates

# fJ(n) function makes a nxn matrix as (I - (1/n) * J(n)) such that J(n) is a nxn  
#  matrix with all entries 1.

# fcel(n,d) generates n points such that the triangles between them have 
#  equal edges

# fcnt(a) = fJ(no. of rows of matrix "a")*a

# fgpa(a3,tol1,tol2) compute rotated and scaled shapes, Gpa statistics (dif.), and number of
#  iterations.

# fopa(a,b) function computes the ordinary procrustes statistics for two matrices  a
#  and b.

# fort(a,b) computes an orthogonal matrix to rotate matrix b such that the  
#  difference between a and b be minimum. 

# fos(a,b)  computes a scalar to scale matrix b such that the difference between a 
#  and b be minimum.

# ftrsq(a,b) tr{(b'aa'b)^(1/2)}
 
# graf(a3) plot a3[ , , i] for all i's in one plot

# msh(a3) compute the mean shape of a3[ , , i]'s

# Enorm(a) compute || a ||=sqrt{ trace( x'x ) }

# rgpa(a3,p) find the new rotated data till dif(old)-dif(new)<p

# sgpa(a3) scaling a3 by bgpa coefficients

# sh(a) compute the co-ordinates of shape point for triangle a. a should be a 3x2 
#  matrix.

# sim1(n, d, s) simulated n points with normal distribution 
#  (mean=fcel(n,d) sd=s )

# vec1(a3) vectorizes matrices a3[ , , i]


#=========================================================================


add<-function(a3)
{
	s <- 0
	for(i in 1:dim(a3)[3]) {
		s <- s + a3[,  , i]
	}
	return(s)
}


bgpa<-function(a3,proc.output=FALSE)
{
#assumes a3 is centred
	h <- 0
#	zd <- cnt3(a3)
	zd<-a3

	s <- 0
n<-dim(a3)[3]
#	for(j in 1:dim(a3)[3]) {
#		s <- s + (Enorm(zd[,  , j])^2)
#	}
	aa<-apply(zd,c(3),Enorm)^2
	s<-sum(aa)
#	for(i in 1:dim(a3)[3]) {
#		h[i] <- sqrt(s/(Enorm(zd[,  , i])^2)) * eigen(zz)$vectors[i, 1]
#	}


#try to speed it up!

      omat<-t(vec1(zd))
 kk<-dim(omat)[2]
 nn<-dim(omat)[1]
     if (nn > kk){
 
#      qq<-diag(cov(vec1(zd)))
       qq<-rep(0,times=nn)     

      for (i in 1:n){
      qq[i]<-var(omat[i,])*(n-1)/n
      omat[i,]<-omat[i,]-mean(omat[i,])
      }
omat<-diag(sqrt(1/qq))%*%omat
n<-kk
Lmat<-t(omat)%*%omat/n


eig<-eigen(Lmat,symmetric=TRUE)
U<-eig$vectors
lambda<-eig$values

V<-omat%*%U

vv<-rep(0,times=n)
for (i in 1:n){
vv[i]<-sqrt(t(V[,i])%*%V[,i])
V[,i]<-V[,i]/vv[i]
}

delta<-sqrt(abs(lambda/n))*vv
od<-order(delta,decreasing=TRUE)
delta<-delta[od]
V<-V[,od]

	h<-sqrt(s/aa)*V[,1]
}
if (kk>=nn){
	zz <- cor(vec1(zd))
	h<-sqrt(s/aa)*eigen(zz)$vectors[, 1]
}

h<-abs(h)
	return(h)
}

close1<-function(a)
{
	a1 <- matrix(0:0, nrow = dim(a)[1] + 1, ncol = dim(a)[2])
	for(i in 1:dim(a)[1]) {
		a1[i,  ] <- a[i,  ]
	}
	a1[dim(a)[1] + 1,  ] <- a[1,  ]
	a1
}


cnt3<-function(a3)
{
	#zz <- array(c(0:0), dim = c(dim(a3)[1], dim(a3)[2], dim(a3)[3]))
	#for(i in 1:dim(a3)[3]) {
		#zz[,  , i] <- fcnt(a3[,  , i])
	#}
	
	zz<-apply(a3,3,fcnt)
	zz<-array(zz,dim(a3))
	return(zz)
}

 del<-function(po, w1)
{
	plot(po, type = "n", xlab = "x", ylab = "y")
	text(po)
	n <- dim(po)[1]
	for(i in 1:n) {
		for(j in i:n) {
			if(w1[i, j] > 0) {
				a1 <- c(po[i, 1], po[j, 1])
				b1 <- c(po[i, 2], po[j, 2])
				lines(a1, b1)
			}
		}
	}
}

dis<-function(a, b, c)
{
	d <- 0
	d[1] <- sqrt((a[1] - b[1])^2 + (a[2] - b[2])^2)
	d[2] <- sqrt((a[1] - c[1])^2 + (a[2] - c[2])^2)
	d[3] <- sqrt((c[1] - b[1])^2 + (c[2] - b[2])^2)
	d
}

dif.old<-function(a3)
{
        s <- 0
        for(i in 1:(dim(a3)[3] - 1)) {
                for(j in (i + 1):dim(a3)[3]) {
                        s <- s + ((Enorm(a3[,  , i] - a3[,  , j]))^2)
                }
        }
        return(s)
}


#dif<-function(a3)
#original (slow) version
#{
#        s <- 0
#n<-dim(a3)[3]
#mshape<-add(a3)/n
#psum<-0
#for (i in 1:n){
#x<-a3[,,i]-mshape
#psum<-psum+sum(diag(t(x)%*%x))
#}
#psum*n
#}

#dif<-function(a3)
##faster version
#{
#x<-sweep(a3,c(1,2),apply(a3,c(1,2),mean))
#z<-Enorm(as.vector(x))^2/dim(a3)[3]
#z
#}

dif<-function (a3) 
{
#version that does not depend on scale of original measurements
# assumes already centred
    cc<-centroid.size(add(a3)/dim(a3)[3])
    x <- sweep(a3, c(1, 2), apply(a3, c(1, 2), mean))
    z <- Enorm(as.vector(x)/cc)^2/dim(a3)[3]
    z
}

 fJ<-function(n)
{
	zz <- matrix(1:1, n, n)
	H <- diag(n) - (1/n) * zz
	H
}

 fcel<-function(n, d)
{
	v <- ceiling(sqrt(n))
	p <- matrix(c(0:0), n, 2)
	for(i in 1:v) {
		for(j in 1:v) {
			if((v * (i - 1) + j) < (n + 1)) {
				p[(v * (i - 1) + j), 1] <- (d/4) * (-1)^i + (d * 
				  j)
				p[(v * (i - 1) + j), 2] <- i * ((d * sqrt(3))/2
				  )
			}
		}
	}
	p
}


 fcnt<-function(a)
{
	aa <- fJ(dim(a)[1]) %*% a
	aa
}


fgpa.singleiteration<-function(a3, p)
{
# Note this is an approximation to GPA - 
# It carries out an initial match by optimally rotating all the data, 
# the rescaling the observations, then rotating the observations
# NB it does not repeat this until convergence, but in practice
# for many real datasets this gives an excellent registration
#
	zd <- list(rot. = 0, r.s.r. = 0, Gpa = 0, I.no. = 0, mshape = 0)
	zd$rot. <- rgpa(a3, p)
	zz <- rgpa(sgpa(zd$rot.$rotated), p)
	zd$r.s.r. <- zz$rotated
	zd$Gpa <- zz$dif
	zd$I.no. <- zz$r.no.
	zd$mshape <- msh(zd$r.s.r.)
	return(zd)
}

fgpa<-function(a3, tol1, tol2,proc.output=FALSE,reflect=FALSE)
{
#
#  Fully iterative fgpa (now assumes a3 is already centred)
#
#
	zd <- list(rot. = 0, r.s.r. = 0, Gpa = 0, I.no. = 0, mshape = 0)
      p<-tol1
if (proc.output){cat(" Step             | Objective function | change \n")}       
if (proc.output){cat("---------------------------------------------------\n")}       
      x1<-dif(a3)
if (proc.output){cat("Initial objective fn",x1,"  -  \n")}       
if (proc.output){cat("-----------------------------------------\n")}       
	zz <- rgpa(a3, p,proc.output=proc.output,reflect=reflect)
      x2<-dif(zz$rotated)
if (proc.output){cat("Rotation step      0",x2,x1-x2," \n")}       
if (proc.output){cat("-----------------------------------------\n")}       
ii<-1
	zz <- rgpa(sgpa(zz$rotated,proc.output=proc.output), 
	p,proc.output=proc.output,reflect=reflect)
x1<-x2
x2<-dif(zz$rotated)
rho<- x1-x2
if (proc.output){cat("Scale/rotate step ",ii,x2,rho," \n")}
if (proc.output){cat("-----------------------------------------\n")}       

if (rho > tol2){
        while (rho > tol2){
x1<-x2
ii<-ii+1
	zz <- rgpa(sgpa(zz$rotated,proc.output=proc.output), 
	p,proc.output=proc.output,reflect=reflect)
     x2<-dif(zz$rotated)
rho<- x1-x2 
if (proc.output){cat("Scale/rotate step ",ii,x2,rho," \n")}
if (proc.output){cat("-----------------------------------------\n")}       

        }
        }
	zd$r.s.r. <- zz$rotated
	zd$Gpa <- zz$dif
	zd$I.no. <- ii
	zd$mshape <- msh(zd$r.s.r.)
	return(zd)
}

fgpa.rot<-function(a3, tol1, tol2,proc.output=FALSE,reflect=FALSE)
{
# Assumes that a3 has been centred already
	zd <- list(rot. = 0, r.s.r. = 0, Gpa = 0, I.no. = 0, mshape = 0)
      p<-tol1
	zz <- rgpa(a3, p,proc.output=proc.output,reflect=reflect)
      x1<-msh(zz$rotated)       
ii<-zz$r.no.
#	zz <- rgpa(zz$rotated, p,proc.output=proc.output,reflect=reflect)
#x2<-msh(zz$rotated)
#rho<-riemdist(x1,x2)
#        while (rho > tol2){
#print(rho)
#x1<-x2
#ii<-ii+1
#	zz <- rgpa(zz$rotated, p,proc.output=proc.output)
#     x2<-msh(zz$rotated)
#rho<-riemdist(x1,x2) 
#        }
	zd$r.s.r. <- zz$rotated
	zd$Gpa <- zz$dif
	zd$I.no. <- ii
	zd$mshape <- msh(zd$r.s.r.)
	return(zd)
}


 fopa<-function(a, b)
{
abar<-fcnt(a)
bbar<-fcnt(b)
	q1 <- sum(diag( abar %*% t( abar )))
	q2 <- fos(a, b)^2 * sum(diag( bbar %*% t( bbar )))
	q3 <- 2 * fos(a, b) * sum(diag(fort(a, b) %*% t( abar ) %*% bbar ))
	gs <- q1 + q2 - q3
	gs
}


fort.ROTATEANDREFLECT<-function(a, b)
{
	x <- t(fcnt(a)) %*% fcnt(b)
        xsvd<-svd(x)
	t <- xsvd$v %*% t(xsvd$u)
	return(t)
}

fos.REFLECT<-function(a, b)
{
abar<-fcnt(a)
bbar<-fcnt(b)
	z <- ftrsq(abar, bbar )/sum(diag(t( bbar ) %*% bbar ))
	z
}


fos <- function (a, b)
{
    z <- cos(riemdist(a,b))*centroid.size(a)/centroid.size(b)
    z
}


ftrsq<-function(a, b)
{z <- sum(sqrt(abs(eigen(t(b) %*% a %*% t(a) %*% b)$values)))
	z
}

 graf<-function(a3)
{
	l <- 0
	xmin <- 0
	xmax <- 0
	ymin <- 0
	ymax <- 0
	for(i in 1:dim(a3)[3]) {
		xmin[i] <- min(a3[, 1, i])
		xmax[i] <- max(a3[, 1, i])
		ymin[i] <- min(a3[, 2, i])
		ymax[i] <- max(a3[, 2, i])
	}
	l <- c(min(xmin), min(ymin), max(xmax), max(ymax))
	plot((min(l) - 1):(max(l) + 1), (min(l) - 1):(max(l) + 1), type = "n")
	for(i in 1:dim(a3)[3]) {
		lines(close1(a3[,  , i]))
	}
}

msh<-function(a3)
{
	s <- 0
#	print("finding mean shape")
	m<-apply(a3,c(1,2),mean)
#	print("found mean shape")	
#	for(i in 1:dim(a3)[3]) {
#		s <- s + a3[,  , i]
#	}
#	m <- (1/dim(a3)[3]) * s
	return(m)
}


#Enorm<-function(a)
#{
#	return(sqrt(sum(diag(t(a) %*% a))))
#}


rgpa<-function(a3, p,reflect=FALSE,proc.output=FALSE)
{
# assumes a3 already centred now
if (reflect==TRUE) 
{ fort<-fort.ROTATEANDREFLECT } 
	zd <- list(rotated = 0, dif = 0, r.no. = 0, inc = 0)
	l <- dim(a3)[3]
	a <- 0
	d <- 0
	n <- 0
#	zz <- cnt3(a3)
	zz <- a3
#        print("Rotations ...")
#       print("Iteration,meanSS before,meanSS after,difference,tolerance")
	d[1] <- 10^12
	d[2] <- dif(zz)
	a[1] <- d[2]
	s <- add(zz)
#       print(c(d[1],d[2]))
if (dif(zz) > p){
	while(d[1] - d[2] > p) {
		n <- n + 1
		d[1] <- d[2]
		for(i in 1:l) {
			old <- zz[,  , i]
			zz[,  , i] <- old %*% fort(((1/(l - 1)) * (s - old)), 
				old)
			s <- s - old + zz[,  , i]
		}
		d[2] <- dif(zz)
		a[n + 1] <- d[2]
#		print(c(n,d[1],d[2],d[1]-d[2],p))
if (proc.output){cat("  Rotation iteration  ",n,d[2],d[1]-d[2]," \n")}
	}
	}
	zd$rotated <- zz
	zd$dif <- a
	zd$r.no. <- n
        zd$inc<-d[1]-d[2]
if (proc.output){cat("-----------------------------------------\n")}       

fort<-fort.ROTATION
	return(zd)
}


# sgpa<-function(a3)
#{
#	zz <- a3
#	a <- bgpa(zz)
#	for(i in 1:dim(a3)[3]) {
#		zz[,  , i] <- a[i] * a3[,  , i]
#	}
#	return(zz)
#}

sgpa<-function(a3,proc.output=FALSE)
{
#assumes a3 is centred
	zz <- a3
	di<-dim(a3)
	a <- bgpa(zz,proc.output=proc.output)
	i<-rep(dim(a3)[1]*dim(a3)[2],dim(a3)[3])
	sequen<-rep(a,i)
	
	zz<-array(as.vector(a3)*sequen,di)
if (proc.output){cat("  Scaling updated \n")}

	return(zz)
}



sh<-function(a)
{
	u1 <- (a[2, 1] - a[1, 1])/sqrt(2)
	u2 <- (a[2, 2] - a[1, 2])/sqrt(2)
	v1 <- (2 * a[3, 1] - a[2, 1] - a[1, 1])/sqrt(6)
	v2 <- (2 * a[3, 2] - a[2, 2] - a[1, 2])/sqrt(6)
	d <- c(0, 0)
	d[1] <- (u1 * v1 + u2 * v2)/(u1^2 + u2^2)
	d[2] <- (u1 * v2 - u2 * v1)/(u1^2 + u2^2)
	d
}

sim1<-function(n, d, s)
{
	a <- fcel(n, d)
	sig <- matrix(c(1:1), n, 1)[, 1]
	sig <- sig * s
	b <- a
	b[, 1] <- rnorm(n, mean = a[, 1], sd = sig)
	b[, 2] <- rnorm(n, mean = a[, 2], sd = sig)
	b
}
vec1<-function(a3)
{
	#zz <- array(c(0:0), dim = c((dim(a3)[1] * dim(a3)[2]), dim(a3)[3]))
	#for(i in 1:dim(a3)[3]) {
		#for(j in 1:dim(a3)[2]) {
			#for(k in 1:dim(a3)[1]) {
				#zz[((j - 1) * dim(a3)[1] + k), i] <- a3[k, j, i
				  #]
			#}
		#}
	#}
	
	zz<-matrix(a3,dim(a3)[1]*dim(a3)[2], dim(a3)[3])
	return(zz)
}

fort.ROTATION<-function(a, b)
{
	x <- t(fcnt(a)) %*% fcnt(b)
        xsvd<- svd(x)
        v<- xsvd$v
        u<- xsvd$u
	tt <- v %*% t(u)
        chk1<-Re(prod(eigen(v)$values))
        chk2<-Re(prod(eigen(u)$values))
	if ( (chk1 < 0) && (chk2 > 0) )
        {
        v[,dim(v)[2]]<-v[,dim(v)[2]]*(-1)
	tt <- v %*% t(u)
}
	if ( (chk2 < 0) && (chk1 > 0) )
        {
        u[,dim(u)[2]]<-u[,dim(u)[2]]*(-1)
	tt <- v %*% t(u)
}
 	return(tt)
}

############end of Mohammad Faghihi's (adapted) routines


#alias functions (all lower-case)

hotelling2d<-Hotelling2D
hotellingtest<-Hotellingtest
procrustesgpa<-procrustesGPA
goodall2d<-Goodall2D
goodalltest<-Goodalltest

# alias
TPSgrid<-tpsgrid

#if you wish the default to *not* include reflection 
#invariance (as is normal in shape analysis) then you need the line below. 

fort<-fort.ROTATION

#######
#exact Gaussian MLE - isotropic distribution

#######not fully tested yet

isomle<-function(x){
	if(is.complex(x)) {
		tem <- array(0, c(nrow(x), 2, ncol(x)))
		tem[, 1,  ] <- Re(x)
		tem[, 2,  ] <- Im(x)
		x <- tem
	}
	k<-dim(x)[1]
	m <- dim(x)[2]
	n<-dim(x)[3]
if (m > 2){
print("Only valid for 2D data")
}
if (m ==2){
pm<-rep(0,times=2*k-3)
tem<-procrustes2d(x)
tem1<-bookstein.shpv(tem$mshape)
sigm<-sum(diag(var(tem$tan)))/(n-1)/2
#cat("Isotropic shape MLE \n")
pm[1:(k-2)]<-tem1[3:k,1]
pm[(k-1):(2*k-4)]<-tem1[3:k,2]
pm[2*k-3]<-10
ans<-nlm(objfuniso,hessian=TRUE,pm,uu=x)
#while (ans$code!=1){
#print("code not equal 1")
#print(pm)
#pm<-pm+rnorm(2*k-3,0,0.1)
#pm[2*k-3]<-abs(pm[2*k-3])
#ans<-nlm(objfuniso,hessian=TRUE,pm,uu=x) #print(ans)
#}
out<-list(code=0,mshape=0,tau=0,kappa=0,varcov=0,gradient=0)
mn<-matrix(0,k,2)
mn[1,1]<--0.5
mn[2,1]<-0.5
mn[3:k,1]<-ans$estimate[1:(k-2)]
mn[3:k,2]<-ans$estimate[(k-1):(2*k-4)]
out$mshape<-mn
out$code<-ans$code
out$loglike <- -ans$minimum
out$gradient<-ans$gradient
out$tau<-sqrt(1/ans$estimate[2*k-3]**2)
out$kappa<-centroid.size(mn)**2/(4*out$tau**2)
out$varcov<-solve(ans$hessian)
out$se<-c(sqrt(diag(out$varcov)))
out$se[2*k-3] <- out$se[2*k-3] *out$tau**2
out
}
}
objfuniso<-function(pm,uu){ k<-dim(uu)[1]
h<-defh(k-1)
zero<-matrix(0,k-1,k)
L1<-cbind(h,zero)
L2<-cbind(zero,h)
L<-rbind(L1,L2)
mustar<-c(-1/2,1/2,pm[1:(k-2)],0,0,pm[(k-1):(2*k-4)])
mu<-L%*%mustar
obj<--loglikeiso2(uu,mu,1/pm[2*k-3])
obj
}

loglikeiso<-function(uu,mu,s){
nsam<-dim(uu)[3]
sum<-0
for (i in 1:nsam){
sum<-sum+log(isodens(uu[,,i],mu,s))
}
sum
}

loglikeiso2<-function(uu,mu,s){
nsam<-dim(uu)[3]
sum<-0
for (i in 1:nsam){
sum<-sum+isologdens(uu[,,i],mu,s)
}
sum
}

isodens<-function(usam,mu,s){
k<-dim(usam)[1]
u<-kendall.shpv(usam)
uuu<-u[,1]
vvv<-u[,2]
up<-c(1,uuu,0,vvv)
vp<-c(0,-vvv,1,uuu)
usu<-t(up)%*%up
beta<-c(t(mu)%*%up,t(mu)%*%vp)
sin2rho<-1-t(beta)%*%beta/(usu*c(t(mu)%*%mu))
kappa<-c(t(mu)%*%mu)/(4*s**2)
#finf<-gamma(k-1)*pi/(pi*usu)**(k-1)
dens<-oneFone(k-2,2*kappa*(1-sin2rho))%*%exp(-2*kappa*sin2rho)
dens
}

isologdens<-function(usam,mu,s){
k<-dim(usam)[1]
u<-kendall.shpv(usam)
uuu<-u[,1]
vvv<-u[,2]
up<-c(1,uuu,0,vvv)
vp<-c(0,-vvv,1,uuu)
usu<-t(up)%*%up
beta<-c(t(mu)%*%up,t(mu)%*%vp)
sin2rho<-1-t(beta)%*%beta/(usu*c(t(mu)%*%mu))
kappa<-c(t(mu)%*%mu)/(4*s**2)
#finf<-lgamma(k-1)+log(pi)-(k-1)*log(pi*usu)
dens<-loneFone(k-2,2*kappa*(1-sin2rho))-2*kappa*sin2rho
c(dens)
}
loneFone<-function(r,x){
#note this is log 1F1(-r,1,-x)
if (x > 1){
sum1<-r*log(x)
sum<-0
for (j in 0:r){
sum<-sum+choose(r,j)*x**(j-r)/gamma(j+1)
}
out<-sum1+log(sum)
}
if (x <=1){
sum<-0
for (j in 0:r){
sum<-sum+choose(r,j)*x**(j)/gamma(j+1)
}
out<-log(sum)
}
out
}


kendall.shpv<-function(x){
k<-dim(x)[1]
h<-defh(k-1)
zz<-h%*%x
kendall<-(zz[2:(k-1),1]+1i*zz[2:(k-1),2])/(zz[1,1]+1i*zz[1,2])
kendall<-cbind(Re(kendall),Im(kendall))
kendall
}

oneFone<-function(r,x){
#note this is 1F1(-r,1,-x)
sum<-0
for (j in 0:r){
sum<-sum+choose(r,j)*x**j/gamma(j+1)
}
sum
}

schizophrenia.dat<-c( 0.345632   , -0.0360314
    , -0.356301   , 0.0234333
   , 0.0119311     , 0.17692
     , 0.37789    , 0.480402
    , 0.719631     , -0.41189
    , 0.397921    , -0.140558
    , 0.351751    , -0.385748
    , 0.333756    , -0.655051
    , 0.032181    , -0.275235
    , 0.112563    , -0.506533
    , -0.233126     , -0.28334
    , -0.667337   , 0.0522613
    , 0.188945    , -0.142714
    , 0.237198    , 0.048306
    , -0.340236   , 0.0997385
   , -0.0161814    , 0.201017
    , 0.301584    , 0.516546
    , 0.510795    , -0.323537
    , 0.269407   , -0.0562209
    , 0.239301    , -0.253218
    , 0.229339     , -0.48236
   , 0.0201328    , -0.122625
    , 0.048306    , -0.341874
    , -0.200997    , -0.122698
     , -0.62316    , 0.120534
    , 0.124688   , -0.0262477
    , 0.341616    , 0.048306
    , -0.408509    , 0.119819
  , -0.00814926    , 0.277322
     , 0.39797     , 0.62498
    , 0.591116    , -0.299441
     , 0.35776   , -0.0361405
     , 0.27143    , -0.249202
    , 0.273516    , -0.530553
    , 0.032181    , -0.110577
   , 0.0724024     , -0.34589
    , -0.188949    , -0.138762
    , -0.707498    , 0.140615
    , 0.124688   , -0.0262477
    , 0.329567     , -0.10832
    , -0.359859   , 0.0341931
    , -0.056342    , 0.152824
    , 0.377668    , 0.474464
    , 0.673363    , -0.462009
    , 0.373825    , -0.228912
    , 0.323639    , -0.450005
     , 0.34179    , -0.735372
   , 0.0924219    , -0.295315
    , 0.100555    , -0.540522
    , -0.181195    , -0.260518
    , -0.651361   , 0.0963372
    , 0.174798    , -0.178741
    , 0.193021   , 0.0683864
    , -0.539516    , 0.210918
     , -0.16074    , 0.333555
    , 0.293246    , 0.639288
    , 0.520753    , -0.301366
    , 0.245311   , -0.0281084
    , 0.142916     , -0.28133
    , 0.128938    , -0.510473
    , -0.152558    , -0.118609
   , -0.0761516    , -0.379879
    , -0.414127    , -0.127988
    , -0.848201    , 0.246974
  , -0.00592517  , -0.00203414
    , 0.337599    , -0.076192
    , -0.356431    , 0.146356
    , 0.068156    , 0.201017
    , 0.520571    , 0.464342
    , 0.601074    , -0.518234
    , 0.321616    , -0.208831
    , 0.251349    , -0.409844
    , 0.197211    , -0.695211
  , -0.00797966    , -0.251139
  , -0.00787855    , -0.500361
     , -0.26726    , -0.201009
    , -0.719602     , 0.25526
    , 0.110557    , -0.140611
    , 0.223261    , 0.228767
    , -0.431293    , 0.230133
   , -0.0439952    , 0.389752
     , 0.23941    , 0.819652
    , 0.641234    , -0.088515
    , 0.377841    , 0.140566
    , 0.339703     , -0.10864
    , 0.393998    , -0.345813
   , 0.0201328   , -0.0583678
    , 0.172844    , -0.283494
    , -0.268058   , -0.0890055
    , -0.777652    , 0.251188
    , 0.154806     , 0.12832
    , 0.243341    , -0.052358
    , -0.451374    , 0.145795
   , -0.0480112    , 0.297382
    , 0.408185    , 0.574616
     , 0.56291    , -0.538339
     , 0.31561    , -0.190749
    , 0.251349    , -0.437957
    , 0.169098    , -0.715291
  , 0.00005         , -0.194914
 , 0.000153631    , -0.492329
    , -0.278085    , -0.117118
    , -0.797732    , 0.255204
    , 0.122674   , -0.0905492
    , 0.287518   , 0.0600918
    , -0.503583    , 0.121699
   , -0.0861578    , 0.307424
    , 0.355977    , 0.614776
    , 0.603071    , -0.277295
    , 0.303562   , -0.0100258
    , 0.259382    , -0.293379
    , 0.241387    , -0.538584
   , -0.0521564    , -0.114593
   , 0.0443303    , -0.343735
    , -0.274068   , -0.0930217
    , -0.805764    , 0.138738
     , 0.12669  , 0.00182023
    , 0.319646     , 0.20467
    , -0.431294    , 0.153827
   , -0.0660775    , 0.391762
     , 0.29172    , 0.793493
    , 0.771746  , -0.00420256
    , 0.407979    , 0.166681
    , 0.387896     , -0.12872
    , 0.393998    , -0.394006
   , 0.0864079   , -0.0382915
    , 0.168828    , -0.279477
    , -0.225876    , -0.109086
    , -0.779633    , 0.126677
    , 0.190947    , 0.126318
    , 0.303582    , 0.208686
    , -0.395149    , 0.302422
  , -0.00985258    , 0.387746
    , 0.388105     , 0.72522
    , 0.611103    , -0.289343
    , 0.327658   , 0.0743116
    , 0.283478    , -0.229121
    , 0.201226    , -0.494408
  , 0.00608663   , -0.0141951
   , 0.0362983    , -0.279477
     , -0.24194  , -0.00466831
    , -0.699312    , 0.339529
     , 0.12669   , 0.0901736
    , 0.287518   , -0.0362937
    , -0.428786    , 0.133507
   , -0.0801519    , 0.233129
    , 0.374063    , 0.614792
    , 0.671344    , -0.478098
    , 0.383883    , -0.134524
    , 0.315606    , -0.417876
    , 0.265483    , -0.675131
   , 0.0342193     , -0.19491
   , 0.0443223    , -0.466226
    , -0.266036    , -0.201455
    , -0.735453    , 0.202987
    , 0.176577   , -0.0906634
    , 0.251373    , 0.212702
     , -0.42477    , 0.149571
    , -0.124329    , 0.345579
    , 0.145148    , 0.767403
    , 0.671344  , -0.00821852
    , 0.424044    , 0.138569
    , 0.407976    , -0.116671
     , 0.44219    , -0.377942
    , 0.122573   , -0.0583636
    , 0.202939    , -0.283647
     , -0.14957   , -0.0970378
    , -0.755533    , 0.114633
    , 0.180593   , 0.0940755
     , 0.29555   , -0.0362937
    , -0.348465 , 0.000976592
   , -0.0480233    , 0.148792
    , 0.273662    , 0.458166
     , 0.65528    , -0.341552
    , 0.371835    , -0.106411
    , 0.327655    , -0.377715
    , 0.333757     , -0.65505
   , 0.0984763    , -0.247119
    , 0.118601    , -0.500514
    , -0.181699    , -0.253664
    , -0.675212   , 0.0182479
    , 0.192642    , -0.122792
    , 0.309487    , 0.265173
      , -0.3438    , 0.289185
   , 0.0199632    , 0.482141
    , 0.355444    , 0.735605
    , 0.613122    , -0.156788
    , 0.357761    , 0.128518
    , 0.315607    , -0.132736
    , 0.305644     , -0.36991
   , 0.0121007   , 0.0380178
    , 0.100555    , -0.223253
    , -0.182924   , 0.0134128
    , -0.623245    , 0.301185
    , 0.162814    , 0.142714
    , 0.193021    , 0.305334
    , -0.472766    , 0.163996
    , -0.172808     , 0.38174
   , 0.0847168    , 0.829799
    , 0.538908    , 0.114214
    , 0.249327    , 0.273096
    , 0.199141   , 0.0319229
    , 0.273516    , -0.165091
   , -0.0360921   , 0.0340017
   , 0.0563381    , -0.141071
    , -0.257222   , -0.0263121
    , -0.791835   , 0.0884059
  , -0.00382644    , 0.178572
     , 0.39784    , -0.011935
    , -0.440185    , 0.132558
  , -0.00814926    , 0.309451
    , 0.470037    , 0.542738
    , 0.580994    , -0.385704
    , 0.389889    , -0.212847
    , 0.283478    , -0.401812
    , 0.237372    , -0.675131
   , 0.0201328    , -0.267203
   , 0.0644106    , -0.508393
     , -0.23742    , -0.196261
    , -0.711598    , 0.184719
   , 0.0944764    , -0.130548
    , 0.269326   , 0.0844506
    , -0.327735    , 0.048221
   , -0.0242135    , 0.253226
    , 0.257186    , 0.643139
    , 0.601074      , -0.1849
    , 0.325632   , 0.0160683
    , 0.307574    , -0.217073
    , 0.301629    , -0.442199
   , 0.0803737    , -0.126641
    , 0.112603    , -0.351767
    , -0.153083     , -0.13602
    , -0.635293   , 0.0280923
    , 0.158733  , -0.00605034
    , 0.217118    , 0.212965
    , -0.416089    , 0.192799
    , -0.104535    , 0.353627
    , 0.229073    , 0.727476
    , 0.601074      , -0.1849
    , 0.293503    , 0.112454
    , 0.267414    , -0.112655
    , 0.285565    , -0.402039
   , 0.0281649   , -0.0101749
   , 0.0844909    , -0.259397
    , -0.221356   , -0.0275863
    , -0.743727    , 0.188735
   , 0.0944764    , 0.118448
    , 0.293423    , -0.076192
    , -0.327736   , 0.0442049
   , -0.0121653    , 0.160856
    , 0.357587    , 0.394143
    , 0.572961    , -0.457993
    , 0.281455    , -0.228912
    , 0.263398    , -0.470085
    , 0.189179    , -0.679147
  , 0.00005         , -0.259171
     , 0.02425    , -0.524458
    , -0.193244    , -0.208309
    , -0.583084   , 0.0883332
    , 0.142669    , -0.178741
    , 0.289407   , 0.0201935
    , -0.351832    , 0.164687
  , 0.00389893    , 0.265274
     , 0.40578    , 0.510609
     , 0.49264    , -0.433896
      , 0.3176    , -0.120478
    , 0.235285    , -0.329523
    , 0.201227    , -0.554649
  , 0.00406852     , -0.17885
   , 0.0081857    , -0.403976
    , -0.245452    , -0.107908
    , -0.647341    , 0.228896
   , 0.0583318   , -0.0622752
    , 0.301455   , 0.0603542
    , -0.391993    , 0.224928
   , -0.0121653    , 0.341579
    , 0.441925    , 0.594946
     , 0.45248    , -0.389719
     , 0.29752   , -0.0963815
     , 0.17906     , -0.30141
    , 0.140986    , -0.570714
     , -0.02806    , -0.106561
   , -0.0560713    , -0.391927
    , -0.245452   , -0.0556988
    , -0.683486    , 0.293153
   , 0.0583318  , -0.00605033
    , 0.237198   , -0.0922563
    , -0.396005  , 0.00807245
    , -0.116583     , 0.11668
     , 0.29333    , 0.402175
     , 0.47256    , -0.502169
     , 0.22523     , -0.22088
     , 0.17906    , -0.433941
    , 0.157051    , -0.634971
   , -0.0320761    , -0.271219
   , -0.0400071    , -0.500361
    , -0.257501    , -0.260518
     , -0.69955   , 0.0521888
   , 0.0382515    , -0.158661
    , 0.321535   , -0.0641438
    , -0.379941    , 0.108474
   , -0.0201975    , 0.192985
    , 0.445941    , 0.450368
     , 0.49264    , -0.506185
    , 0.325632    , -0.184735
    , 0.203157    , -0.365667
    , 0.124922    , -0.675131
   , -0.0119958    , -0.206962
   , -0.0279589    , -0.476265
    , -0.285613    , -0.192245
    , -0.651357    , 0.176687
   , 0.0944764    , -0.126532
    , 0.317519   , -0.0761919
    , -0.219297   , -0.0300639
    , 0.072172    , 0.140776
     , 0.40578    , 0.430287
    , 0.580994    , -0.457993
    , 0.345712    , -0.216863
    , 0.287494    , -0.417876
      , 0.2695    , -0.663083
    , 0.096438    , -0.259171
   , 0.0764588    , -0.496345
    , -0.165131     , -0.24847
    , -0.518831  , -0.00808046
    , 0.162749    , -0.178741
    , 0.363819    , -0.132663
    , -0.295482    , 0.140679
    , 0.178874    , 0.174894
    , 0.568759     , 0.29153
    , 0.593042    , -0.550362
    , 0.345712    , -0.297185
    , 0.267414     , -0.50623
    , 0.136969     , -0.74742
   , -0.0039636    , -0.275235
   , 0.0724427    , -0.520442
    , -0.239946    , -0.153262
    , -0.524623    , 0.229146
    , 0.194967    , -0.196981
    , 0.255389    , 0.160493
    , -0.435309    , 0.214068
    , -0.124316    , 0.357623
    , 0.267522     , 0.65901
    , 0.597058    , -0.136708
    , 0.305552   , 0.0763094
    , 0.255365    , -0.140768
    , 0.257451    , -0.402038
    , 0.032181   , -0.0382875
    , 0.084491    , -0.279478
    , -0.243962   , -0.0167163
    , -0.761587    , 0.231107
   , 0.0945651   , 0.0881595
    , 0.279486  , 0.00788297
    , -0.330892    , 0.141779
   , 0.0222759    , 0.267264
    , 0.386119    , 0.484248
    , 0.574959    , -0.365648
     , 0.33569    , -0.102395
    , 0.263398    , -0.345587
    , 0.233355    , -0.566697
   , -0.0180098    , -0.178854
   , 0.0965392    , -0.387911
    , -0.233908     , -0.12515
    , -0.627023    , 0.206999
   , 0.0681437   , -0.0866473)
schizophrenia.dat<-array(schizophrenia.dat,c(2,13,28))
schizophrenia.dat<-aperm(schizophrenia.dat,c(2,1,3))


permutationtest<-function(A,B,nperms=200){

A1<-A
A2<-B
B <- nperms
nsam1<-dim(A1)[3]
nsam2<-dim(A2)[3]



Gtem<-Goodalltest(A1,A2)
Htem<-Hotellingtest(A1,A2)


Gumc<-Gtem$F
Humc<-Htem$F

Gtabpval<-Gtem$pval
Htabpval<-Htem$pval


if (B > 0){

Apool <- array( 0, c(dim(A1)[1],dim(A1)[2],dim(A1)[3]+dim(A2)[3] ) )
Apool[,,1:nsam1]<-A1
Apool[,,(nsam1+1):(nsam1+nsam2)] <- A2

out<-list(H=0,H.pvalue=0,H.table.pvalue=0,G=0,G.pvalue=0,G.table.pvalue=0)

Gu<-rep(0,times=B)
Hu<-rep(0,times=B)



cat("Permutations - sampling without replacement: ")
cat(c("No of permutations = ",B,"\n"))

for (i in 1:B){
cat(c(i," "))


select<-sample( 1: (nsam1+nsam2) )


Gu[i]<-Goodalltest( Apool[,,select[1:nsam1] ] , Apool[,,select[ (nsam1+1):(nsam2+nsam1) ] ])$F
Hu[i]<-Hotellingtest(Apool[,,select[1:nsam1] ],Apool[,,select[ (nsam1+1):(nsam1+nsam2) ] ] )$F

}

Gu<-sort(Gu)
numbig<-length(Gu[Gumc<Gu])
pvalG<-(1+numbig)/(B+1)
Hu<-sort(Hu)
numbig<-length(Hu[Humc<Hu])
pvalH<-(1+numbig)/(B+1)
cat(" \n")

out$H<-Humc
out$H.pvalue<-pvalH
out$H.table.pvalue<-Htabpval

out$G<-Gumc
out$G.pvalue<-pvalG
out$G.table.pvalue<-Gtabpval

}

if (B == 0){

out<-list(H=0,H.table.pvalue=0,G=0,G.table.pvalue=0)

out$H<-Humc
out$H.table.pvalue<-Htabpval

out$G<-Gumc
out$G.table.pvalue<-Gtabpval
}
out
}

permutationtest <- permutationtest2


frechet<-function(x,mean="intrinsic"){
if (mean=="intrinsic"){
option <- 1
}
if (mean=="partial.procrustes"){
option <- 2
}
if (mean=="full.procrustes"){
option <- 3
}
if (mean=="mle"){
option <- 4
}
if (is.double(mean)){
if (mean>0){
option <- -mean
}
}
n<-dim(x)[3]
for (i in 1:n){
x[,,i]<-x[,,i]/centroid.size(x[,,i])
}
if (option < 4){
pm<-procGPA(x,scale=FALSE,tol1=10^(-8))$mshape
m<-dim(x)[2]
k<-dim(x)[1]
ans<-list(mshape=0,var=0,code=0,gradient=0)
out <- nlm(objfun, hessian=TRUE, c(pm), uu=x,option=option,iterlim=1000)
B<-matrix(out$estimate,k,m)
ans$mshape<-procOPA(pm,B)$Bhat
ans$var<-out$minimum
ans$code<-out$code
ans$gradient<-out$gradient
}
if (option==4){
pm<-procGPA(x,scale=FALSE,tol1=10^(-8))$mshape
m<-dim(x)[2]
k<-dim(x)[1]
if (m==2){
theta<-c(log(centroid.size(pm)**2/(4*0.1**2)),pm)
ans<-list(mshape=0,kappa=0,code=0,gradient=0)
out <- nlm(objfun4, hessian=TRUE, theta, uu=x,iterlim=1000)
B<-matrix(out$estimate[-1],k,m)
ans$mshape<-procOPA(pm,B)$Bhat
ans$kappa<-exp(out$estimate[1])
ans$loglike <- -out$minimum
ans$code<-out$code
ans$gradient<-out$gradient
}
if (m!=2){
print("MLE is only appropriate for planar shapes")
}
}
ans
}




objfun<-function(pm,uu,option){
m<-dim(uu)[2]
k<-dim(uu)[1]
pm<-matrix(pm,k,m)
sum<-0
for (i in 1:dim(uu)[3]){

if (option==1){
sum<-sum+(riemdist(pm,uu[,,i]))**2
}
if (option==2){
sum<-sum+4*sin(riemdist(pm,uu[,,i])/2)**2
}
if (option==3){
sum<-sum+sin(riemdist(pm,uu[,,i]))**2
}
if (option<0){
h <- -option
sum<-sum+  ((1 - cos(riemdist(pm,uu[,,i]))**(2*h) )/h)
}
}
sum
}


objfun4<-function(pm,uu){
m<-dim(uu)[2]
k<-dim(uu)[1]
n<-dim(uu)[3]
kappa<-exp(pm[1])
pm<-matrix(pm[-1],k,m)
sum<-0
for (i in 1:n){
sin2rho<- sin(riemdist(pm,uu[,,i]))**2
sum<-sum+loneFone(k-2,2*kappa*(1-sin2rho))-2*kappa*sin2rho
}
-sum
}

MDSshape<-function(x,alpha=1,projalpha=1/2){
mu<-procGPA(x)$mshape
k<-dim(x)[1]
n<-dim(x)[3]
m<-dim(x)[2]
H<-defh(k-1)
sum<-matrix(0,k-1,k-1)
for (i in 1:n){
Z<-preshape(x[,,i])
if (alpha==1){
sum<-sum+(Z)%*%t((Z))
}
if (alpha==1/2){
ee<-eigen((Z)%*%t((Z)),symmetric=TRUE)
sum<-sum+ee$vectors%*%diag(sqrt(abs(ee$values)))%*%t(ee$vectors)
}
}
eig<-eigen(sum/n,symmetric=TRUE)
lam<-eig$values

if (m == 2){
if (projalpha==1/2){
meanshape<-cbind( t(H)%*%(sqrt(lam[1])*eig$vectors[,1])/sqrt(lam[1]+lam[2]) , -t(H)%*%(sqrt(lam[2])*eig$vectors[,2])/sqrt(lam[1]+lam[2]) )
}
if (projalpha==1){
lambar<-(lam[1]+lam[2])/2
meanshape<-cbind( t(H)%*%(sqrt(lam[1]-lambar+1/m)*eig$vectors[,1]) , -t(H)%*%(sqrt(lam[2]-lambar+1/m)*eig$vectors[,2]) )
}
}

if (m == 3){
if (projalpha==1/2){
meanshape<-cbind( t(H)%*%(sqrt(lam[1])*eig$vectors[,1])/sqrt(lam[1]+lam[2]+lam[3]) , t(H)%*%(sqrt(lam[2])*eig$vectors[,2])/sqrt(lam[1]+lam[2]+lam[3]),t(H)%*%(sqrt(lam[3])*eig$vectors[,3])/sqrt(lam[1]+lam[2]+lam[3]) )
}
if (projalpha==1){
lambar<-(lam[1]+lam[2]+lam[3])/3
meanshape<-cbind( t(H)%*%(sqrt(abs(lam[1]-lambar+1/m))*eig$vectors[,1]) , t(H)%*%(sqrt(abs(lam[2]-lambar+1/m))*eig$vectors[,2]) ,t(H)%*%(sqrt(abs(lam[3]-lambar+1/m))*eig$vectors[,3]))
}
}
if (riemdist(meanshape,mu)>riemdist(meanshape,mu,reflect=TRUE)){
meanshape[,m]<--meanshape[,m]
}
meanshape
}




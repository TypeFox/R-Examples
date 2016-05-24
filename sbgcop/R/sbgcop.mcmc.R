"sbgcop.mcmc" <-
function(Y,S0=diag(dim(Y)[2]),n0=dim(Y)[2]+2, 
 nsamp=100,odens=max(1,round(nsamp/1000)),impute=any(is.na(Y)),
 plugin.threshold=100,
 plugin.marginal=(apply(Y,2,function(x){ length(unique(x))})>plugin.threshold),
 seed=1,verb=TRUE){

########## check input
ok_S0<-all(eigen(S0)$val>0) & dim(S0)[1]==dim(Y)[2] & dim(S0)[2]==dim(Y)[2]
ok_n0<-(n0>=0)

if(!ok_S0) { cat("Error: S0 must be a positive definite p x p matrix \n") }
if(!ok_n0) { cat("Error: n0 must be positive \n") }

if(ok_S0 & ok_n0) {

########## data
vnames<-colnames(Y) 
Y<-as.matrix(Y)
colnames(Y)<-vnames
n<-dim(Y)[1]
p<-dim(Y)[2]
##########

########## starting values
set.seed(seed)
R<-NULL
for(j in 1:p) { R<-cbind(R, match(Y[,j],sort(unique(Y[,j])))) }
Rlevels<-apply(R,2,max,na.rm=TRUE)
Ranks<- apply(Y,2,rank,ties.method="max",na.last="keep")
N<-apply(!is.na(Ranks),2,sum)
U<- t( t(Ranks)/(N+1))
Z<-qnorm(U)
Zfill<-matrix(rnorm(n*p),n,p)
Z[is.na(Y)]<-Zfill[is.na(Y) ]
S<-cov(Z)
##########

########## things to keep track of
Y.pmean<-Y
if(impute){Y.pmean<-matrix(0,nrow=n,ncol=p)}
LPC<-NULL
C.psamp<-array(dim=c(p,p,floor(nsamp/odens)))
Y.imp<-NULL
if(impute){Y.imp<-array(dim=c(n,p,floor(nsamp/odens) )) }
dimnames(C.psamp)<-list(colnames(Y),colnames(Y),1:floor(nsamp/odens))
##########

########## start MCMC
for(ns in 1:nsamp) {

#### update Z
for(j in sample(1:p)) {
Sjc<- S[j,-j]%*%solve(S[-j,-j])
sdj<- sqrt( S[j,j] -S[j,-j]%*%solve(S[-j,-j])%*%S[-j,j]  )
muj<- Z[,-j]%*%t(Sjc)

if(!plugin.marginal[j])
{
  for(r in 1:Rlevels[j]){
  ir<- (1:n)[R[,j]==r & !is.na(R[,j])]
  lb<-suppressWarnings(max( Z[ R[,j]==r-1,j],na.rm=TRUE))
  ub<-suppressWarnings(min( Z[ R[,j]==r+1,j],na.rm=TRUE))
  Z[ir,j]<-qnorm(runif(length(ir),
           pnorm(lb,muj[ir],sdj),pnorm(ub,muj[ir],sdj)),muj[ir],sdj)
                       }
}
ir<-(1:n)[is.na(R[,j])]
Z[ir,j]<-rnorm(length(ir),muj[ir],sdj)
                          }
#### 

#### update S
S<-solve(rwish(solve(S0*n0+t(Z)%*%Z),n0+n))
####

#### save results
if(ns %% odens ==0) {
C<-S/(sqrt(diag(S))%*%t(sqrt(diag(S))))

lpc<-ldmvnorm(Z%*%diag(1/sqrt(diag(S))),C)
LPC<-c(LPC,lpc)
C.psamp[,,ns/odens]<-C

if(impute)
{
Y.imp.s<-Y
for(j in 1:p) {
Y.imp.s[is.na(Y[,j]),j]<-quantile(Y[,j],
                       pnorm(Z[is.na(Y[,j]),j],0,sqrt(S[j,j])),na.rm=TRUE,type=1)
               }
Y.imp[,,ns/odens]<-Y.imp.s
Y.pmean<-  ( (ns/odens-1)/(ns/odens) )*Y.pmean+  (1/(ns/odens) )*Y.imp.s 
}
                  }

if(verb==TRUE & (ns %%(odens*10))==0){
        cat(round(100*ns/nsamp),"percent done ",date(),"\n")
                                     }
####

         }
##########

G.ps<-list(C.psamp=C.psamp,Y.pmean=Y.pmean,Y.impute=Y.imp,LPC=LPC)
class(G.ps)<-"psgc"
return(G.ps)
              }
}


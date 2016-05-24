`stability.par` <-
function(data,rep,MSerror,alpha=0.1,main=NULL,cova=FALSE,
name.cov=NULL,file.cov=0,console=FALSE){
#-----------
KK <- "Environmental index"
y<-data
if (cova) {
x<- as.matrix(file.cov)
KK <- name.cov # cuando hay covariable
}

A <- nrow(y) # Number of genotypes
M <- ncol(y) # Number of environments
N <- rep  # Number of replications
MKE <- MSerror # Pooled error mean square
RR <- main # "Stability" Title of the experiments
#
FM0 <- qf(1-alpha,M-1,M * (A - 1) * (N - 1))
# Data
#dimension for genotype (A)
dimA <- rep(0,A); dim(dimA)<-A
SHV<-dimA; MV<-dimA; SU<-dimA; MV1<-dimA; FF<-dimA; GY<-dimA; U1<-dimA; G<-dimA
SI<-dimA; B<-dimA; SA<-dimA; S<-dimA; FS<-dimA; FSS<-dimA; R<-dimA; GYS<-dimA
F1<-dimA; GY<-dimA; NN<-dimA; X1<-dimA; X1M<-dimA; W<-dimA;GYY<-dimA; MMM<-dimA
#dimension for environment (M)
dimA <- rep(0,M); dim(dimA)<-M
SHM<-dimA; MM<-dimA; II<-dimA; X2<-dimA; X2M<-dimA
#dimension (A,M)
dimA <- rep(0,A*M); dim(dimA)<-c(A,M)
U<-dimA;  G1<-dimA
SV <- 0  ; SM <- 0 ; GG1 <-0  ; L <- 0;  SMES <- 0
SS<- sum(y^2)
SHT<- sum(y)
for ( i in 1: A) {
SHV[i] <- sum(y[i, ])
MV[i] <- SHV[i] / M
SV <- SV + SHV[i] ^ 2 / M
}
FK <- SHT ^ 2 / (A * M)
SKV <- (SV - FK) * N
MKV = SKV / (A - 1)
for ( j in 1: M) {
SHM[j] <- sum(y[, j])
MM[j] <- SHM[j] / A
SM <- SM + SHM[j] ^ 2 / A
II[j] <- MM[j] - SHT / (A * M)
}
SKM <- (SM - FK) * N
MKM <- SKM / (M - 1)
SKT <- (SS - FK) * N
SKVM <-SKT - SKV - SKM; MKVM <- SKVM / ((A - 1) * (M - 1))
for ( i in 1: A) {
for ( j in 1: M) {
U[i, j] <- y[i, j] - MM[j]
SU[i] <- SU[i] + U[i, j]
}
U1[i] <- SU[i] / M
}

b <- 1 / ((M - 1) * (A - 1) * (A - 2))
for ( i in 1: A) {
for ( j in 1: M) {

G[i] <- G[i] + (U[i, j] - U1[i]) ^ 2
}
GG1 <- GG1 + G[i]
}
for ( i in 1: A) {

SI[i] <- b * (A * (A - 1) * G[i] - GG1) * N
}

if ( cova ) ZZ <- sum(x^2)
if (!cova ) IN <- sum(II^2)

for ( i in 1: A) {
B[i]<-0
for ( j in 1: M) {
if ( cova ) B[i] <- B[i] + ((U[i, j] - U1[i]) * x[j] / ZZ)
if (!cova ) B[i] <- B[i] + ((U[i, j] - U1[i]) * II[j] / IN)
}
}

for ( i in 1: A) {
SA[i]<-0
for ( j in 1: M) {
if ( cova ) SA[i] <- SA[i] + (U[i, j] - U1[i] - B[i] * x[j]) ^ 2
if (!cova ) SA[i] <- SA[i] + (U[i, j] - U1[i] - B[i] * II[j]) ^ 2
}
}

for ( i in 1: A) {
L <- L + SA[i]
}
for ( i in 1: A) {
S[i] <- (A / ((A - 2) * (M - 2))) * (SA[i] - L / (A * (A - 1))) * N
}
KI <- ((A - 1) * (M - 2)) / A
SS<- 0
for ( i in 1: A) {
SS <- SS + S[i] * KI
}

SKMB <- SS
MS <- SKMB / ((A - 1) * (M - 2))
SKH <- SKVM - SKMB
MKH <- SKH / (A - 1)
SHKLT <- M * (A - 1) * (N - 1)

T05 <- qt(0.95,SHKLT)
DMV05 <- T05 * sqrt(2 * MKE / (N * M))
for ( i in 1: A) {
SMES <- SMES + MV[i]
}
MES <- SMES / A
for ( i in 1: A) {
if( MV[i] > MES) MV1[i] <- 1
if( MV[i] >= (MES + DMV05) ) MV1[i] <- 2
if( MV[i] >= (MES + 2 * DMV05)) MV1[i] <- 3
if( MV[i] < MES ) MV1[i] <- -1
if( MV[i] <= (MES - DMV05)) MV1[i] <- -2
if( MV[i] <= (MES - 2 * DMV05)) MV1[i] <- -3
}

#-------
FV <- MKV / MKVM; FM <- MKM / MKE
FVM <-MKVM / MKE; FH <- MKH / MS
FMS <- MS / MKE

for ( i in 1: A) {
FS[i] <- SI[i] / MKE
FSS[i] <- S[i] / MKE
}
SHF2 <-(A - 1) * (M - 1); SHF1 = (A - 1)
# value F and t-student
F05 <- qf(0.95,SHF1, SHF2)
F01 <- qf(0.99,SHF1, SHF2)
pvalue <- round(1-pf( FV ,SHF1, SHF2),3)
DD<-paste("",pvalue)
if (pvalue <0.001) DD<-"<0.001"
SHF2 <- M * (A - 1) * (N - 1)
SHF1 <- M - 1
F05 <- qf(0.95,SHF1, SHF2)
F01 <- qf(0.99,SHF1, SHF2)
pvalue<- round(1-pf( FM ,SHF1, SHF2),3)
NNN<-paste("",pvalue)
if (pvalue <0.001) NNN<-"<0.001"
SHF2 <- M * (A - 1) * (N - 1); SHF1 <- (A - 1) * (M - 1)
F05 <- qf(0.95,SHF1, SHF2)
F01 <- qf(0.99,SHF1, SHF2)
pvalue <-  round(1-pf( FVM ,SHF1, SHF2),3)
LL<-paste("",pvalue)
if (pvalue <0.001) LL<-"<0.001"
SHF2 <- (A - 1) * (M - 2); SHF1 <- A - 1
F05 <- qf(0.95,SHF1, SHF2)
F01 <- qf(0.99,SHF1, SHF2)
pvalue <-  round(1-pf( FH ,SHF1, SHF2),3)
HH<-paste("",pvalue)
if (pvalue <0.001) HH<-"<0.001"
SHF2 <- M * (A - 1) * (N - 1); SHF1 <- (A - 1) * (M - 2)
F05 <- qf(0.95,SHF1, SHF2)
F01 <- qf(0.99,SHF1, SHF2)
pvalue <-  round(1-pf( FMS ,SHF1, SHF2),4)
BB<-paste("",pvalue)
if (pvalue <0.001) BB<-"<0.001"
SHF2 <- M * (A - 1) * (N - 1)
SHF1 <- M - 1
F05 <- qf(0.95,SHF1, SHF2)
F01 <- qf(0.99,SHF1, SHF2)
for ( i in 1: A) {
if(FS[i] >= F01 ) NN[i] <- "**"
if( (FS[i] < F01) & (FS[i] >= F05) ) NN[i] <- "*"
if( FS[i] < F05 ) NN[i] <- "ns"
}
SHF2 <- M * (A - 1) * (N - 1)
SHF1 <- M - 2
F05 <- qf(0.95,SHF1, SHF2)
F01 <- qf(0.99,SHF1, SHF2)
for ( i in 1: A) {
if( FSS[i] >= F01 ) MMM[i] <- "**"
if( (FSS[i] < F01) & (FSS[i] >= F05) ) MMM[i] <- "*"
if( FSS[i] < F05 ) MMM[i] <- "ns"
}

#-------
if(console){
cat("\n","INTERACTIVE PROGRAM FOR CALCULATING SHUKLA'S STABILITY VARIANCE AND KANG'S")
cat("\n","                       YIELD - STABILITY (YSi) STATISTICS")
cat("\n",RR,"\n",KK," - covariate \n")
cat("\n","Analysis of Variance\n")
#cat("\n")
#cat("\n","Source         d.f.   Sum of Squares  Mean Squares         F  p.value")
#cat("\n",rep("-",35))
}
fuentes<- c( "TOTAL       ","GENOTYPES", "ENVIRONMENTS","INTERACTION",
             "HETEROGENEITY" , "RESIDUAL","POOLED ERROR")
gl <- c(A*M-1, A-1, M-1, (M - 1) * (A - 1), A - 1,(A - 1) * (M - 2),M * (A - 1) * (N - 1))
SC <- round(c(SKT, SKV, SKM , SKVM ,SKH, SKMB),4)
SC <- c(SC,0)

CM <- round(c(MKV, MKM, MKVM ,MKH ,MS,MKE),4)
CM <- c(0,CM)

Fcal<-round(c(FV,FM,FVM,FH,FMS),2)
Fcal<-c(0,Fcal,0)

resul<-c(" ",DD,NNN,LL,HH,BB," ")

Z<-data.frame(gl,SC,CM,Fcal, resul)
names(Z)<-c("d.f.","Sum of Squares","Mean Squares","F","p.value")
rownames(Z)<- c( "TOTAL       ","GENOTYPES", "ENVIRONMENTS","INTERACTION",
		"HETEROGENEITY" , "RESIDUAL","POOLED ERROR")
Z[7,c(2,4)]<-" ";Z[1,c(3,4)]<-" ";
if(console){
cat("\n")
print(Z)
}
Z1<-Z # Analysis
#Z<-as.matrix(Z)

#Z[1,4:6]   <-" "
#Z[7,c(3,5,6)]<-"       "
#for ( i in 1: 7) {
#cat("\n", Z[i,1],"\t",Z[i,2],"\t", Z[i,3],"\t", Z[i,4], "\t",Z[i,5], Z[i,6])
#}
#cat("\n",rep("-",35))

for ( i in 1: A) {
for ( j in 1: M) {
X1[i] <- X1[i] + y[i, j]
}
X1M[i] <- X1[i] / M
}
for ( j in 1: M) {
for ( i in 1: A) {
X2[j] <- X2[j] + y[i, j]
}
X2M[j] <- X2[j] / A
}
SH <- 0
for ( j in 1: M) {
SH <- SH + X2[j]
}
MM1 <- SH / (A * M)
for ( i in 1: A) {
for ( j in 1: M) {
W[i] <- W[i] + (y[i, j] - X1M[i] - X2M[j] + MM1) ^ 2 * N
}
}
MV<-round(MV,6); SI<-round(SI,6); S<-round(S,6); W<-round(W,6);
Z<-data.frame(MV, SI, NN, S,MMM,W)
names(Z)<-c("Mean","Sigma-square",".","s-square",".","Ecovalence")
rownames(Z)<-rownames(y)
if(console){
cat("\n","Genotype. Stability statistics\n\n")
print(Z)
}
Z2<-Z # statistics
Z<-as.matrix(Z)
#for ( i in 1: A) {
#cat("\n",i, "\t", Z[i,1],"\t", Z[i,2], Z[i,3], "\t",Z[i,4], Z[i,5],"\t", Z[i,6])
#}
FF <- SI / MKE # each genotype
SHF2 <- M * (A - 1) * (N - 1); SHF1 <- M - 1
F05 <- qf(0.95,SHF1, SHF2)
F01 <- qf(0.99,SHF1, SHF2)
for ( i in 1: A) {
if( FF[i] < FM0) F1[i] <- 0
if( FF[i] >= FM0) F1[i] <- -2
if( FF[i] >= F05) F1[i] <- -4
if( FF[i] >= F01) F1[i] <- -8
}
for ( i in 1: A) {
R0 <- 1
for ( j in 1: A) {
if( MV[j] < MV[i]) R0 <- R0 + 1
}
R[i] <- R0
}
for ( i in 1: A) GY[i] <- R[i] + MV1[i]
for ( i in 1: A) GYS[i] <- GY[i] + F1[i]
SGYS<-0
for ( i in 1: A) SGYS <- SGYS + GYS[i]
MGYS <- SGYS / A
for ( i in 1: A) {
GYY[i]<-""
if( GYS[i] > MGYS ) GYY[i] <- "+"
}
names<-c("Yield","Rank","Adj.rank","Adjusted","Stab.var","Stab.rating","YSi","..." )
Z<-data.frame(MV, R, MV1,GY,SI,F1,GYS,GYY)
rownames(Z)<-rownames(y)
names(Z)<-names
Z3<-Z # stability
if(console){
cat("\n\nSignif. codes:  0 '**' 0.01 '*' 0.05 'ns' 1\n\n")
cat("Simultaneous selection for yield and stability  (++)\n\n")
print(Z)
cat("\n","Yield Mean:", MES)
cat("\n","YS    Mean:", MGYS)
cat("\n","LSD (0.05):", DMV05)
cat("\n",rep("-",11))
cat("\n","+   selected genotype"                                            )
cat("\n","++  Reference: Kang, M. S. 1993. Simultaneous selection for yield")
cat("\n","and stability: Consequences for growers. Agron. J. 85:754-757."   )
cat("\n")
}
out<-list(analysis=Z1,statistics =Z2,stability =Z3)
invisible(out)
}


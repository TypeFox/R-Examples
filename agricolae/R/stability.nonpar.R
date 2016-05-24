`stability.nonpar` <-
function(data, variable= NULL, ranking = FALSE,console=FALSE)
{
row.names(data)<-data[,1]
data<-data[,-1]
audpc<-as.matrix(data)
rnames<-rownames(data)
cnames<-colnames(data)
n<-ncol(audpc); k<-nrow(audpc)
lineas<-matrix(0,k)
ldev <-matrix(0,k)
s2 <- matrix(0,k)
s1 <- matrix(0,k)
mr <- matrix(0,k)
ambientes<-matrix(0,1,n)
Maudpc <- mean(audpc)
for(i in 1:k) {
lineas[i]<-mean(audpc[i,])
ldev[i] <-lineas[i] - Maudpc }
cyld<-matrix(0,k,n)
rcyld<-matrix(0,k,n)
raudpc<-matrix(0,k,n)
dimnames(raudpc)<-list(rnames,cnames)
for(j in 1:n) ambientes[j]<-mean(audpc[,j])

for(i in 1:k){
for(j in 1:n) cyld[i,j] <- audpc[i,j] - ldev[i] }
for(j in 1:n) {
rcyld[,j]<-rank(cyld[,j])
raudpc[,j]<-rank(audpc[,j])}
#=========================================
nn=n*(n-1)/2
for (i in 1:k) {
sumq<-0
for(j in 1:(n-1)) for(jj in (j+1):n) sumq <- sumq+abs(rcyld[i,j]-rcyld[i,jj])
s1[i]=sumq/nn
}
#=========================================
for(i in 1:k) {
s2[i]<-var(rcyld[i,])
# mr[i]<-mean(rcyld[i,]) 
}

#Calculation of expectation and variance

k2<-k**2;                     
es1<-(k2-1)/(3*k) 

es2<-(k2-1)/12

vs1<-k2*((k**2-4)*(n+3)+30)/(45*k**2*n*(n-1))
vs2<-(k2-1)*(2*(k2-4)*(n-1)+5*(k2-1))/(360*n*(n-1))
alphas<-1-0.05/k
chi.ind<-qchisq(alphas,1)
chi.sum<-qchisq(0.95,k)
#=========================================
z1<-(s1-es1)**2/vs1
z2<-(s2-es2)**2/vs2
suma.z1 <- sum(z1)
suma.z2 <- sum(z2)
mr<-rank(lineas)
stat1 <- data.frame(lineas,mr,s1,z1,s2,z2)
row.names(stat1)<-rnames
names(stat1)<- c("Mean","Rank", "s1", "Z1", "s2", "Z2")
stat2 <- data.frame(MEAN=Maudpc,es1,es2,vs1,vs2,chi.ind,chi.sum)
if(console){
cat("\n")
cat("Nonparametric Method for Stability Analysis\n")
cat("-------------------------------------------\n\n")
cat("Estimation and test of nonparametric measures\n")
cat("Variable:",variable,"\n\n")
if (ranking) {
cat("Ranking...\n")
print(raudpc)
cat("\n")
}
cat("Statistics...\n")
print(round(stat1,2))
cat("------------------------")
cat("\nSum of Z1: ",suma.z1)
cat("\nSum of Z2: ",suma.z2)
cat("\n------------------------\n\n")
cat("Test...\n")
cat("The Z-statistics are measures of stability. The test for the significance\n")
cat("of the sum of Z1 or Z2 are compared to a Chi-Square value of chi.sum. \n")
cat("individual Z1 or Z2 are compared to a Chi-square value of chi.ind.\n\n")
print(stat2)
cat("---\n")
cat("expectation and variance: es1, es2, vs1, vs2\n")
}
out<-list(ranking=stat1,statistics =stat2)
invisible(out)
}


# ------------------------------------------------------------------------------------- #
# ------------------------------------------------------------------------------------- #
#      REGIONAL FREQUENCY ANALYSIS OF ANNUAL FLOWS IN PIEMONTE AND VALLE D'AOSTA        #
# ------------------------------------------------------------------------------------- #
# ------------------------------------------------------------------------------------- #

# ------------------------------------------------------------------------------------- #
#                        Regionalization of the growth-curve                            #
# ------------------------------------------------------------------------------------- #


# Data loading

data(hydroSIMN)

str(annualflows)


# Plot consistency of data series

y <- annualflows["anno"][,]
cod <- annualflows["cod"][,]
consistencyplot(y,cod)
par(ask = interactive())

# L-moments of data series

D <- annualflows["dato"][,]
lmSIMN <- data.frame(t(sapply(split(D,cod),Lmoments)))
print(lmSIMN)


# L-moments ratio plots

plot(lmSIMN[3:5])


# Choice of sites with more than 15 records

Dlen <- tapply(D,cod,length)
annualflows15 <- annualflows[unsplit(Dlen,cod)>=15,]
parameters15 <- parameters[Dlen>=15,]
D15 <- annualflows15["dato"][,]
cod15 <- annualflows15["cod"][,]


# L-moments ratio plot

D15lcv <- tapply(D15,cod15,LCV)
D15lca <- tapply(D15,cod15,LCA)
plot(D15lca,D15lcv,xlab="L-CA",ylab="L-CV",pch=1,col=1,main="L-moments space",cex.main=1,font.main=1)
grid()


# Choice of the classification variables:
#  create a function using the 'leaps' function of package 'subselect' 
#  to perform all-possible-regressions
# bestregressions <- function(dip,ind) {
#  Y <- as.numeric(dip)
#  X <- ind
#  Sy <- var(Y)
#  Sx <- var(X)
#  Sxy <- var(X,Y)
#  Dm.mat <- Sx
#  Dm.H <- Sxy %*% t(Sxy)/Sy
#  require(subselect)
#  Dm.leaps <- leaps(Dm.mat, kmin=1, kmax=3, H=Dm.H, r=1, nsol=3)
#  Dm.leaps
#  for(i in 3:1) {for(j in 1:3) {print(colnames(X)[Dm.leaps$subsets[j,c(1:3),i]])}}
# }

# bestregressions(as.numeric(D15lcv),parameters15[,3:16])
# [1] "S2000" "Rc"    "IB"
# [1] "Am"    "S2000" "Rc"
# [1] "Am"    "LLDP"  "S2000"
# [1] "S2000" "Ybar"
# [1] "Hm"   "Ybar"
# [1] "Am"    "S2000"
# [1] "S2000"
# [1] "Hm"
# [1] "Ybar"


# Choice of the classification variables:
#  or reasoning with distance matrices:
# bestregressions(as.numeric(AD.dist(D15,cod15)),data.frame(apply(parameters15[,3:16],2,dist)))
# [1] "Am"    "S2000" "Ybar"
# [1] "S2000" "EST"   "Ybar"
# [1] "Pm"    "S2000" "Ybar"
# [1] "S2000" "Ybar"
# [1] "Hm"   "Ybar"
# [1] "S2000" "EST"
# [1] "S2000"
# [1] "Hm"
# [1] "Ybar"


# We choose Hm and Ybar as classification variables.
# Mantel test:

Y <- AD.dist(D15,cod15)
X <- data.frame(apply(parameters15[,c("Hm","Ybar")],2,dist))
datamantel <- cbind(as.numeric(Y),X)
regrmantel <- lm(Y ~ Hm + Ybar, datamantel)
print(summary(regrmantel))

print(mantel.lm(regrmantel, Nperm=100))


# We choose Hm and Ybar as classification variables.
# Cluster formation:

param <- parameters15[c("Hm","Ybar")]
n <- dim(param)[1]; k <- dim(param)[2]
param.norm <- (param - matrix(colMeans(param),nrow=n,ncol=k,byrow=TRUE))/matrix(sapply(param, sd),nrow=n,ncol=k,byrow=TRUE)
clusters <- traceWminim(param.norm,4);
names(clusters) <- parameters15["cod"][,]
print(clusters)


# Plot of clusters:

par(mfrow=c(2,2))
plot(parameters15[c("Hm","Ybar")],col=clusters,pch=clusters,cex=0.6,
  main="Clusters in morphologic parameters space",cex.main=1,font.main=1)
points(tapply(parameters15["Hm"][,],clusters,mean),tapply(parameters15["Ybar"][,],clusters,mean),
  col=c(1:4),pch=c(1:4))
legend("topleft",paste("clust ",c(1:4)),col=c(1:4),pch=c(1:4),bty="n")
grid()

plot(parameters15[c("Xbar","Ybar")],col=clusters,pch=clusters,
  main="Clusters in geographical space",cex.main=1,font.main=1)
#legend("topleft",paste("clust ",c(1:4)),col=c(1:4),pch=c(1:4),bty="n")
grid()

plot(D15lca,D15lcv,xlab="L-CA",ylab="L-CV",pch=clusters,col=clusters,
  main="Clusters in L-moments space",cex.main=1,font.main=1)
#legend("topleft",paste("clust ",c(1:4)),col=c(1:4),pch=c(1:4),bty="n")
grid()

D15lkur <- tapply(D15,cod15,Lkur)
plot(D15lca,D15lkur,xlab="L-CA",ylab="L-kur",pch=clusters,col=clusters,
  main="Clusters in L-moments space",cex.main=1,font.main=1)
#legend("topleft",paste("clust ",c(1:4)),col=c(1:4),pch=c(1:4),bty="n")
grid()
par(mfrow=c(1,1))


# Homogeneity test (Cluster 1):

fac <- factor(annualflows15["cod"][,],levels=names(clusters[clusters==1]))
D15.1 <- annualflows15[!is.na(fac),"dato"]
cod15.1 <- factor(annualflows15[!is.na(fac),"cod"])
D15.1med <- tapply(D15.1,cod15.1,mean)
D15.1adim <- D15.1/unsplit(D15.1med,cod15.1)
regionalplotpos(D15.1adim,cod15.1,xlab="cluster1",main="Empirical distributions",cex.main=1,font.main=1)
HWs <- HW.tests(D15.1,cod15.1); HWs
mtext(paste("HW = ",round(HWs[1],2)),3,-1.5,adj=0.05)
nomi <- row.names(parameters15[clusters==1,])
legend("bottomright",nomi,col=c(1:length(nomi)),pch=c(1:length(nomi)),bty="n")


# Homogeneity test (Cluster 2):

fac <- factor(annualflows15["cod"][,],levels=names(clusters[clusters==2]))
D15.2 <- annualflows15[!is.na(fac),"dato"]
cod15.2 <- factor(annualflows15[!is.na(fac),"cod"])
D15.2med <- tapply(D15.2,cod15.2,mean)
D15.2adim <- D15.2/unsplit(D15.2med,cod15.2)
regionalplotpos(D15.2adim,cod15.2,xlab="cluster2",main="Empirical distributions",cex.main=1,font.main=1)
HWs <- HW.tests(D15.2,cod15.2); HWs
mtext(paste("HW = ",round(HWs[1],2)),3,-1.5,adj=0.05)
nomi <- row.names(parameters15[clusters==2,])
legend("bottomright",nomi,col=c(1:length(nomi)),pch=c(1:length(nomi)),bty="n")


# Homogeneity test (Cluster 3):

fac <- factor(annualflows15["cod"][,],levels=names(clusters[clusters==3]))
D15.3 <- annualflows15[!is.na(fac),"dato"]
cod15.3 <- factor(annualflows15[!is.na(fac),"cod"])
D15.3med <- tapply(D15.3,cod15.3,mean)
D15.3adim <- D15.3/unsplit(D15.3med,cod15.3)
regionalplotpos(D15.3adim,cod15.3,xlab="cluster3",main="Empirical distributions",cex.main=1,font.main=1)
HWs <- HW.tests(D15.3,cod15.3); HWs
mtext(paste("HW = ",round(HWs[1],2)),3,-1.5,adj=0.05)
nomi <- row.names(parameters15[clusters==3,])
legend("bottomright",nomi,col=c(1:length(nomi)),pch=c(1:length(nomi)),bty="n")


# Homogeneity test (Cluster 4):

fac <- factor(annualflows15["cod"][,],levels=names(clusters[clusters==4]))
D15.4 <- annualflows15[!is.na(fac),"dato"]
cod15.4 <- factor(annualflows15[!is.na(fac),"cod"])
D15.4med <- tapply(D15.4,cod15.4,mean)
D15.4adim <- D15.4/unsplit(D15.4med,cod15.4)
regionalplotpos(D15.4adim,cod15.4,xlab="cluster4",main="Empirical distributions",cex.main=1,font.main=1)
HWs <- HW.tests(D15.4,cod15.4); HWs
mtext(paste("HW = ",round(HWs[1],2)),3,-1.5,adj=0.05)
nomi <- row.names(parameters15[clusters==4,])
legend("bottomright",nomi,col=c(1:length(nomi)),pch=c(1:length(nomi)),bty="n")


# Homogeneity test (justification of the test choise):

lmom4reg <- data.frame(rbind(Lmoments(D15.1adim),Lmoments(D15.2adim),Lmoments(D15.3adim),Lmoments(D15.4adim)))
Lspace.HWvsAD()
points(lmom4reg[,c(4,3)])
text(lmom4reg[,c(4,3)],row.names(lmom4reg),adj=c(-.5,-.5))


# Model selection (L-moments ratio diagram):

Lmoment.ratio.diagram()
points(Lmoments(D15.1adim)[4],Lmoments(D15.1adim)[5],pch=1,col=1)
points(Lmoments(D15.2adim)[4],Lmoments(D15.2adim)[5],pch=2,col=2)
points(Lmoments(D15.3adim)[4],Lmoments(D15.3adim)[5],pch=3,col=3)
points(Lmoments(D15.4adim)[4],Lmoments(D15.4adim)[5],pch=4,col=4)
legend("bottomleft",paste("clust ",c(1:4)),col=c(1:4),pch=c(1:4),bty="n")


# Model selection (Cluster 1):

Lmom.1 <- Lmoments(D15.1adim); print(Lmom.1)

par1 <- par.gamma(Lmom.1["l1"],Lmom.1["l2"],Lmom.1["lca"]);
print(par1)

F1 <- F.gamma(D15.1adim,par1$xi,par1$beta,par1$alfa)
regionalplotpos(D15.1adim,cod15.1,xlab="cluster1",main="Empirical distributions",cex.main=1,font.main=1)
lines(sort(D15.1adim),sort(F1))
AD <- gofP3test(D15.1adim,Nsim=100)
ADb <- A2_GOFlaio(D15.1adim, dist="GAM")
AD2 <- gofGEVtest(D15.1adim,Nsim=100)
AD2b <- A2_GOFlaio(D15.1adim, dist="GEV")
AD3 <- gofGENLOGIStest(D15.1adim,Nsim=100)
AD4 <- gofGENPARtest(D15.1adim,Nsim=100)
AD5 <- gofLOGNORMtest(D15.1adim,Nsim=100)
mtext(paste("P(AD) P3 = ",round(AD[2],3),", ",round(ADb[2],3)),3,-1.5,adj=0.05)
mtext(paste("P(AD) GEV = ",round(AD2[2],3),", ",round(AD2b[2],3)),3,-2.5,adj=0.05)
mtext(paste("P(AD) GL = ",round(AD3[2],3)),3,-3.5,adj=0.05)
mtext(paste("P(AD) GP = ",round(AD4[2],3)),3,-4.5,adj=0.05)
mtext(paste("P(AD) LN = ",round(AD5[2],3)),3,-5.5,adj=0.05)


# Model selection (Cluster 2):

Lmom.2 <- Lmoments(D15.2adim); print(Lmom.2)

par2 <- par.gamma(Lmom.2["l1"],Lmom.2["l2"],Lmom.2["lca"])
print(par2)

F2 <- F.gamma(D15.2adim,par2$xi,par2$beta,par2$alfa)
regionalplotpos(D15.2adim,cod15.2,xlab="cluster2",main="Empirical distributions",cex.main=1,font.main=1)
lines(sort(D15.2adim),sort(F2))
AD <- gofP3test(D15.2adim,Nsim=100)
ADb <- A2_GOFlaio(D15.2adim, dist="GAM")
AD2 <- gofGEVtest(D15.2adim,Nsim=100)
AD2b <- A2_GOFlaio(D15.2adim, dist="GEV")
AD3 <- gofGENLOGIStest(D15.2adim,Nsim=100)
AD4 <- gofGENPARtest(D15.2adim,Nsim=100)
AD5 <- gofLOGNORMtest(D15.2adim,Nsim=100)
mtext(paste("P(AD) P3 = ",round(AD[2],3),", ",round(ADb[2],3)),3,-1.5,adj=0.05)
mtext(paste("P(AD) GEV = ",round(AD2[2],3),", ",round(AD2b[2],3)),3,-2.5,adj=0.05)
mtext(paste("P(AD) GL = ",round(AD3[2],3)),3,-3.5,adj=0.05)
mtext(paste("P(AD) GP = ",round(AD4[2],3)),3,-4.5,adj=0.05)
mtext(paste("P(AD) LN = ",round(AD5[2],3)),3,-5.5,adj=0.05)


# Model selection (Cluster 3):

Lmom.3 <- Lmoments(D15.3adim); print(Lmom.3)

par3 <- par.gamma(Lmom.3["l1"],Lmom.3["l2"],Lmom.3["lca"])
print(par3)

F3 <- F.gamma(D15.3adim,par3$xi,par3$beta,par3$alfa)
regionalplotpos(D15.3adim,cod15.3,xlab="cluster3",main="Empirical distributions",cex.main=1,font.main=1)
lines(sort(D15.3adim),sort(F3))
AD <- gofP3test(D15.3adim,Nsim=100)
ADb <- A2_GOFlaio(D15.3adim, dist="GAM")
AD2 <- gofGEVtest(D15.3adim,Nsim=100)
AD2b <- A2_GOFlaio(D15.3adim, dist="GEV")
AD3 <- gofGENLOGIStest(D15.3adim,Nsim=100)
AD4 <- gofGENPARtest(D15.3adim,Nsim=100)
AD5 <- gofLOGNORMtest(D15.3adim,Nsim=100)
mtext(paste("P(AD) P3 = ",round(AD[2],3),", ",round(ADb[2],3)),3,-1.5,adj=0.05)
mtext(paste("P(AD) GEV = ",round(AD2[2],3),", ",round(AD2b[2],3)),3,-2.5,adj=0.05)
mtext(paste("P(AD) GL = ",round(AD3[2],3)),3,-3.5,adj=0.05)
mtext(paste("P(AD) GP = ",round(AD4[2],3)),3,-4.5,adj=0.05)
mtext(paste("P(AD) LN = ",round(AD5[2],3)),3,-5.5,adj=0.05)


# Model selection (Cluster 4):

Lmom.4 <- Lmoments(D15.4adim); print(Lmom.4)

par4 <- par.gamma(Lmom.4["l1"],Lmom.4["l2"],Lmom.4["lca"])
print(par4)

F4 <- F.gamma(D15.4adim,par4$xi,par4$beta,par4$alfa)
regionalplotpos(D15.4adim,cod15.4,xlab="cluster4",main="Empirical distributions",cex.main=1,font.main=1)
lines(sort(D15.4adim),sort(F4))
AD <- gofP3test(D15.4adim,Nsim=100)
ADb <- A2_GOFlaio(D15.4adim, dist="GAM")
AD2 <- gofGEVtest(D15.4adim,Nsim=100)
AD2b <- A2_GOFlaio(D15.4adim, dist="GEV")
AD3 <- gofGENLOGIStest(D15.4adim,Nsim=100)
AD4 <- gofGENPARtest(D15.4adim,Nsim=100)
AD5 <- gofLOGNORMtest(D15.4adim,Nsim=100)
mtext(paste("P(AD) P3 = ",round(AD[2],3),", ",round(ADb[2],3)),3,-1.5,adj=0.05)
mtext(paste("P(AD) GEV = ",round(AD2[2],3),", ",round(AD2b[2],3)),3,-2.5,adj=0.05)
mtext(paste("P(AD) GL = ",round(AD3[2],3)),3,-3.5,adj=0.05)
mtext(paste("P(AD) GP = ",round(AD4[2],3)),3,-4.5,adj=0.05)
mtext(paste("P(AD) LN = ",round(AD5[2],3)),3,-5.5,adj=0.05)


# Comparison between regional growth-curves:

Fs <- seq(0.001,0.999,by=.001)
q1 <- invF.gamma(Fs,par1$xi,par1$beta,par1$alfa)
q2 <- invF.gamma(Fs,par2$xi,par2$beta,par2$alfa)
q3 <- invF.gamma(Fs,par3$xi,par3$beta,par3$alfa)
q4 <- invF.gamma(Fs,par4$xi,par4$beta,par4$alfa)

lognormplot(c(q1,q2,q3,q4),line=FALSE,type="n")
normpoints(q1,type="l",lty=1,col=1)
normpoints(q2,type="l",lty=2,col=2)
normpoints(q3,type="l",lty=3,col=3)
normpoints(q4,type="l",lty=4,col=4)
legend("bottomright",paste("cluster ",c(1:4)),col=c(1:4),lty=c(1:4),bty="n")


# Final tests of goodness-of-fit. 

# The regionalized distribution adaptation, with mean 
#   D.est = exp(7.86 + 0.000291*Hm + 0.0722*NORD - 1.70*IB)
# and parameters of Pearson type III distribution
#   xi1 = 0.5458;  beta1 = 0.09534; alfa1 = 4.764;
#   xi2 = 0.09843; beta2 = 0.08508; alfa2 = 10.60;
#   xi3 = 0.2801;  beta3 = 0.1237;  alfa3 = 5.817;
#   xi4 = 0.3496;  beta4 = 0.2093;  alfa4 = 3.107
# is tested throught the Anderson-Darling and the Cramer-von Mises
# goodness-of-fit tests.
# In this case (fully specified distribution) the percentage points
# are:
#        .25    .10    .05    .01   .001
#   A2 1.248  1.933  2.492  3.880  6.000
#   W2 0.209  0.347  0.461  0.743  1.167  \n

A2fullspecif <- function (x,xi,beta,alfa) {
 x <- sort(x)
 n <- length(x)
 F <- pgamma((x - xi)/beta, alfa)
 F[F<=0] <- 0.000001
 F[F>=1] <- 0.999999
 -n-(1/n)*sum((seq(1,2*n-1,by=2))*log(F) + (seq(2*n-1,1,by=-2))*log(1-F))
} 
W2fullspecif <- function (x,xi,beta,alfa) {
 x <- sort(x)
 n <- length(x)
 F <- pgamma((x - xi)/beta, alfa)
 sum((F - seq(1,2*n-1,by=2)/(2*n))^2) + 1/(12*n)
}
par(mfrow=c(3,2))
W2s <- rep(NA,47)
A2s <- rep(NA,47)
em <- rep(NA,47)
gruppi <- rep(NA,47)
ferma=6
for(i in 1:47) {
 codice <- unique(cod)[i]
 Hm <- parameters[parameters["cod"]==codice,"Hm"]
 NORD <- parameters[parameters["cod"]==codice,"NORD"]
 IB <- parameters[parameters["cod"]==codice,"IB"]
 Ybar <- parameters[parameters["cod"]==codice,"Ybar"]
 baricentri.clusterD <- cbind(c(2327.455,1404.200,1892.214,812.125),c(45.77200,46.01480,44.58071,44.47513))
 baricentri.clusterD[,1] <- (baricentri.clusterD[,1]-1694)/652.1
 baricentri.clusterD[,2] <- (baricentri.clusterD[,2]-45.12)/0.7229
 HmYbar <- c((Hm-1694)/652.1,(Ybar-45.12)/0.7229)
 distanze <- as.matrix(dist(rbind(HmYbar,baricentri.clusterD)))[-1,1]
 gruppo <- which.min(distanze)
 em[i] <- exp(7.86 + 0.000291*Hm + 0.0722*NORD - 1.70*IB)
 if(gruppo==1){
  xi=0.5458; beta=0.09534; alfa=4.764;
 }
 else if(gruppo==2){
  xi=0.09843; beta=0.08508; alfa=10.60;
 }
 else if(gruppo==3){
  xi=0.2801; beta=0.1237; alfa=5.817;
 }
 else if(gruppo==4){
  xi=0.3496; beta=0.2093; alfa=3.107;
 }
 gruppi[i] <- gruppo
 campione <- annualflows[annualflows["cod"]==codice,3]
 plotpos(campione,main=row.names(parameters[parameters["cod"]==codice,]),cex.main=1,font.main=1)
 lines(sort(campione),pgamma((sort(campione)-em[i]*xi)/(em[i]*beta),alfa))
 A2s[i] <- A2fullspecif(campione,em[i]*xi,em[i]*beta,alfa)
 W2s[i] <- W2fullspecif(campione,em[i]*xi,em[i]*beta,alfa)
 mtext(paste("A2 = ",round(A2s[i],3)),3,-1.5,adj=0.05) 
 mtext(paste("W2 = ",round(W2s[i],3)),1,-1.5,adj=0.95)
 if(i==ferma) {
  ferma <- ferma+6
  #readline("Press Return to continue to the next graphs")
 }
}
#tabella <- cbind(parameters[1],round(tapply(D,cod,mean)),round(em),gruppi,A2s,W2s)
tabella <- cbind(parameters[1:2],round(em),gruppi,A2s,W2s)
names(tabella) <- c("cod","tildeDm","hatDm","group","A2","W2")
print(tabella) 
par(mfrow=c(1,1))


# Dependency on the estimation of the mean.

errors <- abs((tabella[,2] - tabella[,3])/tabella[,2])
codici <- tabella[,1]
plot(A2s,errors,main="",xlab=expression(A^2),ylab=expression(paste("|(",hat(D)[m] - tilde(D)[m],")/",tilde(D)[m],"|")))
grid()
abline(v=2.492,lty=2,col=2)
abline(v=3.880,lty=2,col=2)
text(A2s[A2s>3.880],errors[A2s>3.880],labels=codici[A2s>3.880]-2,pos=1)
axis(1,at=c(2.492,3.880),labels=c("5\\%","1\\%"),col=2,col.axis=2,cex.axis=.8)

# ------------------------------------------------------------------------------------- #
# ------------------------------------------------------------------------------------- #
#                                       THE END                                         #
# ------------------------------------------------------------------------------------- #
# ------------------------------------------------------------------------------------- #


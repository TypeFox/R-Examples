CAvariants <-
function(
  Xtable, mj=NULL,mi=NULL, printdims=4,firstaxis=1,lastaxis=2,
  catype = "CA",digits=4) { 
if (printdims<1) stop(paste("Attention: number of dims for output must be at least 1\n\n"))
if (lastaxis<2) stop(paste("Attention: last axis must be at least 2\n\n"))
if (!any(catype==c("CA","SOCA","DOCA","NSCA","SONSCA","DONSCA"))) stop(paste("Must be CA, DOCA, SOCA, NSCA, SONSCA or DONSCA"))

#if (!any(is.wholenumber(Xtable))) stop(paste("Must be integer values in contingency table"))

# READ DATA FILE
#  assume for now header and row names exist

#Xtable <- read.table(file = datafile, header=header)

#if (header==FALSE) { 
#for (i in 1:dim(Xtable)[1]) rownames(Xtable)[i] <- paste("r",i,sep="")
#for (i in 1:dim(Xtable)[2]) colnames(Xtable)[i] <- paste("c",i,sep="")
#}

X <- as.matrix(Xtable)
#M<-min(nrow(Xtable),ncol(Xtable))-1
  rowlabels <- rownames(Xtable)
  collabels <- colnames(Xtable) 

rows <- dim(X)[1]
cols <- dim(X)[2]
n <- sum(X)
if (is.null(mj)){ 
   mj <- c(1:cols)}
else 
mj<-c(mj) #natural scores for columns  
if (is.null(mi)){ 
   mi <- c(1:rows)}
else 
mi<-c(mi)  #natural scores for rows
maxaxes <- min(rows,cols)-1
r<-maxaxes
S <- switch(catype, "CA"=cabasic(X),  "SOCA"=socabasic(X,mj),"DOCA"=docabasic(X,mi,mj),"NSCA"=nscabasic(X),"SONSCA"=sonscabasic(X,mj),
"DONSCA"=donscabasic(X,mi,mj))

##########################------CA
if(catype=="CA"){
Fmat <- S@RX %*% S@Rweights %*% S@Raxes
Gmat <- S@CX %*% S@Cweights %*% S@Caxes
#dmum1 <- diag( (S@mu + (S@mu==0)) * (1-(S@mu==0)) )
Fbi <- S@Cweights %*%  S@Caxes # no orthonormal
Gbi <-S@Rweights %*%   S@Raxes # no orthonormal
pcc <- t(S@CX)
#dimnames(pcc)<-dimnames(X)
tau=NULL
tauden=NULL
inertia <- (S@mu[1:r]*S@mu[1:r])
comps<-diag(inertia)
Trend<-(Fmat[,firstaxis:lastaxis]%*%t(Gbi[,firstaxis:lastaxis]))
Z<-Trend
}
#########################################################---------DOCA
if(catype=="DOCA"){
Z<-S@Z/sqrt(n)
pcc<-S@RX #centered column profile matrix
Gbi <- S@Raxes
Fbi<- S@Caxes 
Gmat <-  S@CX %*%  S@Caxes #row principal coordinates
Fmat <- S@RX  %*% S@Raxes #column principal coordinates
if  ((Z[1,1]<0) & (Z[1,2]>0)||(Z[1,1]>0) & (Z[1,2]<0)){
Gmat<-(-1)*Gmat
Fmat<-(-1)*Fmat
}

reconstruction<-t(Gmat%*%t(S@Cweights%*%Fbi))
dimnames(reconstruction)<-dimnames(X)
inertia <- S@mu[1:r]/n #number of inertia of row poly
inertia2<-S@mu2[1:r]/n #number of inertias of column poly
Z1<-S@Z
tau=NULL
tauden=NULL
comps <- compstable.exe(Z1) 
Icompnames <- c( "Location", "Dispersion", "Cubic","Error", "** Chi-squared Statistic **")
Jcompnames <- c("Component Value", "P-value")
dimnames(Z) <- list(paste("u", 1:(rows - 1),sep=""), paste("v", 1:(cols - 1),sep=""))
dimnames(comps$compsR) <- list(paste(Icompnames), paste(Jcompnames))
dimnames(comps$compsC) <- list(paste(Icompnames), paste(Jcompnames))

Trend<-(Fmat[,firstaxis:lastaxis]%*%t(S@Rweights%*%Gbi[,firstaxis:lastaxis]))
#browser()
}
#########################################################################---------SOCA
if(catype=="SOCA"){
pcc<-S@RX
Z<-S@Z/sqrt(n)
dimnames(pcc)<-dimnames(X)
Gmat <- S@CX %*%  S@Caxes #row principal coordinates
Fmat <- S@RX %*% S@Rweights %*% S@Raxes #column principal coordinates
if ((Z[1,1]<0) & (Z[1,2]>0)||(Z[1,1]>0) & (Z[1,2]<0)){Gmat<-(-1)*Gmat}
Gbi <-S@Raxes
Fbi <- S@Caxes 
inertia <- (S@mu[1:r]^2)/n
inertia2<-(S@mu2[1:r])/n
#inertia2<-(S@mu2[-1])/n
#comps<-diag(inertia2)
tauden=NULL
tau=NULL
Z1<-S@Z
comps <- compsonetable.exe(Z1) 
#Icompnames <- c( "** Column Components **", "Location", "Dispersion", "Cubic","Error", "** Chi-squared Statistic **")
Icompnames <- c( "Location", "Dispersion","Cubic", "Error", "** Chi-squared Statistic **")
Jcompnames <- c("Component Value", "P-value")
dimnames(Z) <- list(paste("m", 1:nrow(Z),sep="" ), paste("v", 1:(cols - 1),sep=""))
dimnames(comps) <- list(paste(Icompnames), paste(Jcompnames))
Trend<-(Fmat[,firstaxis:lastaxis]%*%t(S@Rweights%*%Gbi[,firstaxis:lastaxis]))
#browser()
}

####################################-------------NSCA
if(catype=="NSCA"){
Fbi <-  S@Caxes
Gbi <-   S@Raxes 
#dmum1 <-diag( (S@mu + (S@mu==0)) * (1-(S@mu==0)) )
dmum1 <-diag( S@mu [1:r])
pcc<-S@RX
dimnames(pcc)<-dimnames(X)
Gmat <-  S@Raxes[,1:r] %*% dmum1
Fmat <- S@Caxes[,1:r] %*% dmum1
tauden<-S@tauDen
inertia <- S@mu[1:r]*S@mu[1:r]
tau<-sum(inertia)/tauden
comps<-diag(inertia)
Trend<-(Fmat[,firstaxis:lastaxis]%*%t(S@Rweights%*%Gbi[,firstaxis:lastaxis]))
Z<-Trend

}
##################################-------------DONSCA
if(catype=="DONSCA"){
Fbi <- S@Caxes
Gbi<-S@Raxes
pcc<-S@RX
Z<-S@Z
dimnames(pcc)<-dimnames(X)
Gmat <-  S@CX %*% S@Cweights %*% S@Caxes #row principal coordinates
Fmat <- S@RX  %*% S@Rweights %*% S@Raxes #column principal coordinates
if  ((Z[1,1]<0) & (Z[1,2]>0)||(Z[1,1]>0) & (Z[1,2]<0)){Gmat<-(-1)*Gmat}
inertia <- S@mu[1:r]
inertia2<-S@mu2[1:r] 
tauden<-S@tauDen
Z2<-1/sqrt(tauden)*sqrt((n-1)*(rows-1))*S@Z
#Z2<-sqrt((n-1)*(rows-1))*S@Z #when tau
tau<-sum(inertia)/tauden
comps <- compstable.exe(Z2) 
Icompnames <- c( "Location", "Dispersion","Cubic", "Error", "** C-Statistic **")
Jcompnames <- c("Component Value", "P-value")
dimnames(Z) <- list(paste("u", 1:(rows-1 ),sep=""), paste("v", 1:(cols -1),sep=""))
#dimnames(Z) <- list(paste("u", 1:(rows ),sep=""), paste("v", 1:(cols -1),sep=""))
dimnames(comps$compsR) <- list(paste(Icompnames), paste(Jcompnames))
dimnames(comps$compsC) <- list(paste(Icompnames), paste(Jcompnames))

Trend<-(Fmat[,firstaxis:lastaxis]%*%t(Gbi[,firstaxis:lastaxis]))

#browser()
}
############################################------------------SONSCA
if(catype=="SONSCA"){
pcc<-S@RX
dimnames(pcc)<-dimnames(X)
Z<-S@Z
Gmat <- S@CX %*% S@Cweights %*% S@Caxes #column principal coordinates with principal axes
Fmat <- S@RX %*% (S@Rweights) %*% S@Raxes #row principal coordinates with polys
Gbi <- S@Raxes
#Gbi<-(sqrt(S@Rweights))%*%S@Raxes
Fbi <- S@Caxes 
if ((Z[1,1]<0) & (Z[1,2]>0)||(Z[1,1]>0) & (Z[1,2]<0)){Gmat<-(-1)*Gmat
#Fbi<-(-1)*Fbi
}
inertia <- S@mu[1:r]
inertia2 <- S@mu2[1:r]
tauden<-S@tauDen
tau<-sum(inertia)/tauden
Z1<-1/sqrt(tauden)*sqrt((n-1)*(rows-1))*S@Z
comps <- compsonetable.exe(Z1) 
Icompnames <- c( "Location", "Dispersion", "Cubic","Error", "** C-Statistic **")
Jcompnames <- c("Component Value", "P-value")
dimnames(Z) <- list(paste("m", 1:nrow(S@Z),sep=""), paste("v", 1:(cols - 1),sep=""))
dimnames(comps) <- list(paste(Icompnames), paste(Jcompnames))
Trend<-t(Gmat[,firstaxis:lastaxis]%*%t(Fbi[,firstaxis:lastaxis]))

#browser()
}

##################################################################################################################
# OTHER CALCULATIONS

# Calc inertia sum
#dmum2 <- diag( 1/(inertia + (S@mu==0)) * (1-(S@mu==0)) )
inertiasum <- sum(inertia)
inertiapc <- 100*inertia/inertiasum
cuminertiapc <- cumsum(inertiapc)
inertiapc <- round(100*inertiapc)/100
cuminertiapc <- round(100*cuminertiapc)/100
inertias <- round(cbind(inertia,inertiapc,cuminertiapc),digits=digits)
#browser()

##########################################################
if((catype=="SOCA")|(catype=="SONSCA")|(catype=="DOCA")|(catype=="DONSCA")){
inertiasum2 <- sum(inertia2) #for column categories diag(Z'Z)
inertiapc2 <- 100*inertia2/inertiasum2
cuminertiapc2 <- cumsum(inertiapc2)
inertiapc2 <- round(100*inertiapc2)/100
cuminertiapc2 <- round(100*cuminertiapc2)/100
inertias2 <- round(cbind(inertia2,inertiapc2,cuminertiapc2),digits=digits)
}
else inertias2<-inertias
# Calc contributions and correlations

Xstd <- X/sum(X)
if ((catype=="CA")|(catype=="SOCA")|(catype=="DOCA")){
dr <- diag(rowSums(Xstd))}
else {uni<-rep(1,rows)
dr<-diag(uni)}
dc <- diag(colSums(Xstd))
dimnames(Trend)<-list(rowlabels,collabels)
#############################
#cacorpo<- new("cacorporateplus",S=S,
#DataMatrix=X, rows=rows, cols=cols, 
#rowlabels=rowlabels, collabels=collabels,
#Rprinccoord=Fmat, Cprinccoord=Gmat, Rstdcoord=Fbi, Cstdcoord=Gbi,
# inertiasum=inertiasum, inertias=inertias, inertias2=inertias2,comps=comps,
#  maxaxes=maxaxes,catype=catype,printdims=printdims,mj=mj,mi=mi,pcc=pcc,Jmass=dc,Imass=dr,
#Trend=Trend,Z=Z)

list(
DataMatrix=X, rows=rows, cols=cols, r=r,
rowlabels=rowlabels, collabels=collabels,
Rprinccoord=Fmat[,1:r], Cprinccoord=Gmat[,1:r], Rstdcoord=Fbi[,1:r], Cstdcoord=Gbi[,1:r],tauden=tauden,tau=tau,
 inertiasum=inertiasum, inertias=inertias, inertias2=inertias2,comps=comps,
  maxaxes=maxaxes,catype=catype,printdims=printdims,mj=mj,mi=mi,pcc=pcc,Jmass=dc,Imass=dr,
Trend=Trend,Z=Z)

}

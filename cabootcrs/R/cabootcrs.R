
cabootcrs <- function(
  xobject=NULL, datafile=NULL, 
  groupings=NULL, grouplabels=NULL, 
  plotsymbolscolours=c(19,"alldifferent",18,"alldifferent"), 
  othersmonochrome="black", crpercent=95,
  nboots=1000, resampledistn="Poisson", multinomialtype="whole",
  poissonzeronewmean=0, newzeroreset=0, printdims=4, lastaxis=2, 
  catype = "sca", bootstdcoords=FALSE ) { 

# please ignore forward references to MCA functions

if (nboots==0) {
  cat(paste("\n Standard Correspondence Analysis results only, no confidence regions calculated\n\n"))
} else { 
if (nboots<999) cat(paste("\n WARNING", nboots, "is too few bootstrap replicates for reliable results\n\n"))
if (!any(resampledistn==c("Poisson","multinomial","nonparametric","balanced"))) stop(paste("Resampling must be Poisson, multinomial, nonparametric or balanced\n\n"))
if (!any(multinomialtype==c("whole","rowsfixed","columnsfixed"))) 
    stop(paste("Multinomial resampling is either whole, rowsfixed or columnsfixed\n\n"))
if ( (resampledistn=="Poisson") & any(multinomialtype==c("rowsfixed","columnsfixed")) )
    cat(paste("\n WARNING can't fix rows or columns with Poisson resampling\n\n"))
if (poissonzeronewmean<0) stop(paste("Poisson mean for zero entries must be non-negative\n\n"))
if (!any(newzeroreset==c(0,1))) stop(paste("Zeros in sample can only be set to zero or one in bootstrap\n\n"))
}
if (printdims<2) stop(paste("number of dims for output must be at least 2\n\n"))
if (lastaxis<2) stop(paste("last axis must be at least 2\n\n"))
if (!any(catype==c("sca","mca","omca"))) stop(paste("Must be sca, mca or omca"))

# Most axes looked at by rearrange routine
maxrearrange <- 6

# READ DATA FILE OR DATA OBJECT
# add row and column names if needed, or if currently has default "1" and "V1" etc

if (is.null(xobject)) { 
Xtable <- read.table(file = datafile)
if ((row.names(Xtable)[[1]]=="1")&(names(Xtable)[[1]]=="V1")) { 
  for (i in 1:dim(Xtable)[1]) rownames(Xtable)[i] <- paste("r",i,sep="")
  for (i in 1:dim(Xtable)[2]) colnames(Xtable)[i] <- paste("c",i,sep="")
}
} else { 
Xtable <- as.data.frame(xobject)
if ((is.null(rownames(xobject)))|(row.names(Xtable)[[1]]=="1")) {
  for (i in 1:dim(Xtable)[1]) rownames(Xtable)[i] <- paste("r",i,sep="")
}
if ((is.null(colnames(xobject)))|(names(Xtable)[[1]]=="V1")) {
  for (i in 1:dim(Xtable)[2]) colnames(Xtable)[i] <- paste("c",i,sep="")
}
}

X <- as.matrix(Xtable)

if ( (dim(X)[2]>dim(X)[1]) & (catype=="sca") ) {
  X <- t(X) 
  rowlabels <- colnames(Xtable)
  collabels <- rownames(Xtable) 
} else {
  rowlabels <- rownames(Xtable)
  collabels <- colnames(Xtable) 
}

rows <- dim(X)[1]
cols <- dim(X)[2]
n <- sum(X)
axisvariances <- min(lastaxis,rows-1,cols-1) # most axes for which variances can be calculated

if ((lastaxis >= maxrearrange)&(cols >=maxrearrange+2)) 
  cat(paste("\n WARNING variances for the ", maxrearrange, "th axis and above may be unreliable\n\n"))
if (lastaxis>axisvariances) 
  cat(paste("\n WARNING there are only", axisvariances, "axes in the solution\n\n"))

# S <- switch(catype, "sca"=sca(X), "mca"=mca(X), "omca"=omca(X))

S <- switch(catype, "sca"=sca(X), "mca"=sca(X), "omca"=sca(X))

# zero rows/cols to be omitted from bootstrapping and summaries
# rows/cols with one non-zero element, will have variance zero
# rows/cols with two non-zero elements, will use chi^2 with df 1

zerorowS <- rowSums(X>0)==0
zerocolS <- colSums(X>0)==0
onerowS <- rowSums(X>0)==1
onecolS <- colSums(X>0)==1
tworowS <- rowSums(X>0)==2
twocolS <- colSums(X>0)==2


# BOOTSTRAPPING

# Initialise sums of squares matrices
#  (computationally poor and may cause problems if large, 
#   but should be OK as values will tend to be close to zero)

RSBsum <- matrix(0,rows,axisvariances)
RSBsumsq <- matrix(0,rows,axisvariances)
RSBsumcp <- array(0,c(rows,axisvariances,axisvariances))
CSBsum <- matrix(0,cols,axisvariances)
CSBsumsq <- matrix(0,cols,axisvariances)
CSBsumcp <- array(0,c(cols,axisvariances,axisvariances))
rownB <- rep(nboots,rows)
colnB <- rep(nboots,cols)
RSBvar <- RSBsum
CSBvar <- CSBsum
RSBcov <- RSBsumcp
CSBcov <- CSBsumcp
sameaxisorder <- 0

# Used when bootstrapping an indicator matrix in MCA

if (resampledistn=="nonparametric") {
bno <- matrix( sample( rep(1:rows,nboots), replace=TRUE ), rows, nboots)
} else {
if (resampledistn=="balanced") {
bno <- matrix( sample( rep(1:rows,nboots) ), rows, nboots)
}
}

for (b in 1:nboots) { 

if (resampledistn=="multinomial") {
  Xr <- vector("numeric",rows*cols)
    if (multinomialtype=="whole") { Xr <- rmultinom(1,n,X) }
    if (multinomialtype=="rowsfixed") { for (i in 1:rows) { Xr[ seq(i,(cols-1)*rows+i,by=rows) ] <- rmultinom(1,sum(X[i,]),X[i,]) } }
    if (multinomialtype=="columnsfixed") { for (i in 1:cols) { Xr[((i-1)*rows+1):(i*rows)] <- rmultinom(1,sum(X[,i]),X[,i]) } } 
} else {
if (resampledistn=="Poisson") {
  Xr <- rpois(rows*cols,X+poissonzeronewmean*(X==0))
} else { 
  Xr <- X[bno[,b],]
}
}

XB <- matrix(data=Xr,nrow=rows,ncol=cols)

if (newzeroreset==1) { XB <- (XB==0 & X>0)+XB } # set new zeros to 1 

# B <- switch(catype, "sca"=sca(XB), "mca"=mca(XB), "omca"=omca(XB))

B <- switch(catype, "sca"=sca(XB), "mca"=sca(XB), "omca"=sca(XB))

# Check for zero rows/cols

zerorowB <- rowSums(XB>0)==0
zerocolB <- colSums(XB>0)==0
rownB <- rownB-zerorowB
colnB <- colnB-zerocolB

# Rearrange bootstrap axes if needed, currently only up to 6

Re <- rearrange( S@Raxes, B@Raxes, S@Caxes, B@Caxes, B@r  )

RaxesBRe <- B@Raxes
CaxesBRe <- B@Caxes

RaxesBRe[ ,1:Re$numrearranged] <- RaxesBRe[ ,1:Re$numrearranged] %*% Re$T
CaxesBRe[ ,1:Re$numrearranged] <- CaxesBRe[ ,1:Re$numrearranged] %*% Re$T

RSB <- ( B@Rprofile - S@Rprofile ) %*% B@Rweights %*% RaxesBRe[ ,1:axisvariances]
CSB <- ( B@Cprofile - S@Cprofile ) %*% B@Cweights %*% CaxesBRe[ ,1:axisvariances]
RSB <- RSB * (1-zerorowS) * (1-zerorowB)
CSB <- CSB * (1-zerocolS) * (1-zerocolB)
sameaxisorder <- sameaxisorder + Re$same

# Unlikely to want to bootstrap standard coordinates in practice,
# but may be useful in research

if (bootstdcoords==TRUE) {
dmum1 <- diag( 1/(B@mu + (B@mu==0)) * (1-(B@mu==0)) ) 
dmum1[1:Re$numrearranged,1:Re$numrearranged] <- dmum1[1:Re$numrearranged,1:Re$numrearranged] %*% Re$T
dmum1 <- dmum1[1:axisvariances,1:axisvariances]
RSB <- RSB %*% dmum1
CSB <- CSB %*% dmum1
}

RSBsum <- RSBsum + RSB
RSBsumsq <- RSBsumsq + RSB*RSB
CSBsum <- CSBsum + CSB
CSBsumsq <- CSBsumsq + CSB*CSB
if (axisvariances>1) {
for (a1 in 1:(axisvariances-1)) {
  for (a2 in (a1+1):axisvariances) { 
    RSBsumcp[ ,a1,a2] <- RSBsumcp[ ,a1,a2] + RSB[ ,a1]*RSB[ ,a2]
    CSBsumcp[ ,a1,a2] <- CSBsumcp[ ,a1,a2] + CSB[ ,a1]*CSB[ ,a2]
} } }

} # boots

RSBmean <- diag( ( 1/( rownB + (rownB==0) ) ) * (1-(rownB==0)) ) %*% RSBsum
CSBmean <- diag( ( 1/( colnB + (colnB==0) ) ) * (1-(colnB==0)) ) %*% CSBsum

rbm1 <- diag( ( 1/( rownB-1 + 2*(rownB<=1) ) ) * (1-(rownB<=1)) )
cbm1 <- diag( ( 1/( colnB-1 + 2*(colnB<=1) ) ) * (1-(colnB<=1)) )

RSBvar <- rbm1 %*% ( RSBsumsq - diag(rownB) %*% RSBmean * RSBmean )
CSBvar <- cbm1 %*% ( CSBsumsq - diag(colnB) %*% CSBmean * CSBmean )
if (axisvariances>1) {
for (a1 in 1:(axisvariances-1)) {
  for (a2 in (a1+1):axisvariances) { 
  RSBcov[ ,a1,a2] <- rbm1 %*% ( RSBsumcp[ ,a1,a2] - diag(rownB) %*% RSBmean[ ,a1] * RSBmean[ ,a2] )
  CSBcov[ ,a1,a2] <- cbm1 %*% ( CSBsumcp[ ,a1,a2] - diag(colnB) %*% CSBmean[ ,a1] * CSBmean[ ,a2] )
} } } 

# OTHER CALCULATIONS

Fmat <- S@Rprofile %*% S@Rweights %*% S@Raxes
Gmat <- S@Cprofile %*% S@Cweights %*% S@Caxes

dmum1 <- diag( 1/(S@mu + (S@mu==0)) * (1-(S@mu==0)) )
Gbi <- Gmat %*% dmum1
Fbi <- Fmat %*% dmum1

# Calc inertia sum

inertia <- S@mu*S@mu
dmum2 <- diag( 1/(inertia + (S@mu==0)) * (1-(S@mu==0)) )
inertiasum <- sum(inertia)
inertiapc <- 100*inertia/inertiasum
cuminertiapc <- cumsum(inertiapc)
inertiapc <- round(100*inertiapc)/100
cuminertiapc <- round(100*cuminertiapc)/100
inertias <- cbind(inertia,inertiapc,cuminertiapc)

# Calc contributions and correlations

Xstd <- X/sum(X)
dr <- diag( as.vector(rowSums(Xstd)) )
dc <- diag( as.vector(colSums(Xstd)) )

Fsq <- Fmat * Fmat
RowCTR <- dr %*% Fsq %*% dmum2
Frs <- diag(1/rowSums(Fsq))
RowREP <- Frs %*% Fsq

Gsq <- Gmat * Gmat
ColCTR <- dc %*% Gsq %*% dmum2
Grs <- diag(1/rowSums(Gsq))
ColREP <- Grs %*% Gsq

bootca <- new("cabootcrsresults", br=S, 
  DataMatrix=X, rows=rows, columns=cols, 
  rowlabels=rowlabels, collabels=collabels,
  Rowprinccoord=Fmat, Colprinccoord=Gmat, Rowstdcoord=Fbi, Colstdcoord=Gbi,
  RowCTR=RowCTR, RowREP=RowREP, ColCTR=ColCTR, ColREP=ColREP, 
  RowVar=RSBvar, RowCov=RSBcov, ColVar=CSBvar, ColCov=CSBcov,
  inertiasum=inertiasum, inertias=inertias, 
  nboots=nboots, resampledistn=resampledistn, 
  multinomialtype=multinomialtype, sameaxisorder=sameaxisorder,
  poissonzeronewmean=poissonzeronewmean, newzeroreset=newzeroreset, 
  printdims=printdims, axisvariances=axisvariances )

summaryca(bootca, datasetname=as.character(datafile))

plotca(bootca, datasetname=as.character(datafile), 
               groupings=groupings, grouplabels=grouplabels, 
               plotsymbolscolours=plotsymbolscolours, 
               othersmonochrome=othersmonochrome, crpercent=crpercent)

bootca

}


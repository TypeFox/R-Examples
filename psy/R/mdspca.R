mdspca <-
function(datafile, supvar="no", supsubj="no", namesupvar=colnames(supvar,do.NULL=FALSE), namesupsubj=colnames(supsubj,do.NULL=FALSE), dimx=1, dimy=2, cx=0.75) 

{


#***********************************************
# missing values are omitted, normalization of var and supvar
#***********************************************

if ((is.na(supvar) || (supvar!="no")) && (!is.na(supsubj) && (supsubj=="no")))
{
svar <- 1
ssubj <- 0
p <- dim(datafile)[2]
supvar <- as.matrix(supvar)
pp <- dim(supvar)[2]
interm <- cbind(datafile,supvar)
interm <- na.omit(interm)
mat <- as.matrix(interm[,1:p])
matp <- as.matrix(interm[,(p+1):(p+pp)])
n <- dim(mat)[1]
for(i in 1:p) mat[,i] <- (mat[,i]-mean(mat[,i]))/(sqrt(n-1)*sd(mat[,i]))
for(i in 1:pp) matp[,i] <- (matp[,i]-mean(matp[,i]))/(sqrt(n-1)*sd(matp[,i]))
}

if ((is.na(supvar) || (supvar!="no")) && (is.na(supsubj) || (supsubj!="no")))
{
svar <- 1
ssubj <- 1
p <- dim(datafile)[2]
supvar <- as.matrix(supvar)
pp <- dim(supvar)[2]
supsubj <- as.matrix(supsubj)
ppp <- dim(supsubj)[2]
interm <- cbind(datafile,supvar,supsubj)
interm <- na.omit(interm)
mat <- as.matrix(interm[,1:p])
matp <- as.matrix(interm[,(p+1):(p+pp)])
supsubj <- as.matrix(interm[,(p+pp+1):(p+pp+ppp)])
n <- dim(mat)[1]
for(i in 1:p) mat[,i] <- (mat[,i]-mean(mat[,i]))/(sqrt(n-1)*sd(mat[,i]))
for(i in 1:pp) matp[,i] <- (matp[,i]-mean(matp[,i]))/(sqrt(n-1)*sd(matp[,i]))
}

if ((!is.na(supvar) && (supvar=="no")) && (is.na(supsubj) || (supsubj!="no")))
{
svar <- 0
ssubj <- 1
p <- dim(datafile)[2]
supsubj <- as.matrix(supsubj)
ppp <- dim(supsubj)[2]
interm <- cbind(datafile,supsubj)
interm <- na.omit(interm)
mat <- as.matrix(interm[,1:p])
supsubj <- as.matrix(interm[,(p+1):(p+ppp)])
n <- dim(mat)[1]
for(i in 1:p) mat[,i] <- (mat[,i]-mean(mat[,i]))/(sqrt(n-1)*sd(mat[,i]))
}

if ((!is.na(supvar) && (supvar=="no")) && (!is.na(supsubj) && (supsubj=="no")))
{
svar <- 0
ssubj <- 0
p <- dim(datafile)[2]
mat <- as.matrix(na.omit(datafile))
n <- dim(mat)[1]
for(i in 1:p) mat[,i] <- (mat[,i]-mean(mat[,i]))/(sqrt(n-1)*sd(mat[,i]))
}

#**********************************************
# definitions
#**********************************************

n <- dim(mat)[1]

names <- matrix(nrow=p+2)
one2 <- matrix(1, nrow=p)
load <- matrix(nrow=p, ncol=p)
loady <- matrix(nrow=p+2, ncol=p)

#***********************************************
# correlations and loadings
#***********************************************

names[1:p] <- attributes(datafile)$names
un <- matrix(1, ncol=p, nrow=n)
one <- matrix(1, nrow=n)

matcorp <- cor(mat)
decomp <- eigen(matcorp, symmetric=TRUE)
eigenval <- decomp$values
eigenvect <- decomp$vectors
eigenval <- pmax(0.00001*one2, eigenval)
load <- eigenvect*sqrt(kronecker(one2,t(eigenval)))

loady[1:p,] <- load[1:p,1:p]

#***********************************************
# supplementary variables
#***********************************************

if(svar==1)
{
one3 <- matrix(1,nrow=pp)
namesp <- matrix(nrow=pp)
namesp[1:pp] <- namesupvar
loadp <- matrix(nrow=pp, ncol=p)
loadyp <- matrix(nrow=pp, ncol=p)
loadp <- (t(matp)%*%mat%*%eigenvect)*kronecker(one3,t(1/sqrt(eigenval)))

loadyp[1:pp,] <- loadp[1:pp,1:p]
}


#***********************************************
# supplementary subjects
#***********************************************


if(ssubj==1)
{

nn <- ppp

mod <- matrix(nrow=nn)
for(i in 1:nn)
{
	factsub <- as.factor(supsubj[,i])
	mod[i] <- nlevels(factsub)
}
nmod <- sum(mod)

names2 <- matrix(nrow=nmod)
mat2 <- matrix(nrow=nmod, ncol=p)
load2 <- matrix(nrow=nmod, ncol=p)
loady2 <- matrix(nrow=nmod, ncol=p)

compt <- 0
mat <- as.data.frame(mat)
for(i in 1:nn) for(j in 1:mod[i])
{
	compt <- compt+1
	factsub <- as.factor(supsubj[,i])
	names2[compt] <- paste(namesupsubj[i],levels(factsub)[j])
	mat2[compt,] <- sapply(split(mat,factsub)[[j]],mean)
}

load2 <- (mat2%*%eigenvect)*sqrt(n/p)
loady2[1:nmod,] <- load2[1:nmod,1:p]



}


#****************************************************************************
#****************************************************************************
# Drawing
#****************************************************************************
#****************************************************************************

#****************** two more points for a non truncated drawing ********

loady[p+1,] <- 1.5
loady[p+2,] <- -1.5	
names[p+1] <- "."
names[p+2] <- "."

#****************************************************************************
#********************************  plots ***********************************
#****************************************************************************


par(pty="s")
if (is.na(supsubj) || (supsubj!="no"))
{
par(mfrow=c(1,2))
par(oma=c(0,0,0,0))
par(mar=c(0,0,0,0))
dimmax <- max(abs(loady2[,dimx])+0.2,abs(loady2[,dimy])+0.2,1.2)
#*************** new axes (centered) ********
par(mar=rep(0,4))
plot(x=c(-1*dimmax+0.1,dimmax-0.1),y=c(0,0),type="l",axes=FALSE,frame.plot=FALSE,ann=FALSE,xlim=c(-1*dimmax,dimmax),ylim=c(-1*dimmax,dimmax),col="grey")
lines(x=c(0,0),y=c(-1*dimmax+0.1,dimmax-0.1),type="l",col="grey")

#****************** plot of correlations ***********
symbols(x=loady2[,dimx], y=loady2[,dimy], squares=rep(.03,length(loady2[,dimy])), inches=FALSE, add=TRUE,fg="blue",bg="blue")

#***************** name plot *****************
text(x=loady2[,dimx],y=loady2[,dimy]-0.05,labels=names2,cex=cx)

}


#*************** new axes (centered) ********
if (!is.na(supsubj) && (supsubj=="no")) {par(mar=rep(0,4))}
plot(x=c(-1.1,1.1),y=c(0,0),type="l",axes=FALSE,frame.plot=FALSE,ann=FALSE,xlim=c(-1.2,1.2),ylim=c(-1.2,1.2),col="grey")
lines(x=c(0,0),y=c(-1.1,1.1),type="l",col="grey")

#************** circle (r=1)****************
symbols(x=0, y=0, circles=1, inches=FALSE, add=TRUE, lwd=2)

#****************** plot of correlations ***********
symbols(x=loady[,dimx], y=loady[,dimy], circles=rep(.01*cx*2,length(loady[,dimy])), inches=FALSE, add=TRUE,fg="grey",bg="red")

#***************** name plot *****************
text(x=loady[,dimx],y=loady[,dimy]-0.05,labels=names,cex=cx)



if (is.na(supvar) || (supvar!="no"))
{
#****************** plot of correlations sup var ***********
symbols(x=loadyp[,dimx], y=loadyp[,dimy], circles=rep(.01*cx*2,length(loadyp[,dimy])), inches=FALSE, add=TRUE,fg="grey",bg="green")

#***************** name plot *****************
text(x=loadyp[,dimx],y=loadyp[,dimy]-0.05,labels=namesp,cex=cx)
}

#****************** name of factors ***********
pf1 <- floor(100*eigenval[dimx]/sum(eigenval))
pf2 <- floor(100*eigenval[dimy]/sum(eigenval))
annotate1 <- paste("x = F",dimx," : ",pf1,"% var", sep="")
annotate2 <- paste("y = F",dimy," : ",pf2,"% var", sep="")
text(x=1,y=1,labels=annotate1,cex=cx)
text(x=1,y=1-0.05*cx*2,labels=annotate2,cex=cx)

par(mfrow=c(1,1))

#****************************************************************************
#******************************** end of plots ****************************
#****************************************************************************

}

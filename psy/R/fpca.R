fpca <-
function(formula=NULL, y=NULL, x=NULL,data, cx=0.75,pvalues="No", partial="Yes", input="data", contraction="No", sample.size=1)

#**********************************************
#
# datafile is the name of the dataframe that contains the data
# y is the number of the column related to the dependant variable (ex: y = 6)
# x is the vector of the number of the columns related to the independant variables
# (ex: x = c(1,2,3,4,5,7,8,9,10)
#
# option: pvalues is a vector of determined pvalues that will replace the correlations
# to the focus variable (if pvalues != "No" then partial = "No") (pvalues="No" by default)
# 
# q <- 1 (q>1 for a future version)
#
# option: partial is an option to present a focused PCA that is a simple renormalization
# of conventional PCA (partial="Yes" by default)
#
# option: input indicates wether the input correspond to data (default) or to a
# correlation matrix (if input != "data" then partial = "No") (input="data" by default)
#
# option: contraction change the appearance of the figure
# (if contraction="Yes" then pvalues="No") (contraction="No" by default)
#
# option: sample.size, size of the sample when input!="data"
#
#**********************************************

{


if (is.null(y))
    {
    call <- match.call()

        mf <- match.call(expand.dots = FALSE)
        m <- match(c("formula", "data", "subset"), names(mf), 0L)
    mf <- mf[c(1L, m)]
    mf$drop.unused.levels <- TRUE
    mf[[1L]] <- as.name("model.frame")
    mf <- eval(mf, parent.frame())


    Y <- model.extract(mf, "response")
    z=dim(mf)[2]
    x <- mf[,2:z]
    namey=names(mf)[1]
    namex=names(mf[,2:z])

    }

 datafile<-data
 
if (pvalues[1]!="No") partial <- "No"

#**********************************************
# definitions
#**********************************************

	p <- ifelse(is.null(y),ncol(x),length(x))
	if (input=="data") n <-ifelse(is.null(y),nrow(x),dim(datafile)[1]) else n <- sample.size
	if (input=="data") {if(is.null(y)) mat<-mf else mat<-matrix(ncol=p+1, nrow=n)  }
	names <- matrix(nrow=p+2)
	one2 <- matrix(1, nrow=p)
	load <- matrix(nrow=p, ncol=p)
	norm <- matrix(nrow=p, ncol=p-1)
	loadx <- matrix(nrow=p+2, ncol=p-1)
	loadyp <- matrix(nrow=p+2, ncol=p-1)
	loadym <- matrix(nrow=p+2, ncol=p-1)
	loady <- matrix(nrow=p+2, ncol=p-1)
	q <- 1

if (input=="data")
{

	#***********************************************
	# missing values are NOT omitted
	# input x (independant) and y (dependant variable)
	# x and y are normalized
	# correlations of (x,y)
	#***********************************************
	
if (is.null(formula)){
	q <- min(q,p-1)
	mat[,1] <- datafile[,c(y)]
	
  namey<-ifelse(is.numeric(y),attributes(datafile)$names[y],y)

	for(i in 1:p)	
		{
		mat[,i+1] <- datafile[,c(x[i])]
		names[i] <- ifelse(is.numeric(x[i]),attributes(datafile)$names[x[i]],x[i])
		}


	mat <- na.omit(mat) #following command used to work for data.frames containing NA, didn't work if no NAs
	}
	
	n <- dim(mat)[1]
	xv <- matrix(ncol=p, nrow=n)
	yv <- matrix(nrow=n)
	un <- matrix(1, ncol=p, nrow=n)
	one <- matrix(1, nrow=n)

	for(i in 1:p)	
		{
		xv[,i] <- (mat[,i+1]-mean(mat[,i+1]))/(sqrt(var(mat[,i+1]))*sqrt(n-1))
		}
	yv <- mat[,1]
	yv <- (yv-mean(yv))/(sqrt(var(yv))*sqrt(n-1))

	matcor <- cor(mat)
} else

{
  if (is.null(formula)){

	namey<-ifelse(is.numeric(y),attributes(datafile)$names[y],y)

		for(i in 1:p)
		{
	names[i] <- ifelse(is.numeric(x[i]),attributes(datafile)$names[x[i]],x[i])
	}
		matcor <- matrix(nrow=p+1,ncol=1)
	matcor[1,1] <- 1
	matcor[2:(p+1),1] <- datafile[x,y]
	matcorp <- datafile[x,x]
	}
	
	else {
	
	namex=names(mf[names(mf)!=namey])
  X<-as.data.frame(t(x))[,c(namex)]
  matcor<- t(as.data.frame(t(Y))[,c(namey,namex)])
  matcorp <- X
  
	}
	
	
	decomp <- eigen(matcorp, symmetric=TRUE)
	eigenval <- decomp$values
	eigenvect <- decomp$vectors
	eigenval <- pmax(0*one2, eigenval)
	load <- eigenvect*sqrt(kronecker(one2,t(eigenval)))
}

#***********************************************
#traditional FPCA
#***********************************************

if (pvalues[1]=="No")
	{
	if (input=="data")
		{
		if (partial=="Yes")
			{
			#***********************************************
			# xp is x related to y (projection on y orthog)
			#***********************************************	
			scal <- t(t(xv)%*%yv)
			xp <- xv - (un*yv)*kronecker(one,scal)
			}else

			{
			xp <- xv	
			}
		#***********************************************
		# decomposition of the covariance matrix of xp
		#***********************************************
		
		#matcorp <- var(xp, na.method="omit")
                matcorp <- var(xp)
		decomp <- eigen(matcorp, symmetric=TRUE)
		eigenval <- decomp$values
		eigenvect <- decomp$vectors
		eigenval <- pmax(0*one2, eigenval)
		load <- eigenvect*sqrt(kronecker(one2,t(eigenval)))
		}
	
	#***********************************************
	# renormalization of the loadings
	#***********************************************
	
	if (contraction=="No")
		{
		for(i in 1:p)
			{
			for(j in 1:p-1)
				{
				norm[i,j] <- sqrt(load[i,1]*load[i,1] + load[i,j+1]*load[i,j+1])
				loadx[i,j] <- load[i,1]*sqrt(2-2*abs(matcor[i+1,1]))/norm[i,j]
				if (matcor[i+1,1] > 0) loadyp[i,j] <- load[i,j+1]*sqrt(2-2*matcor[i+1,1])/norm[i,j] else loadym[i,j] <- load[i,j+1]*sqrt(2+2*matcor[i+1,1])/norm[i,j]
				loady[i,j] <- load[i,j+1]*sqrt(2-2*abs(matcor[i+1,1]))/norm[i,j] - 0.05
				}
			}
		}

	if (contraction!="No")
		{
		for(i in 1:p)
			{
			for(j in 1:p-1)
				{
				norm[i,j] <- sqrt(load[i,1]*load[i,1] + load[i,j+1]*load[i,j+1])
				loadx[i,j] <- 1.5*load[i,1]*(1-abs(matcor[i+1,1]))/norm[i,j]
				if (matcor[i+1,1] > 0) loadyp[i,j] <- 1.5*load[i,j+1]*(1-matcor[i+1,1])/norm[i,j] else loadym[i,j] <- 1.5*load[i,j+1]*(1+matcor[i+1,1])/norm[i,j]
				loady[i,j] <- 1.5*load[i,j+1]*(1-abs(matcor[i+1,1]))/norm[i,j] - 0.05
				}
			}
		}
		

	}
#***********************************************
#pvalued FPCA
#***********************************************

if (pvalues[1]!="No")

	{
	matcorp <- var(xv, na.rm=TRUE)
	decomp <- eigen(matcorp, symmetric=TRUE)
	eigenval <- decomp$values
	eigenvect <- decomp$vectors
	eigenval <- pmax(0*one2, eigenval)
	load <- eigenvect*sqrt(kronecker(one2,t(eigenval)))
	
	#***********************************************
	# renormalization of the loadings (1.5 is for drawing convenience)
	#***********************************************
	
	pnorm <- matrix(nrow=p)
	pvaluesabs <- abs(pvalues)
	pvaluesabs <- pmax(pvaluesabs,0.001)
	for (i in 1:p) if (pvaluesabs[i]==0) pnorm[i] <- 0 else pnorm[i] <- pvaluesabs[i]^(log(pvaluesabs[i])/-50)
	
	
	for(i in 1:p)
		{
		for(j in 1:p-1)
			{
			norm[i,j] <- sqrt(load[i,1]*load[i,1] + load[i,j+1]*load[i,j+1])
			loadx[i,j] <- 1.5*load[i,1]*pnorm[i]/norm[i,j]
			if (pvalues[i] > 0) loadyp[i,j] <- 1.5*load[i,j+1]*pnorm[i]/norm[i,j] else loadym[i,j] <- 1.5*load[i,j+1]*pnorm[i]/norm[i,j]
			loady[i,j] <- 1.5*load[i,j+1]*pnorm[i]/norm[i,j] - 0.05
			}
		}

	}


#****************************************************************************
#****************************************************************************
# Drawing
#****************************************************************************
#****************************************************************************

#****************** two more points for a non truncated drawing ********

if (is.null(formula)){
for(j in 1:p-1)
{
loadx[p+1,j] <- 1.5
loady[p+1,j] <- 1.5
loadx[p+2,j] <- -1.5
loady[p+2,j] <- -1.5 
}
names[p+1] <- "."
names[p+2] <- "."
}

else {names=namex}
#****************************************************************************
#******************************** q plots ***********************************
#****************************************************************************

j <- 1

{

#*************** new axes (centered) ********
par(pty="s")
par(mar=rep(0,4))
plot(x=c(-1.6,1.6),y=c(0,0),type="l",axes=FALSE,frame.plot=FALSE,ann=FALSE,xlim=c(-1.7,1.7),ylim=c(-1.7,1.7),col="grey")
lines(x=c(0,0),y=c(-1.6,1.6),type="l",col="grey")

#********************************************************
#traditionnal FPCA
#********************************************************

if (pvalues[1]=="No")
{

#************** circles (r=0, r=0.2, ...)****************

radius <- matrix(nrow=5)

if (contraction=="No")
{
radius[1] <- 1.414
radius[2] <- 1.265
radius[3] <- 1.095
radius[4] <- 0.894
radius[5] <- 0.632
} 
else
{
radius[1] <- 1.5
radius[2] <- 1.2
radius[3] <- 0.9
radius[4] <- 0.6
radius[5] <- 0.3
}

symbols(x=0, y=0, circles=radius[1], inches=FALSE, add=TRUE, lwd=2)
symbols(x=c(0,0,0,0), y=c(0,0,0,0), circles=radius[2:5], inches=FALSE, add=TRUE, lwd=1,fg="grey")

}

#********************************************************
#pvalued FPCA
#********************************************************

if (pvalues[1]!="No")
{
#************** circles (p=0.1, p=0.05, ...)****************
symbols(x=0, y=0, circles=1.5, inches=FALSE, add=TRUE, lwd=2)
symbols(x=c(0,0,0), y=c(0,0,0), circles=c(1.35,1,.557), inches=FALSE, add=TRUE, lwd=1,fg="grey")
symbols(x=0, y=0, circles=1.254, inches=FALSE, add=TRUE, lwd=1,fg="red")
}

#************** dependant variable *****************
symbols(x=0, y=0, circles=0.03, bg="black", inches=FALSE, add=TRUE, lwd=1)

#********************************************************
#traditionnal FPCA
#********************************************************

if (pvalues[1]=="No")
{
#************** circle with p = 5% *****************
if (input=="data")
{
	e <- exp(1.96*2/sqrt(n-3))
}else 
{
	if (sample.size<4) e <- 0 else e <- exp(1.96*2/sqrt(sample.size-3))
}
if (contraction=="No") rayonsign <- sqrt(2-2*(e-1)/(1+e)) else rayonsign <- (1-(e-1)/(1+e))*1.5
symbols(x=0, y=0, circles=rayonsign, inches=FALSE, add=TRUE, lwd=1, fg="red")

#**************** legends : r=0, r=0.2, ... ****************
text(x=c(rep(0.01,5)),y=radius+.04,
     labels=c("r = 0","r = 0.2","r = 0.4","r = 0.6","r = 0.8"),cex=0.5) 
}

#********************************************************
#pvalued FPCA
#********************************************************

if (pvalues[1]!="No")
{
#**************** legends : p=0, p=0.1, ... ****************
text(x=c(rep(0.01,5)),y=c(.563,.982,1.239,1.335,1.48),
     labels=c("p < 0.001","p = 0.01","p = 0.05","p = 0.1","p = 1"),cex=cx)
}

#****************** plot of positive correlations ***********
symbols(x=loadx[,j], y=loadyp[,j], circles=rep(.03,length(loadyp[,j])), inches=FALSE, add=TRUE,fg="blue",bg="green")

#****************** plot of negative correlations ***********
symbols(x=loadx[,j], y=loadym[,j], circles=rep(.03,length(loadym[,j])), inches=FALSE, add=TRUE,fg="red",bg="yellow")

#***************** names ******************
text(x=-0.18,y=-.12,labels=namey, cex=cx+0.25)#focus variable 
text(x=loadx[,j],y=loady[,j],labels=names,cex=cx)#other variables

#****************** name of factors ***********
#annotate <- paste("Factors : 1,", j, sep="")
#text(x=1,y=1.3,labels=annotate,cex=cx)

#****************************************************************************
#******************************** end of q plots ****************************
#****************************************************************************

}

#******************* end **********************
}

library('geostatsp')
matrix(NNmat(7, 7)[,25], 7, 7)

myr = raster(extent(0,600,0,300), nrows=25,ncols=30)
theNN = NNmat(myr)


params=c(range = 6*xres(myr),
		cellSize=xres(myr),
		shape=2,
		variance=1600)


# precision matrix without adjusting for edge effects
precMat =maternGmrfPrec(theNN, param=params,
		adjustEdges=FALSE,
		adjustParam=FALSE) 


# find better parameters using a numerical optimizer
precMatAdj =maternGmrfPrec(theNN, param=params, 
		adjustEdges=TRUE,adjustParam=TRUE) 

# and with the adjustment
precMatCorr =maternGmrfPrec(theNN, param=params, 
		adjustEdges=TRUE,adjustParam=FALSE) 


Nx = attributes(theNN)$Nx
Ny = attributes(theNN)$Ny
midcell = Nx*round(Ny/2) + round(Nx/2) # the middle cell
edgecell = Nx*5 + 5 # cell near corner


# show precision of mid cell
precMid=matrix(precMat[,midcell], Ny, Nx, byrow=TRUE)
precMid[round(Ny/2)+seq(-3, 3), round(Nx/2)+seq(-3, +3)]

if(Sys.info()['user'] =='patrick') {

# invert to get variance matrices
precMatAdjInv = solve(precMatAdj)
precMatCorrInv = solve(precMatCorr)
precMatInv = solve(precMat)

# check range of diagonal entries, should be
params['variance']
range(diag(precMatInv))
range(diag(precMatCorrInv))
range(diag(precMatAdjInv))

# map marginal variances
tempA=tempC=tempU = attributes(precMat)$raster
values(tempC) = diag(precMatCorrInv)
values(tempU) = diag(precMatInv)
values(tempA) = diag(precMatAdjInv)

if(!interactive()){
	pdf("maternGmrfMarginalVarRasters.pdf",height=10,width=5)
}
par(mfrow=c(3,1))
plot(tempU)
plot(tempC)
plot(tempA)
if(!interactive()){
	dev.off()
}



# variance matrices
Ncell = Nx*Ny
midVec = sparseMatrix(midcell,1,x=1,dims=c(Ncell,1))
edgeVec = sparseMatrix(edgecell,1,x=1,dims=c(Ncell,1))

therast = attributes(precMat)$raster
values(therast)=NA
varRast = brick(therast, therast,therast,therast,therast,therast)
names(varRast) = c('mid','edge','midCor','edgeCor','midAdj','edgeAdj')
values(varRast[['mid']]) = as.vector(Matrix::solve(precMat, midVec))
values(varRast[['edge']]) = as.vector(Matrix::solve(precMat, edgeVec))
values(varRast[['midCor']]) =  as.vector(Matrix::solve(precMatCorr, midVec))
values(varRast[['edgeCor']]) = as.vector(Matrix::solve(precMatCorr, edgeVec))
values(varRast[['midAdj']]) =  as.vector(Matrix::solve(precMatAdj, midVec))
values(varRast[['edgeAdj']]) = as.vector(Matrix::solve(precMatAdj, edgeVec))

projection(varRast) = "+init:units=m"


if(!interactive()){
	pdf("maternGmrfPredRasters.pdf",height=11,width=4)
}
par(mfrow=c(7,2))

mMid = matern(varRast[[1]],y=xyFromCell(varRast,midcell),param=params)
mEdge = matern(varRast[[1]],y=xyFromCell(varRast,edgecell),param=params)
plot(mMid)
plot(mEdge)


plot(varRast[['mid']])
plot(varRast[['edge']])

plot(varRast[['mid']]-mMid)
plot(varRast[['edge']]-mEdge)

plot(varRast[['midCor']])
plot(varRast[['edgeCor']])
plot(varRast[['midCor']]-mMid)
plot(varRast[['edgeCor']]-mEdge)

plot(varRast[['midAdj']])
plot(varRast[['edgeAdj']])
plot(varRast[['midAdj']]-mMid)
plot(varRast[['edgeAdj']]-mEdge)


if(!interactive()){
dev.off()
}


# compare covariance matrix to the matern

xseqFull = seq(-1.5, 1.5, len=401)*params["range"]

diagdist = sqrt(2)

xymid = xyFromCell(varRast, midcell)

xseq = seq(xmin(varRast)+xres(varRast)/2,xmax(varRast), by=xres(varRast))
yseq = seq(ymin(varRast)+yres(varRast)/2,ymax(varRast), by=yres(varRast))

vx= extract(varRast, 
		SpatialPoints(cbind(xseq,
						rep(xymid[2],length(xseq)))
		))
vy = extract(varRast, 
		SpatialPoints(cbind(
						rep(xymid[1],length(yseq)),
						yseq)
		)
)

aseq = seq(0,100,by=xres(varRast))
aseq = sort(unique(c(aseq, -aseq)))

lur = SpatialPoints(cbind(xymid[1]+aseq,xymid[2]+aseq))
lll = SpatialPoints(cbind(xymid[1]-aseq,xymid[2]+aseq))

vur = extract(varRast, 
		lur		)


vll= extract(varRast, 
		lll		)
xyedg = xyFromCell(varRast, edgecell)

xseqe = seq(xmin(varRast)+xres(varRast)/2,xmax(varRast), by=xres(varRast))
yseqe = seq(ymin(varRast)+yres(varRast)/2,ymax(varRast), by=yres(varRast))

vxe= extract(varRast, 
		SpatialPoints(cbind(xseqe,
						rep(xyedg[2],length(xseqe)))
		))
vye = extract(varRast, 
		SpatialPoints(cbind(
						rep(xyedg[1],length(yseqe)),
						yseqe)
		)
)

aseqe = seq(0,100,by=xres(varRast))
aseqe = sort(unique(c(aseqe, -aseqe)))

lure = SpatialPoints(cbind(xyedg[1]+aseqe,xyedg[2]+aseqe))
llle = SpatialPoints(cbind(xyedg[1]-aseqe,xyedg[2]+aseqe))

vure = extract(varRast, 
		lure		)

vlle= extract(varRast, 
		llle		)




theedg = c('edge','edgeAdj','edgeCor')



themid = c('mid','midAdj','midCor')


thecol =c(mid='red',edge='red',midCor='green',edgeCor='green',
		midAdj='blue',edgeAdj='blue')
thelwd = c(mid=3,midCor=1,midAdj=1,edge=3,edgeCor=1,edgeAdj=1)

if(!interactive()){
	pdf("maternGmrfPred.pdf",height=5,width=10)
}
par(mfrow=c(2,2),mar=c(3,2,0,0))
# middle cell, full plot
plot(xseqFull, matern(xseqFull, param=params),
			type = 'l',ylab='cov', xlab='dist',
			main="matern v gmrf",
			ylim=c(0,params['variance']*1.1),col='grey',lwd=9)

matlines(xseq-xymid[1],vx[,themid],col=thecol[themid],lty=1,
		type='o',lwd=thelwd[themid])
matlines(yseq-xymid[2],vy[,themid],col=thecol[themid],lty=1,
		type='o',lwd=thelwd[themid])


matlines(aseq*diagdist, 
		vur[,themid],col=thecol[themid],lty=2,type='o',
		lwd=thelwd[themid])
matlines(aseq*diagdist, 
		vll[,themid],col=thecol[themid],lty=2,type='o',lwd=thelwd[themid])






	legend("topright", lty=1, 
			col=thecol[c('edge','edgeCor','edgeAdj')],
			legend=c('vanilla','edge','adj')
	)

	
# edge cell, full plot
	plot(xseqFull, matern(xseqFull, param=params),
			type = 'l',ylab='cov', xlab='dist',
			main="matern v gmrf",
			ylim=c(0,params['variance']*1.1))

	
	
	matlines(xseqe-xyedg[1],vxe[,theedg],col=thecol[theedg],lty=1,
			lwd=thelwd[theedg])
	matlines(yseqe-xyedg[2],vye[,theedg],col=thecol[theedg],lty=1,
			lwd=thelwd[theedg])
	
	
	matlines(aseqe*diagdist, vure[,theedg],col=thecol[theedg],lty=2,
			lwd=thelwd[themid])
	
	matlines(aseqe*diagdist, vlle[,theedg],col=thecol[theedg],lty=2,
			lwd=thelwd[theedg])
	
	# middle cell detail plot
	
plot(xseqFull, matern(xseqFull, param=params),
		type = 'l',ylab='cov', xlab='dist',
		main="matern v gmrf",
		ylim=params['variance']*c(0.5,1.05),
		xlim=c(-0.5, 0.5)*params['range'],col='grey',lwd=4)



	matlines(xseq-xymid[1],vx[,themid],col=thecol[themid],
			lty=1,
			lwd=thelwd[themid])
	matlines(yseq-xymid[2],vy[,themid],col=thecol[themid],lty=1,
			lwd=thelwd[themid])
	

	matlines(aseq*diagdist, vur[,themid],col=thecol[themid],lty=2,
			lwd=thelwd[themid])
	
	matlines(aseq*diagdist, vll[,themid],col=thecol[themid],lty=2,
			lwd=thelwd[themid])
	
	
	
	# edge cells
	
plot(xseqFull, matern(xseqFull, param=params),
		type = 'l',ylab='cov', xlab='dist',
		main="matern v gmrf",
		ylim=params['variance']*c(0.5,1.05),
		xlim=c(-0.5, 0.5)*params['range'])


matlines(xseqe-xyedg[1],vxe[,theedg],col=thecol[theedg],
		lty=1,
		lwd=thelwd[theedg])
matlines(yseqe-xyedg[2],vye[,theedg],col=thecol[theedg],
		lty=1,
		lwd=thelwd[theedg])


matlines(aseqe*diagdist, vure[,theedg],col=thecol[theedg],
		lty=2,
		lwd=thelwd[themid])

matlines(aseqe*diagdist, vlle[,theedg],col=thecol[theedg],
		lty=2,
		lwd=thelwd[theedg])


	
	if(!interactive()){
		dev.off()	
	}

# matern variance matrix
	covMatMatern = 
			matern(attributes(precMat)$raster, 
					param=params[c('range','shape','variance')])
	
	covMatMaternAdj = 
			matern(attributes(precMat)$raster, 
					param=attributes(precMatAdj)$param$optimal[
							c('range','shape','variance')])

	
	
	prodUncor = covMatMatern %*% (precMat)
	prodCor = covMatMatern %*% (precMatCorr)
	prodAdj = covMatMatern %*% (precMatAdj)
	prodAdjOpt = covMatMaternAdj %*% (precMatAdj)
	
	
	
	thebreaks = c(-100, -0.25, 0.25, 0.75, 1.5, 5, 10, 100)

	thecol = terrain.colors(length(thebreaks)-1)
	
	if(!interactive()){
		pdf("maternGmrfDiagRast.pdf",height=8,width=5)
	}
	par(mfrow=c(3,1))
	temp = attributes(precMat)$raster
	values(temp) = diag(prodCor)
	plot(temp,breaks=thebreaks,col=thecol,legend=F)
#	mapmisc::legendBreaks('right',breaks=thebreaks, col=thecol)

	temp = attributes(precMat)$raster
	values(temp) = diag(prodUncor)
	plot(temp,breaks=thebreaks,col=thecol,legend=FALSE)

	temp = attributes(precMat)$raster
	values(temp) = diag(prodAdj)
	plot(temp,breaks=thebreaks,col=thecol,legend=FALSE)
	
	if(!interactive()){
		dev.off()
	}
	
	if(!interactive()){
		pdf("maternGmrfHist.pdf",height=9,width=7)
	}
	par(mfrow=c(4,2),mar=c(4,4,0,0))	
	
	hist(Matrix::diag(prodUncor),breaks=60,ylab='vanilla')
	abline(v=1)
	hist(prodUncor[lower.tri(prodUncor,diag=FALSE)],breaks=60)	
	abline(v=0)
	
	thebreaks = c(-10,seq(0.4, 1.1, len=61),10)
	
	if(interactive()){
	hist(Matrix::diag(prodCor),breaks=thebreaks,ylab='edge')
	} else {
		hist(Matrix::diag(prodCor),breaks=60,ylab='edge')
	}
	abline(v=1)
	hist(prodCor[lower.tri(prodCor,diag=FALSE)],breaks=60)	
	abline(v=0)

	if(interactive()){
	hist(Matrix::diag(prodAdj),breaks=thebreaks,ylab='adj')
} else {
	hist(Matrix::diag(prodAdj),breaks=60,ylab='edge')
}
abline(v=1)
	hist(prodAdj[lower.tri(prodAdj,diag=FALSE)],breaks=60)	
	abline(v=0)

	if(interactive()){
	hist(Matrix::diag(prodAdjOpt),breaks=thebreaks,ylab='adj optimal')
} else {
	hist(Matrix::diag(prodAdjOpt),breaks=60,ylab='edge')
}
abline(v=1)
	hist(prodAdjOpt[lower.tri(prodAdjOpt,diag=FALSE)],breaks=60)	
	abline(v=0)
	
	
	if(!interactive()){
		dev.off()	
	}
	
	
	# testing marginal variance
	if(FALSE){
	
# order one
	Sa = seq(from=4.01,to=5,len=10)
	maxvar=NULL
	for(a in Sa) {
		myp = theNN
		myp@x = c(4+a^2,-2*a,2,1,0,0,0,0)[myp@x]
		myv = solve(myp)
		values(myr) = diag(myv)
		maxvar = c(maxvar,max(myr[15,]))
		
	}
	
	phi = 4/Sa
	scale = sqrt(Sa-4)
	range = sqrt(8)/scale
	
	margvar = (4/pi)^(scale^0.5/2)/(4*pi*(Sa-4))
	margvar = (4/pi)^(scale^2/2)/(4*pi*(Sa-4))
	
	margvarOld = 1/(4*pi*(Sa-4))
	
	
	par(mfrow=c(4,1),mar=c(2,4,0,0))
	they= c(-0.1, 0.1)
	plot(range, maxvar/margvar-1,log='x',ylab='range',ylim=they)
	lines(range, maxvar/margvarOld-1,col='grey')
	abline(h=0)
	
	plot(1-phi, maxvar/margvar-1,ylab='phi',ylim=they)
	lines(1-phi, maxvar/margvarOld-1,col='grey')
	abline(h=0)
	plot(scale, maxvar/margvar-1,ylab='scale',ylim=they)
	lines(scale, maxvar/margvarOld-1,col='grey')
	abline(h=0)
	plot(Sa-4, maxvar/margvar-1,log='x',ylab='Sa-4',ylim=they)
	lines(Sa-4, maxvar/margvarOld-1,col='grey')
	abline(h=0)
	max(maxvar/margvar)
	them=which.max(maxvar/margvar)
	c(range[them],phi[them],scale[them],Sa[them])
	
	
# order two
	
	Sa2 = seq(from=4.01,to=5,len=40)
	maxvar2=NULL
	for(a in Sa2) {
		myp = theNN
		myp@x = c("1" = a*(a*a+12),
				"2" = -3*(a*a+3),
				"3" = 6*a,
				"4" = 3*a, 
				"5" =  -3,
				"6" = -1)[myp@x]
		myv = solve(myp)
		values(myr) = diag(myv)
		maxvar2 = c(maxvar2,max(myr[15,]))
		
	}
	
	phi2 = 4/Sa2
	scale2 = sqrt(Sa2-4)
	range2 = sqrt(8*2)/scale2
	

	margvar = (4/pi)^(scale^2/2)/(4*pi*2*(Sa2-4)^2)
	margvarOld = 1/(4*pi*2*(Sa2-4)^2)
	

	par(mfrow=c(4,1),mar=c(2,4,0,0))
	they= c(-0.1, 0.1)
	plot(range, maxvar2/margvar-1,log='x',ylab='range',ylim=they)
	lines(range, maxvar2/margvarOld-1,col='grey')
	abline(h=0)
	
	plot(1-phi, maxvar2/margvar-1,ylab='phi',ylim=they)
	lines(1-phi, maxvar2/margvarOld-1,col='grey')
	abline(h=0)
	plot(scale, maxvar2/margvar-1,ylab='scale',ylim=they)
	lines(scale, maxvar2/margvarOld-1,col='grey')
	abline(h=0)
	plot(Sa-4, maxvar2/margvar-1,log='x',ylab='Sa-4',ylim=they)
	lines(Sa-4, maxvar2/margvarOld-1,col='grey')
	abline(h=0)
	max(maxvar2/margvar)
	them=which.max(maxvar2/margvar)
	c(range2[them],phi2[them],scale2[them],Sa2[them])
	

	
plot(xseqFull, matern(xseqFull, param=params),
			type = 'l',ylab='cov', xlab='dist',
			main="matern v gmrf",
			ylim=params['variance']*c(0.25,1.05),
			xlim=c(-0.75, 0.75)*params['range'])
	
	# middle cell
	themid='mid'
	
	plot(xseq-xymid[1],
			vx[,themid] - matern(xseq-xymid[1], param=params),
			col='red',lty=thelty[themid],
			lwd=thelwd[themid],type='o',pch=1,ylim=c(-55,12))
	matlines(yseq-xymid[2],
			vy[,themid]- matern(yseq-xymid[2], param=params),
			col='orange',lty=thelty[themid],
			lwd=thelwd[themid],type='o',pch=2)
	
	
	matlines(aseq*diagdist, 
			vur[,themid]- matern(aseq*diagdist, param=params),
			col='yellow',lty=thelty[themid],
			lwd=thelwd[themid],type='o',pch=3)
	
	matlines(aseq*diagdist, vll[,themid]- matern(aseq*diagdist, param=params),
			col='pink',lty=thelty[themid],
			lwd=thelwd[themid],type='o',pch=4)

xseqFull = seq(-400,400,len=1000)

	mbase =  matern(xseqFull, param=
					c(params[c('variance','range')],
							params['shape'])
	)
	lines(xseqFull, 
			matern(xseqFull, 
					param=
							c(params['variance'],
									params['range']/1.025,
									params['shape']/1.23)
			)-mbase,
			col='blue'
)
	lines(xseqFull, 
		matern(xseqFull, 
		param=
		c(params['variance'],
					params['range']/1.0225,
					params['shape']/1.22)
						)-mbase,
						col='grey'
)
}

# test providing ar parameter

paramsAR = c(shape=1,oneminusar=0.1)
N=theNN;param=paramsAR;adjustEdges=TRUE;adjustParam=FALSE 

precMatAR =maternGmrfPrec(theNN, param=paramsAR,
		adjustEdges=FALSE,
		adjustParam=FALSE) 


# find better parameters using a numerical optimizer
precMatAdjAR =maternGmrfPrec(theNN, param=paramsAR, 
		adjustEdges=TRUE,adjustParam=TRUE) 


# and with the adjustment
precMatCorrAR =maternGmrfPrec(theNN, param=paramsAR, 
		adjustEdges=TRUE,adjustParam=FALSE) 
}
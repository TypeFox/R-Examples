getPrec = function(shape,a)	{
	a = as.numeric(a)
	if(shape == 0){
		precEntries = c(
				"1" = a,
				"2" = -1,
				"3" = 0,
				"4" = 0, 
				"5" =  0,
				"6" = 0)
	} else if(shape==1) {	
		
		precEntries = c("1" = 4 + a^2,
				"2" = -2*a,
				"3" = 2,
				"4" = 1, 
				"5" =  0,
				"6" = 0)
	} else if(shape==2) {
		precEntries = c("1" = a*(a*a+12),
				"2" = -3*(a*a+3),
				"3" = 6*a,
				"4" = 3*a, 
				"5" =  -3,
				"6" = -1)
	} else {
		stop("shape parameter must be 0, 1 or 2")			
	}
	precEntries
}

maternGmrfPrec = function(N,...) {
	UseMethod("maternGmrfPrec")	
}

maternGmrfPrec.default = function(N, Ny=N,
		param=c(variance=1, range=1, shape=1, cellSize=1),
		 ...) {

	if(!'shape' %in% names(param)){
		warning("shape parameter appears to be missing")
	}
	
	theNNmat  = NNmat(N, Ny, nearest=param['shape']+1)

	maternGmrfPrec(theNNmat, param,...)
	
}


maternGmrfPrec.matrix = function(N, ...) {
	
	N = as(N, "dsCMatrix")
	maternGmrfPrec(N, ...)
	
}


maternGmrfPrec.dsCMatrix = function(N, 	
		param=c(variance=1, range=1, shape=1, cellSize=1),
		adjustEdges=FALSE,
		...) {

	param['nugget']=0
	
	names(param) = gsub("^var$", "variance", names(param))
	if(!any(names(param)=='variance')) {		
		if ('conditionalVariance' %in% names(param)){
			param['variance']=NA
		} else {
			param['variance']=1
		}
	} else {
		if ('conditionalVariance' %in% names(param)){
			if(any(names(param)=='oneminusar')) {
				warning(
						"both conditionalVariance and variance were supplied, ignoring variance"
				)	
				param['variance']=NA
			} else {
				warning(
						"both conditionalVariance and variance were supplied, ignoring conditionalVariance"
				)	
				param['conditionalVariance']=NA
			}
		}
	}


	
	if(is.null(attributes(N)$raster)) {
		if(!any(names(param)=='cellSize')){
			param['cellSize']=1
		}
		if(!all(c("Nx","Ny")%in% names(attributes(N)))) {
			Nx = Ny = sqrt(ncol(N))
			if(Nx!=round(Nx)){
				warning("N should have Nx and Ny attributes")
			} 
		} else {
				Nx = attributes(N)$Nx
				Ny = attributes(N)$Ny
		}
		theraster = list(
				nrows=Ny,
				ncols=Nx, 
				xmn=0,xmx=Nx*
						param['cellSize'],
				ymn=0, ymx=Ny*
						param['cellSize']
		)	
		if(any(installed.packages()[,'Package'] == 'raster')) {
		theraster = do.call(raster::raster,
				theraster)
		}
		
	} else {
		
		theraster = attributes(N)$raster
		Nx = attributes(theraster)$ncols
		Ny = attributes(theraster)$nrows
		param['cellSize'] = 
				( attributes(theraster)$extent@xmax - 
					attributes(theraster)$extent@xmin) /
				Nx
		
	}
	
	
	
	paramInfo = list(
		original=param,
		theo = NULL,
		optimal = NULL
)
	
	
	midcellCoord = c(round(Nx*0.4),round(Ny*0.4)) # the middle cell
	midcell = c(Nx*(Ny-midcellCoord[2]) + midcellCoord[1])
	midVec = sparseMatrix(midcell,1,x=1,
			dims=c(ncol(N),1))

	distVecFull =  expand.grid(x=seq(1, Nx),
			y=seq(Ny,1))
	buffer = param['shape']+2
	
	cellVec = seq(1, nrow(distVecFull)) 
	
	isInner =  
			distVecFull[,'x']> buffer & 
					distVecFull[,'y']>buffer &
							distVecFull[,'x']< (Nx-buffer) &
							distVecFull[,'y'] < (Ny-buffer) 

	distVecFull = t(distVecFull)-midcellCoord 
	distVecFull = sqrt(apply(distVecFull^2,2,sum))	
		
	
	
	# if 1 - AR parameter supplied
	if(all(c("oneminusar","shape") %in% names(param))){
		a = 4/(1-param['oneminusar'])	
		
		paramInfo$theo = c(param[c('shape','cellSize','oneminusar')],
				rangeInCells=as.numeric(sqrt(8*param['shape'])/sqrt(a-4)),
				a=as.vector(a), 
				scale = as.numeric(sqrt(a-4))
			)
			
			
		if ('conditionalVariance' %in% names(param)){
			if(param['shape'] != 0) {
				paramInfo$theo['variance'] =   
						param['conditionalVariance']/
						(4*pi*param['shape'] *(a-4)^(param['shape'] ))
			} else {
				paramInfo$theo['variance'] =  
						param['conditionalVariance']/(4*pi)
			}
			paramInfo$theo['conditionalVariance'] = 
					param['conditionalVariance']	
		} else {
			paramInfo$theo['variance'] = param['variance']
			if(param['shape'] != 0) {
				paramInfo$theo['conditionalVariance']  =   
					paramInfo$theo['variance']*
						(4*pi*param['shape'] *(a-4)^(param['shape'] ))
			} else {
				paramInfo$theo['conditionalVariance'] =   
					paramInfo$theo['variance']*(4*pi)
			}
			
	 
		}
			
			
		if(min(c(Nx,Ny)<3*paramInfo$theo['rangeInCells'])){
				warning("grid is ", Nx, " by ", Ny,
						"which may be too small for range ",
						paramInfo$theo['rangeInCells'], " cells, 1- Ar param",
						param['oneminusar'])
			}
			
		paramInfo$theo['range'] = as.numeric(
				paramInfo$theo['rangeInCells']*param['cellSize']
			)
			
	
	#####################	
	# else range supplied
	##########################
	} else if(all(c('range','shape') %in% names(param))){
		param['rangeInCells'] = as.numeric(param['range']/param['cellSize'])

		scale = as.numeric(
				sqrt(8*param['shape'])/param["rangeInCells"]
		)
		
		
		a = (scale^2 + 4) 
		
		paramInfo$theo = c(param, 	a=as.numeric(a),
				oneminusar = as.numeric(1-4/a))
 
		
		if(param['shape'] != 0) {
			paramInfo$theo['conditionalVariance']  =   
					paramInfo$theo['variance']*
					(4*pi*param['shape'] *(a-4)^(param['shape'] ))
		} else {
			paramInfo$theo['conditionalVariance'] =   
					paramInfo$theo['variance']*(4*pi)
		}
		
		
		
		
		
		if(min(c(Nx,Ny)<3*param['rangeInCells'])){
			warning("grid is ", Nx, " by ", Ny,
					"which may be too small for range ",
					paramInfo$theo['rangeInCells'], " cells")
		}
		
 

	##############
	## not range or oneminusar
	#######################
	} else {
		warning("param must have elements named shape  and either oneminusar or range")
		print(param)
		a = NULL
	}
 	
	# build the matrix
	precEntries=getPrec(param['shape'],a)
		
	precEntries =  precEntries /
			as.numeric(paramInfo$theo['conditionalVariance']) 
 

 theNNmat = N
 
 
 theN = theNNmat@x
 theN = precEntries[theN]
 theNNmat@x = theN
 
 
 ######### optimization to see if there are better matern
 # parameters than the theoretical values
 
 varMid = solve(theNNmat,midVec)
 
 theX = distVecFull * paramInfo$theo['cellSize']
 toKeep = which(theX<
				 1.75*paramInfo$theo['range'])
 ev = data.frame(
		 x=theX[toKeep], 
		 y=as.vector(varMid[toKeep]))
 
 ev = tapply(ev$y, ev$x, mean)
 paramInfo$empirical = data.frame(x=
				 as.numeric(names(ev)),
		 y = ev)
 
 paramInfo$empirical = paramInfo$empirical[
		 order(paramInfo$empirical$x),
 ]

 # optimal parameters so product of gmrf precision with matern is identity
 
 ndist = c('1'=0, '2'=1, '3'=sqrt(2), '4'=2, '5'=sqrt(5), '6'=3) 
 numbn = c('1'=1,'2'=4,'3'=4,'4'=4,'5'=8,'6'=4)
 
 nbcells = as.vector(outer((-3):3, 1i*(-3):3,"+"))
 nbcells = nbcells[signif(Mod(nbcells),3)%in%signif(ndist,3)]
 nbn = names(ndist)[match(signif(Mod(nbcells),3), signif(ndist,3))]
 
 nbprec = precEntries[nbn]
 
 othercells = as.vector(outer(0:15, 1i*1:15,"+"))
 othercells = othercells[Re(othercells) <= Im(othercells)]
 otherdist = Mod(outer(othercells, nbcells, '-')) * param['cellSize']
 ndist = ndist * param['cellSize']
 
 # inverse only vary shape
 startparam = paramInfo$theo['shape']
 scaleparam = pmax(startparam,0.0001)
 
 objfun = function(qq) {
	 qq = c(qq,paramInfo$theo[c('variance','range')])
	 
	 othermatern = matern(otherdist, param= qq)
	 
	 nmatern = matern(ndist, param=qq)
	 
	 nbprod = as.vector(othermatern %*% nbprec)
	 thediag = sum(nmatern * precEntries * numbn)
	 sum(nbprod^2) + 2*length(othercells)*(thediag-1)^2
	 
 }
 
 optInv = optim(startparam, objfun,
		 lower=startparam/4,
		 upper=4*startparam,
		 control=list(parscale=scaleparam)
		 ,
		 method='L-BFGS-B')
 
 paramInfo$optimalShape = 
		 c(optInv$par, 
				 paramInfo$theo[c('variance','range',
								 'cellSize')])
 
 
 # now estimate variance separately
 objfun = function(qq) {
	 qq = c(qq,paramInfo$theo[c('variance','range')])
	 othermatern = matern(otherdist, param= qq)
	 nbprod = as.vector(othermatern %*% nbprec)
	 sum(nbprod^2) 
 }
 
 optInv = optim(startparam, objfun,
		 lower=startparam/4,
		 upper=4*startparam,
		 control=list(parscale=scaleparam)
		 ,
		 method='L-BFGS-B')

 nmatern = matern(ndist,
		 param=c(optInv$par, 
				 paramInfo$theo[c('range','variance')])
 )
 thediag = sum(nmatern * precEntries * numbn)
 
 paramInfo$optimal = c(optInv$par,
		 paramInfo$theo[c('cellSize','range')],
		 paramInfo$theo['variance']/as.numeric(thediag))

 
 paramInfo$empirical$optimal =
		 matern(paramInfo$empirical$x,
				 param=paramInfo$optimal)
 
 for(D in c('theo', 
		 grep("^optimal", names(paramInfo),value=TRUE))
  ) {

	paramInfo$empirical[,D] =
		matern(paramInfo$empirical$x,
				param=paramInfo[[D]])

}
 
 
 if(FALSE) { # unused parameter optimizations
 startparam = paramInfo$theo[c('shape',
				 'range',
				 'variance')]
 scaleparam = pmax(startparam,0.0001)
 scaleparam['shape'] = 0.1
 
 newPar = optim(
		 startparam, 
		 function(param){
			 sum((
								 (matern(paramInfo$empirical$x,param=param)) -
								 (paramInfo$empirical$y)
								 )^2)
		 },
		 lower=startparam*0.75,
		 upper=startparam*1.5,
		 method='L-BFGS-B',
		 control=list(parscale=
						 scaleparam)
 )
 
 paramInfo$optimal = 
		 newPar$par
 
 paramInfo$empirical$optimal =
		 matern(paramInfo$empirical$x,
				 param=paramInfo$optimal)
 
 
 # optimal with same shape
startparam = paramInfo$theo[c('range',
				'variance')]
scaleparam = pmax(startparam,0.0001)
newPar = optim(
		startparam, 
		function(param){
			param = c(param, paramInfo$theo['shape'])
			sum((
								(matern(paramInfo$empirical$x,param=param)) -
								(paramInfo$empirical$y)
								)^2)
		},
		lower=startparam*0.75,
		upper=startparam*1.5,
		method='L-BFGS-B',
		control=list(parscale=
						scaleparam)
)

paramInfo$optimalSameShape = 
		c(newPar$par, paramInfo$theo['shape'])


paramInfo$empirical$optimalSameShape =
		matern(paramInfo$empirical$x,
				param=paramInfo$optimalSameShape)


 
 
 thezero = which(paramInfo$empirical$x!=0)
 thisx = paramInfo$empirical$x[thezero]
 thisy = paramInfo$empirical$y[thezero]
 startparam = paramInfo$optimal
 scaleparam = pmax(startparam,0.0001)
 scaleparam['shape'] = 0.1
 
 newPar = optim(
		 startparam, 
		 function(param){
			 sum((
								 matern(thisx,param=param) -
								 thisy)^2)
		 },
		 lower=startparam*0.75,
		 upper=startparam*1.5,
		 method='L-BFGS-B',
		 control=list(parscale=
						 scaleparam)
 )
 
 paramInfo$optimalWithNugget = 
		 newPar$par
 paramInfo$optimalWithNugget['nugget'] =
		 mean(paramInfo$empirical[paramInfo$empirical$x==0,'y']-
						 paramInfo$optimalWithNugget['variance']								
		 )
 
 paramInfo$empirical$nugget =
		 matern(paramInfo$empirical$x,
				 param=paramInfo$optimalWithNugget)
 
 paramInfo$empirical[paramInfo$empirical$x==0,'nugget'] =
		 sum(paramInfo$optimalWithNugget[c('variance','nugget')])

 
 
 # optimal parameters so product of gmrf precision with matern is identity

ndist = c('1'=0, '2'=1, '3'=sqrt(2), '4'=2, '5'=sqrt(5), '6'=3) 
numbn = c('1'=1,'2'=4,'3'=4,'4'=4,'5'=8,'6'=4)

nbcells = as.vector(outer((-3):3, 1i*(-3):3,"+"))
nbcells = nbcells[signif(Mod(nbcells),3)%in%signif(ndist,3)]
nbn = names(ndist)[match(signif(Mod(nbcells),3), signif(ndist,3))]

nbprec = precEntries[nbn]

othercells = as.vector(outer(0:12, 1i*1:12,"+"))
othercells = othercells[Re(othercells) <= Im(othercells)]
otherdist = Mod(outer(othercells, nbcells, '-')) * param['cellSize']
ndist = ndist * param['cellSize']


objfun = function(qq) {

	othermatern = matern(otherdist, param= qq)
	
	nmatern = matern(ndist, param=qq)
 
	nbprod = as.vector(othermatern %*% nbprec)
	thediag = sum(nmatern * precEntries * numbn)
 	sum(nbprod^2) + 2*length(othercells)*(thediag-1)^2
	
}

startparam = paramInfo$optimal[c('variance','range','shape')]
scaleparam = pmax(startparam,0.0001)
scaleparam['shape'] = 0.1

optInv = optim(startparam, objfun,
		lower=0.75*startparam,
		upper=1.5*startparam,
		control=list(parscale=scaleparam)
						,
				method='L-BFGS-B')
		
paramInfo$optimalInverse = optInv$par

paramInfo$empirical$optimalInverse =
		matern(paramInfo$empirical$x,
				param=paramInfo$optimalInverse)

# inverse, keep shape

startparam = paramInfo$optimalSameShape[c('variance','range')]
scaleparam = pmax(startparam,0.0001)
 

objfun = function(qq) {
	qq['shape'] = paramInfo$theo['shape']
	
	othermatern = matern(otherdist, param= qq)
	
	nmatern = matern(ndist, param=qq)
	
	nbprod = as.vector(othermatern %*% nbprec)
	thediag = sum(nmatern * precEntries * numbn)
	sum(nbprod^2) + 2*length(othercells)*(thediag-1)^2
	
}

optInv = optim(startparam, objfun,
		lower=0.75*startparam,
		upper=1.5*startparam,
		control=list(parscale=scaleparam)
		,
		method='L-BFGS-B')

paramInfo$optimalInverseSameShape = 
		c(optInv$par, paramInfo$theo['shape'])

paramInfo$empirical$optimalInverseSameShape =
		matern(paramInfo$empirical$x,
				param=paramInfo$optimalInverseSameShape)

# inverse keep range
# inverse, keep shape

startparam = paramInfo$theo[c('variance','shape')]
scaleparam = pmax(startparam,0.0001)


objfun = function(qq) {
	qq['range'] = paramInfo$theo['range']
	
	othermatern = matern(otherdist, param= qq)
	
	nmatern = matern(ndist, param=qq)
	
	nbprod = as.vector(othermatern %*% nbprec)
	thediag = sum(nmatern * precEntries * numbn)
	sum(nbprod^2) + 2*length(othercells)*(thediag-1)^2
	
}

optInv = optim(startparam, objfun,
		lower=0.5*startparam,
		upper=2*startparam,
		control=list(parscale=scaleparam)
		,
		method='L-BFGS-B')

paramInfo$optimalInverseSameRange = 
		c(optInv$par, paramInfo$theo['range'])

paramInfo$empirical$optimalInverseSameRange =
		matern(paramInfo$empirical$x,
				param=paramInfo$optimalInverseSameRange)

# inverse only vary shape
startparam = paramInfo$theo['shape']
scaleparam = pmax(startparam,0.0001)

objfun = function(qq) {
	qq = c(qq,paramInfo$theo[c('variance','range')])
	
	othermatern = matern(otherdist, param= qq)
	
	nmatern = matern(ndist, param=qq)
	
	nbprod = as.vector(othermatern %*% nbprec)
	thediag = sum(nmatern * precEntries * numbn)
	sum(nbprod^2) + 2*length(othercells)*(thediag-1)^2
	
}
		
optInv = optim(startparam, objfun,
		lower=0.5*startparam,
		upper=2*startparam,
		control=list(parscale=scaleparam)
		,
		method='L-BFGS-B')

paramInfo$optimalInverseShape = 
		c(optInv$par, paramInfo$theo[c('variance','range')])

paramInfo$empirical$optimalInverseShape =
		matern(paramInfo$empirical$x,
				param=paramInfo$optimalInverseShape)


# only shape
startparam = paramInfo$theo['shape']
scaleparam = pmax(startparam,0.0001)
newPar = optim(
		startparam, 
		function(param){
			param = c(param, 
					paramInfo$theo[c('variance','range')])
			sum((
								(matern(paramInfo$empirical$x,param=param)) -
								(paramInfo$empirical$y)
								)^2)
		},
		lower=startparam*0.5,
		upper=startparam*2,
		method='L-BFGS-B',
		control=list(parscale=
						scaleparam)
)

paramInfo$optimalShape = 
		c(newPar$par, 
				paramInfo$theo[c('variance','range')])


paramInfo$empirical$optimalShape =
		matern(paramInfo$empirical$x,
				param=paramInfo$optimalShape)

}

 # edge correction
		
		
	if(adjustEdges!=FALSE){
		if(is.logical(adjustEdges)) {
			adjustEdges='optimalShape'
		}
		
		distVecFull =  expand.grid(x=seq(1, Nx),
				y=seq(Ny,1))
		buffer = param['shape']+1
		
		cellVec = seq(1, nrow(distVecFull)) 
		
		whichInner = 
				distVecFull[,'x']> buffer & 
						distVecFull[,'y']>buffer &
						distVecFull[,'x']< (Nx-buffer) &
						distVecFull[,'y'] < (Ny-buffer) 
			
		innerCells = cellVec[whichInner]
		outerCells = cellVec[!whichInner]

		outerCoordsCartesian = SpatialPoints(
			distVecFull[outerCells,]*param['cellSize']
		)				
 		
		paramForM = paramInfo[[adjustEdges]]

    covMatInv = matern(
        outerCoordsCartesian,
				param= paramForM,
        type='precision')

		InnerPrecision = theNNmat[innerCells, innerCells]
			
			#A = x[allCells,-allCells]
			#InnerPrecInvChol = Matrix::solve(Matrix::chol(InnerPrecision))
			#Aic = A %*% InnerPrecInvChol
			# AQinvA = Aic %*% t(Aic)
		A = theNNmat[innerCells,outerCells]
		cholInnerPrec =Cholesky(InnerPrecision,LDL=FALSE,perm=TRUE)

		Aic = solve(cholInnerPrec, 
					solve(cholInnerPrec,A,system='P'),
					system='L')

    AQinvA = crossprod(Aic)
    AQinvA = as(AQinvA, class(covMatInv))
    
    precOuter = covMatInv + AQinvA
    precOuter = forceSymmetric(precOuter)
    precOuter = as.matrix(precOuter)
    
		theNNmat[outerCells,outerCells] = precOuter
		theNNmat = forceSymmetric(theNNmat)
		} # end edge adjustment 
		paramInfo$adjustEdges=adjustEdges
		
		paramInfo$precisionEntries = precEntries

		theNNmat = drop0(theNNmat)

		
		attributes(theNNmat)$param=
				paramInfo

		attributes(theNNmat)$raster= theraster
		
	return(theNNmat)
}
	

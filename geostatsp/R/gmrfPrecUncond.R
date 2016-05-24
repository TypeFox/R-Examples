gmrfPrecUncond = function(x, 
			N = attributes(x)$Nx, Ny=attributes(x)$Ny,
			param = attributes(x)$model,
			border=param["shape"]+1){
		if(!length(N) | !length(Ny) | !length(border) | !length(param))	
			warning("border, N, Ny, or param were not supplied")
		if(border < 1) warning("border should be >1")
		
		if(!all(c("scale","prec","shape","cellSize") %in% names(param))) {
			warning("param needs scale, prec, shape, cellSize")
			print(param)
		}
		
	Nx=N
		
		
	Ncells = Nx*Ny
	
	topCoords = rep(1:Nx, border) + 1i * rep(1:border, rep(Nx, border))
	topCells = Re(topCoords) + (Im(topCoords)-1)*Nx
	
	bottomCoords = rep(1:Nx, border) + 1i * rep(seq(Ny-border+1, Ny), rep(Nx, border))
	bottomCells = Re(bottomCoords) + (Im(bottomCoords)-1)*Nx
	
	Nonside = Ny-2*border
	leftCoords = rep(1:border, Nonside) + 
			1i * rep(seq(border+1, len=Nonside), rep(border,Nonside))
	leftCells = Re(leftCoords) + (Im(leftCoords)-1)*Nx
	rightCoords = rep(seq(Nx-border+1, Nx), Nonside) + 
			1i * rep(seq(border+1, len=Nonside), rep(border,Nonside))
	rightCells = Re(rightCoords) + (Im(rightCoords)-1)*Nx
	
	allCoords = c(topCoords, leftCoords, rightCoords, bottomCoords)
	covMat = matern(
			SpatialPoints(cbind(Re(allCoords),Im(allCoords))),
			param=c(param[c('shape','variance')],
					param['range']/param['cellSize'])
			)
	
	

	#	covChol = chol(covMat)
	#	covInvChol = solve(covChol)
	#	precOuter = solve(covMat)
	
	allCells = c(topCells, leftCells, rightCells, bottomCells)
	InnerPrecision = x[-allCells, -allCells]
	
	#A = x[allCells,-allCells]
	#InnerPrecInvChol = Matrix::solve(Matrix::chol(InnerPrecision))
	#Aic = A %*% InnerPrecInvChol
	# AQinvA = Aic %*% t(Aic)
	
	cholInnerPrec = chol(InnerPrecision)
	A = x[-allCells,allCells]
	Aic = solve(cholInnerPrec, A)
	AQinvA = crossprod(Aic,Aic)

	precOuter = Matrix::solve(covMat) + AQinvA

	x[allCells,allCells] = precOuter
	x
}

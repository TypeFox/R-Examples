


########################################################################################################################
########################################  Class definition #############################################################
########################################################################################################################

################ Auxiliar


setClass(
	Class= "SpatC",
	representation=representation(
		Matalpha="list",
		Matbeta	="list",
		Matphi	="list",
		Mattheta="list"
	)
)

setClass(
	Class= "SpatH",
	representation=representation(
		west		="list",
		east		="list",	
		north		="list",
		south		="list",		
		southwest	="list",
		northeast	="list",
		southeast	="list",
		northwest	="list"
	)
)


setClass(
	Class= "CosSinMatrix",
	representation=representation(
		cosMat		="matrix",
		sinMat		="matrix",
		littlecosMat="matrix",
		littlesinMat="matrix"
	)
)

setClass(
	Class= "MtAux",
	representation=representation(
		seas="list"
	)
)





################ Parameters

setClass(
	Class= "SpatParam",
	representation=representation(
			alpha	="vector",
			beta	="vector",
			phi		="vector",
			theta	="vector",
			Cmat	="matrix",
			lags	="vector",
			dirs	="vector"
	)
)


setClass(
	Class= "VectSubdiag",
	representation=representation(
			east		="vector",
			west		="vector",
			north		="vector",
			south		="vector",
			southeast	="vector",
			northwest	="vector",
			southwest	="vector",
			northeast	="vector",
			lags		="vector",
			dirs		="vector"
	)
)

setClass(
	Class= "Seas",
	representation=representation(
			w			="numeric",
			fvect		="vector",
			f0L			="matrix",
			gvect		="vector",
			g0L			="matrix"
	)
)

setClass(
	Class= "Autoregressive",
	representation=representation(
			avect		="matrix",
			a0vect		="matrix",
			a0L			="matrix",
			spatialA	="SpatParam",
			sigma2A		="numeric",
			H			="matrix",
			subdiag		="VectSubdiag"
	)
)

setClass(
	Class= "Mu",
	representation=representation(
			muvect		="matrix",
			mu0vect		="matrix",
			mu0L		="matrix",
			sigma2Mu		="numeric",
			spatialMu	="SpatParam"
	)
)

setClass(
	Class= "Mt",
	representation=representation(
			Mt	="matrix",
			seas="list"
	)
)

setClass(
	Class= "Xt",
	representation=representation(
			Xt		="matrix",
			X0		="matrix",
			sigma2N	="numeric",
			AR		="list",
			templags="vector"
	)
)

setClass(
	Class= "Parameters",
	representation=representation(
			Mu		="Mu",
			Mt		="Mt",
			Xt		="Xt",
			Yt		="matrix",
			sigma2Y	="numeric",
			sigma2E	="numeric"
	)
)


################ Hyperpriors

setClass(
	Class= "SpatParam0",
	representation=representation(
			alpha0	="matrix",
			beta0	="matrix",
			phi0	="matrix",
			theta0	="matrix"
	)
)

setClass(
	Class= "VectSubdiag0",
	representation=representation(
			east0		="matrix",
			west0		="matrix",
			north0		="matrix",
			south0		="matrix",
			southeast0	="matrix",
			northwest0	="matrix",
			southwest0	="matrix",
			northeast0	="matrix"
	)
)

setClass(
	Class= "Seas0",
	representation=representation(
			f0L0 	="vector",
			sigf0L0	="matrix",
			g0L0	="vector",
			sigg0L0	="matrix"
	)
)

setClass(
	Class= "Autoregressive0",
	representation=representation(
			a0L0 		="vector",
			siga0L0		="matrix",
			sigma2A0		="vector",
			spatialA0	="SpatParam0",
			subdiag0	="VectSubdiag0"
	)
)

setClass(
	Class= "Mu0",
	representation=representation(
			mu0L0 		="vector",
			sigmu0L0	="matrix",
			sigma2Mu0	="vector",
			spatialMu0	="SpatParam0"
	)
)

setClass(
	Class= "Mt0",
	representation=representation(
			seas0	="list"
	)
)

setClass(
	Class= "Xt0",
	representation=representation(
			X00		="matrix",
			sigma2X00="matrix",
			AR0		="list",
			sigma2N0	="vector"
	)
)

setClass(
	Class= "Hyperpriors",
	representation=representation(
			Mu0		="Mu0",
			Mt0		="Mt0",
			Xt0		="Xt0",
			sigma2Y0	="vector",
			sigma2E0	="vector"
	)
)


######################

setClass(
	Class= "HBSTM",
	representation=representation(
			Parameters	="Parameters",
			Hyperpriors	="Hyperpriors",
			mse			="vector",
			seed		="numeric",
			iterations 	="numeric",
			newGrid		="matrix",
			K			="matrix",
			Zt			="matrix",
			fitted		="array",
			residuals	="matrix",
			MCMCsamp	="list",
			MCMCclass 	="character"
	)
)

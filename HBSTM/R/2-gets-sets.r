
############################################################################################################
################################################ Gets and sets #############################################
############################################################################################################

################ Auxiliar

###
setMethod(
	f="[",
	signature=c("SpatC","character","missing","missing"),
	def = function(x,i,j,drop){
		switch(EXP=i,
			Matalpha	= return(x@Matalpha),
			Matbeta		= return(x@Matbeta),
			Matphi 		= return(x@Matphi),
			Mattheta 	= return(x@Mattheta),
			stop("Error: ",i," is not a SpatC slot")
		)
	}
)

setMethod(
	f="[<-",
	signature=c("SpatC","character","missing","ANY"),
	def = function(x,i,j,value){
		switch(EXP=i,
			Matalpha	= x@Matalpha <- value,
			Matbeta		= x@Matbeta <- value,
			Matphi		= x@Matphi <- value,
			Mattheta	= x@Mattheta <- value
		)
		return(x)
	}
)

###
setMethod(
	f="[",
	signature=c("SpatH","character","missing","missing"),
	def = function(x,i,j,drop){
		switch(EXP=i,
			west 		= return(x@west),
			east 		= return(x@east),
			north		= return(x@north),
			south		= return(x@south),
			southwest 	= return(x@southwest),	
			northeast 	= return(x@northeast),
			southeast 	= return(x@southeast),
			northwest 	= return(x@northwest),
			stop("Error: ",i," is not a SpatH slot")
		)
	}
)

setMethod(
	f="[<-",
	signature=c("SpatH","character","missing","ANY"),
	def = function(x,i,j,value){
		switch(EXP=i,
			west		= x@west <- value,
			east		= x@east <- value,
			north		= x@north <- value,
			south		= x@south <- value,
			southwest	= x@southwest <- value,
			northeast	= x@northeast <- value,
			southeast	= x@southeast <- value,
			northwest	= x@northwest <- value
		)
		return(x)
	}
)

###
setMethod(
	f="[",
	signature=c("CosSinMatrix","character","missing","missing"),
	def = function(x,i,j,drop){
		switch(EXP=i,
			cosMat			= return(x@cosMat),
			sinMat			= return(x@sinMat),
			littlecosMat 	= return(x@littlecosMat),
			littlesinMat 	= return(x@littlesinMat),
			stop("Error: ",i," is not a CosSinMatrix slot")
		)
	}
)

setMethod(
	f="[<-",
	signature=c("CosSinMatrix","character","missing","ANY"),
	def = function(x,i,j,value){
		switch(EXP=i,
			cosMat			= x@cosMat <- value,
			sinMat			= x@sinMat <- value,
			littlecosMat	= x@littlecosMat <- value,
			littlesinMat	= x@littlesinMat <- value
		)
		return(x)
	}
)

###
setMethod(
	f="[",
	signature=c("MtAux","character","missing","missing"),
	def = function(x,i,j,drop){
		switch(EXP=i,
			seas	= return(x@seas),
			stop("Error: ",i," is not a MtAux slot")
		)
	}
)

setMethod(
	f="[<-",
	signature=c("MtAux","character","missing","ANY"),
	def = function(x,i,j,value){
		switch(EXP=i,
			seas	= x@seas <- value
		)
		return(x)
	}
)

###
setMethod(
	f="[",
	signature=c("HBSTM","character","missing","missing"),
	def = function(x,i,j,drop){
		switch(EXP=i,
			Parameters	= return(x@Parameters),
			Hyperpriors	= return(x@Hyperpriors),
			seed		= return(x@seed),
			mse			= return(x@mse),
			iterations	= return(x@iterations),
			newGrid		= return(x@newGrid),
			K 			= return(x@K),
			Zt 			= return(x@Zt),
			fitted 		= return(x@fitted),
			residuals	= return(x@residuals),
			MCMCsamp 	= return(x@MCMCsamp),
			MCMCclass 	= return(x@MCMCclass),
			stop("Error: ",i," is not a STmodel slot")
		)
	}
)

setMethod(
	f="[<-",
	signature=c("HBSTM","character","missing","ANY"),
	def = function(x,i,j,value){
		switch(EXP=i,
			Parameters	= x@Parameters <- value,
			Hyperpriors	= x@Hyperpriors <- value,
			seed		= x@seed <- value,
			mse			= x@mse <- value,
			iterations	= x@iterations <- value,
			newGrid		= x@newGrid <- value,
			K			= x@K <- value,
			Zt			= x@Zt <- value,
			fitted		= x@fitted <- value,
			residuals	= x@residuals <- value,
			MCMCsamp	= x@MCMCsamp <- value,
			MCMCclass	= x@MCMCclass <- value
		)
		return(x)
	}
)

################ Parameters

###
setMethod(
	f="[",
	signature=c("SpatParam","character","missing","missing"),
	def = function(x,i,j,drop){
		switch(EXP=i,
			alpha	= return(x@alpha),
			beta	= return(x@beta),
			phi 	= return(x@phi),
			theta 	= return(x@theta),
			Cmat 	= return(x@Cmat),
			lags 	= return(x@lags),
			dirs 	= return(x@dirs),
			stop("Error: ",i," is not a SpatParam slot")
		)
	}
)

setMethod(
	f="[<-",
	signature=c("SpatParam","character","missing","ANY"),
	def = function(x,i,j,value){
		switch(EXP=i,
			alpha	= x@alpha <- value,
			beta	= x@beta <- value,
			phi		= x@phi <- value,
			theta	= x@theta <- value,
			Cmat	= x@Cmat <- value,
			lags	= x@lags <- value,
			dirs	= x@dirs <- value
		)
		return(x)
	}
)

###
setMethod(
	f="[",
	signature=c("VectSubdiag","character","missing","missing"),
	def = function(x,i,j,drop){
		switch(EXP=i,
			east		= return(x@east),
			west		= return(x@west),
			north 		= return(x@north),
			south 		= return(x@south),
			southeast	= return(x@southeast),
			northwest	= return(x@northwest),
			southwest 	= return(x@southwest),
			northeast 	= return(x@northeast),
			lags 		= return(x@lags),
			dirs 		= return(x@dirs),
			stop("Error: ",i," is not a VectSubdiag slot")
		)
	}
)

setMethod(
	f="[<-",
	signature=c("VectSubdiag","character","missing","ANY"),
	def = function(x,i,j,value){
		switch(EXP=i,
			east		= x@east <- value,
			west		= x@west <- value,
			north		= x@north <- value,
			south		= x@south <- value,
			southeast	= x@southeast <- value,
			northwest	= x@northwest <- value,
			southwest	= x@southwest <- value,
			northeast	= x@northeast <- value,
			lags		= x@lags <- value,
			dirs		= x@dirs <- value
		)
		return(x)
	}
)

###
setMethod(
	f="[",
	signature=c("Seas","character","missing","missing"),
	def = function(x,i,j,drop){
		switch(EXP=i,
			w		= return(x@w),
			fvect	= return(x@fvect),
			f0L 	= return(x@f0L),
			gvect	= return(x@gvect),
			g0L		= return(x@g0L),
			stop("Error: ",i," is not a Seas slot")
		)
	}
)

setMethod(
	f="[<-",
	signature=c("Seas","character","missing","ANY"),
	def = function(x,i,j,value){
		switch(EXP=i,
			w		= x@w <- value,
			fvect	= x@fvect <- value,
			f0L		= x@f0L <- value,
			gvect	= x@gvect <- value,
			g0L		= x@g0L <- value
		)
		return(x)
	}
)

###
setMethod(
	f="[",
	signature=c("Autoregressive","character","missing","missing"),
	def = function(x,i,j,drop){
		switch(EXP=i,
			avect	= return(x@avect),
			a0vect	= return(x@a0vect),
			a0L 	= return(x@a0L),
			spatialA= return(x@spatialA),
			sigma2A	= return(x@sigma2A),
			H		= return(x@H),
			subdiag = return(x@subdiag),
			stop("Error: ",i," is not a Autoregressive slot")
		)
	}
)

setMethod(
	f="[<-",
	signature=c("Autoregressive","character","missing","ANY"),
	def = function(x,i,j,value){
		switch(EXP=i,
			avect	= x@avect <- value,
			a0vect	= x@a0vect <- value,
			a0L		= x@a0L <- value,
			spatialA= x@spatialA <- value,
			sigma2A	= x@sigma2A <- value,
			H		= x@H <- value,
			subdiag	= x@subdiag <- value
		)
		return(x)
	}
)

###
setMethod(
	f="[",
	signature=c("Mu","character","missing","missing"),
	def = function(x,i,j,drop){
		switch(EXP=i,
			muvect		= return(x@muvect),
			mu0vect		= return(x@mu0vect),
			mu0L 		= return(x@mu0L),
			sigma2Mu 	= return(x@sigma2Mu),
			spatialMu	= return(x@spatialMu),
			stop("Error: ",i," is not a Mu slot")
		)
	}
)

setMethod(
	f="[<-",
	signature=c("Mu","character","missing","ANY"),
	def = function(x,i,j,value){
		switch(EXP=i,
			muvect		= x@muvect <- value,
			mu0vect		= x@mu0vect <- value,
			mu0L		= x@mu0L <- value,
			sigma2Mu		= x@sigma2Mu <- value,
			spatialMu	= x@spatialMu <- value
		)
		return(x)
	}
)

###
setMethod(
	f="[",
	signature=c("Mt","character","missing","missing"),
	def = function(x,i,j,drop){
		switch(EXP=i,
			Mt		= return(x@Mt),
			seas	= return(x@seas),
			stop("Error: ",i," is not a Mt slot")
		)
	}
)

setMethod(
	f="[<-",
	signature=c("Mt","character","missing","ANY"),
	def = function(x,i,j,value){
		switch(EXP=i,
			Mt		= x@Mt <- value,
			seas	= x@seas <- value
		)
		return(x)
	}
)

###
setMethod(
	f="[",
	signature=c("Xt","character","missing","missing"),
	def = function(x,i,j,drop){
		switch(EXP=i,
			Xt		= return(x@Xt),
			X0		= return(x@X0),
			sigma2N 	= return(x@sigma2N),
			AR 		= return(x@AR),
			templags= return(x@templags),
			stop("Error: ",i," is not a Xt slot")
		)
	}
)

setMethod(
	f="[<-",
	signature=c("Xt","character","missing","ANY"),
	def = function(x,i,j,value){
		switch(EXP=i,
			Xt		= x@Xt <- value,
			X0		= x@X0 <- value,
			sigma2N	= x@sigma2N <- value,
			AR		= x@AR <- value,
			templags= x@templags <- value
		)
		return(x)
	}
)

###
setMethod(
	f="[",
	signature=c("Parameters","character","missing","missing"),
	def = function(x,i,j,drop){
		switch(EXP=i,
			Mu		= return(x@Mu),
			Mt		= return(x@Mt),
			Xt		= return(x@Xt),
			Yt 		= return(x@Yt),
			sigma2Y 	= return(x@sigma2Y),
			sigma2E	= return(x@sigma2E),
			stop("Error: ",i," is not a Parameters slot")
		)
	}
)

setMethod(
	f="[<-",
	signature=c("Parameters","character","missing","ANY"),
	def = function(x,i,j,value){
		switch(EXP=i,
			Mu		= x@Mu <- value,
			Mt		= x@Mt <- value,
			Xt		= x@Xt <- value,
			Yt		= x@Yt <- value,
			sigma2Y	= x@sigma2Y <- value,
			sigma2E	= x@sigma2E <- value
		)
		return(x)
	}
)


################ Hyperpriors

###
setMethod(
	f="[",
	signature=c("SpatParam0","character","missing","missing"),
	def = function(x,i,j,drop){
		switch(EXP=i,
			alpha0	= return(x@alpha0),
			beta0	= return(x@beta0),
			phi0	= return(x@phi0),
			theta0 	= return(x@theta0),
			stop("Error: ",i," is not a SpatParam0 slot")
		)
	}
)

setMethod(
	f="[<-",
	signature=c("SpatParam0","character","missing","ANY"),
	def = function(x,i,j,value){
		switch(EXP=i,
			alpha0	= x@alpha0 <- value,
			beta0	= x@beta0 <- value,
			phi0	= x@phi0 <- value,
			theta0	= x@MatTheta <- value
		)
		return(x)
	}
)

###
setMethod(
	f="[",
	signature=c("VectSubdiag0","character","missing","missing"),
	def = function(x,i,j,drop){
		switch(EXP=i,
			east0		= return(x@east0),
			west0		= return(x@west0),
			north0 		= return(x@north0),
			south0		= return(x@south0),
			southeast0	= return(x@southeast0),
			northwest0	= return(x@northwest0),
			southwest0 	= return(x@southwest0),
			northeast0 	= return(x@northeast0),
			stop("Error: ",i," is not a VectSubdiag0 slot")
		)
	}
)

setMethod(
	f="[<-",
	signature=c("VectSubdiag0","character","missing","ANY"),
	def = function(x,i,j,value){
		switch(EXP=i,
			east0		= x@east0 <- value,
			west0		= x@west0 <- value,
			north0		= x@north0 <- value,
			south0		= x@south0 <- value,
			southeast0	= x@southeast0 <- value,
			northwest0	= x@northwest0 <- value,
			southwest0	= x@southwest0 <- value,
			northeast0	= x@northeast0 <- value
		)
		return(x)
	}
)

###
setMethod(
	f="[",
	signature=c("Seas0","character","missing","missing"),
	def = function(x,i,j,drop){
		switch(EXP=i,
			sigf0L0	= return(x@sigf0L0),
			f0L0	= return(x@f0L0),
			sigg0L0 = return(x@sigg0L0),
			g0L0	= return(x@g0L0),
			stop("Error: ",i," is not a Seas0 slot")
		)
	}
)

setMethod(
	f="[<-",
	signature=c("Seas0","character","missing","ANY"),
	def = function(x,i,j,value){
		switch(EXP=i,
			sigf0L0	= x@sigf0L0 <- value,
			f0L0	= x@f0L0 <- value,
			sigg0L0	= x@sigg0L0 <- value,
			g0L0	= x@g0L0 <- value
		)
		return(x)
	}
)

###
setMethod(
	f="[",
	signature=c("Autoregressive0","character","missing","missing"),
	def = function(x,i,j,drop){
		switch(EXP=i,
			a0L0 		= return(x@a0L0),
			siga0L0 	= return(x@siga0L0),
			sigma2A0		= return(x@sigma2A0),
			spatialA0 	= return(x@spatialA0),
			subdiag0	= return(x@subdiag0),
			stop("Error: ",i," is not a Autoregressive0 slot")
		)
	}
)

setMethod(
	f="[<-",
	signature=c("Autoregressive0","character","missing","ANY"),
	def = function(x,i,j,value){
		switch(EXP=i,
			a0L0		= x@a0L0 <- value,
			siga0L0		= x@siga0L0 <- value,
			sigma2A0		= x@sigma2A0 <- value,
			spatialA0	= x@spatialA0 <- value,
			subdiag0	= x@subdiag0 <- value
		)
		return(x)
	}
)

###
setMethod(
	f="[",
	signature=c("Mu0","character","missing","missing"),
	def = function(x,i,j,drop){
		switch(EXP=i,
			mu0L0 		= return(x@mu0L0),
			sigmu0L0 	= return(x@sigmu0L0),
			sigma2Mu0	= return(x@sigma2Mu0),
			spatialMu0 	= return(x@spatialMu0),
			stop("Error: ",i," is not a Mu0 slot")
		)
	}
)

setMethod(
	f="[<-",
	signature=c("Mu0","character","missing","ANY"),
	def = function(x,i,j,value){
		switch(EXP=i,
			mu0L0		= x@mu0L0 <- value,
			sigmu0L0	= x@sigmu0L0 <- value,
			sigma2Mu0	= x@sigma2Mu0 <- value,
			spatialMu0	= x@spatialMu0 <- value
		)
		return(x)
	}
)

###
setMethod(
	f="[",
	signature=c("Mt0","character","missing","missing"),
	def = function(x,i,j,drop){
		switch(EXP=i,
			seas0 	= return(x@seas0),
			stop("Error: ",i," is not a Mt0 slot")
		)
	}
)

setMethod(
	f="[<-",
	signature=c("Mt0","character","missing","ANY"),
	def = function(x,i,j,value){
		switch(EXP=i,
			seas0	= x@seas0 <- value
		)
		return(x)
	}
)

###
setMethod(
	f="[",
	signature=c("Xt0","character","missing","missing"),
	def = function(x,i,j,drop){
		switch(EXP=i,
			X00 	= return(x@X00),
			sigma2X00= return(x@sigma2X00),
			AR0 	= return(x@AR0),
			sigma2N0 = return(x@sigma2N0),
			stop("Error: ",i," is not a Xt0 slot")
		)
	}
)

setMethod(
	f="[<-",
	signature=c("Xt0","character","missing","ANY"),
	def = function(x,i,j,value){
		switch(EXP=i,
			X00		= x@X00 <- value,
			sigma2X00= x@sigma2X00 <- value,
			AR0		= x@AR0 <- value,
			sigma2N0	= x@sigma2N0 <- value
		)
		return(x)
	}
)

###
setMethod(
	f="[",
	signature=c("Hyperpriors","character","missing","missing"),
	def = function(x,i,j,drop){
		switch(EXP=i,
			Mu0 	= return(x@Mu0),
			Mt0 	= return(x@Mt0),
			Xt0		= return(x@Xt0),
			sigma2Y0 = return(x@sigma2Y0),
			sigma2E0	= return(x@sigma2E0),
			stop("Error: ",i," is not a Hyperpriors slot")
		)
	}
)

setMethod(
	f="[<-",
	signature=c("Hyperpriors","character","missing","ANY"),
	def = function(x,i,j,value){
		switch(EXP=i,
			Mu0		= x@Mu0 <- value,
			Mt0		= x@Mt0 <- value,
			Xt0		= x@Xt0 <- value,
			sigma2Y0	= x@sigma2Y0 <- value,
			sigma2E0	= x@sigma2E0 <- value
		)
		return(x)
	}
)





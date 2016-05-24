################################## Eup slope inference ####################
## Input:
#     dat.matrix	= the data in matrix form where the first colomn
#				  contain NT vector of Y second one the NT vector 
#				  of the first regressor X
# 		dat.dim	= the dimension of the data N and T
# 		used.d	= used dimension d in the interative procedure
#	 	beta.Eup	= the estimated slope estimator for given d
# 		factors	= common factors after scaling according the used 
#				  restriction
# 		loadings	= individual loading parameters after scaling
#				 (according restriction)
#		residuals	= the residual terms
## Output:
#		Eup slope Estimate
#		std
#		Pr(>|z|)
##########################################################################

Eup.inference <- function(Eup.Obj, error.type, kernel.weights){

## collect informations from the Eup.Obj

	y 	<- Eup.Obj$dat.matrix[, 1, drop = FALSE]
	x 	<- Eup.Obj$dat.matrix[,-1, drop = FALSE]
	nr  <- Eup.Obj$dat.dim[1]
	nc	<- Eup.Obj$dat.dim[2]
	P	 <- Eup.Obj$dat.dim[3]
  	d   <- Eup.Obj$used.dim
	slope <- Eup.Obj$slope.para
  	Intercept <- matrix(Eup.Obj$Intercept, 1, 1)
  	is.intercept <- Eup.Obj$is.intercept
  	OvMeans <- Eup.Obj$OvMeans
	F	 <- Eup.Obj$unob.factors
	A	 <- Eup.Obj$ind.loadings
	sig2.hat <- Eup.Obj$sig2.hat
  	resid  <- Eup.Obj$residuals

  ## weights 
  if(error.type[1] == 5| error.type[1] == 8){
	if(is.null(kernel.weights)) {
		Ltrank = round(0.7*sqrt(nr))
		Weights <- 1 - seq(1:nr)/(Ltrank + 1)
		Weights[Weights < 0] = 0
	}
	else{
	if(!is.numeric(kernel.weights)| length(kernel.weights) >= nr)
	stop("The arugment 'kernel.weights' should be numeric with a lenght no larger than 'T-1'.")
	if(length(kernel.weights) < (nr-1)) {
		Weights <- rep.int(0, nr)
		Weights[1:length(kernel.weights)] <- kernel.weights
		}
	else Weights = kernel.weights
	Ltrank <- length(Weights[Weights > 0])
	}
   }


## if intercept
	if(is.intercept) {
		x <- cbind(rep.int(1, (nr*nc)), x)
		P = P +1 
		rownames(Intercept) <- "(Intercept)"
		colnames(Intercept) <- " "
		slope <- rbind(Intercept, slope)
		} 

## Projection matrix of the factors F
  
  I.TxT <- diag(1, ncol= nr, nrow= nr)
  if(d==0) M.F   <- I.TxT
  else M.F   <- I.TxT - tcrossprod(F)/nr

## Projection matrix of the loadings A

  I.NxN <- diag(1, ncol= nc, nrow= nc)
  if(d==0) M.A   <- I.NxN
  else{
    crossdiagA <- diag(crossprod(A))
    S  <- diag(crossdiagA^{-0.5}, d)
    AS <- tcrossprod(A,S)
    P.A   <- tcrossprod(AS)
    M.A   <- I.NxN - P.A
  }
    

## write the x matrices in a list: each regressor is written in a 
## list component
	X.mat.list <- NULL
	for(p in 1:P) X.mat.list[[p]] <- matrix(x[,p], nr, nc)

  # Z_i = M.F * X_i - sum{M.F * X_k*a_ik}/n
	if(d==0) Z.list <- X.mat.list
	else{
		Z.list <- sapply (X.mat.list, function(X) M.F %*% X %*% M.A 
			    , simplify = FALSE)
	}

  # construct the matrix D= sum Z_i'Z_i/NT
	Z	 <- sapply (Z.list, function(Z) c(Z), simplify = TRUE)
	D0     <- crossprod(Z)/(nr*nc)
	inv.D0 <- solve(D0)


 # construct the asym variance

	Error.typ = error.type[1]
#1 # "iid"
	if(Error.typ == 1){
		asy.var <- inv.D0*sig2.hat/(nr*nc)
		}
#2 # " c-heteros.T/N-->0" 
	if(Error.typ == 2){
		ind.sig2.hat <- diag(sqrt(diag(crossprod(resid))/nr)) # (nc x nc)
		Z.ind.sig2 <- sapply(Z.list, function(X) c(X%*%ind.sig2.hat), simplify = TRUE)
		D1 <- crossprod(Z.ind.sig2)/(nr*nc)
		asy.var <- inv.D0%*%D1%*%inv.D0/(nr*nc)
		}
#3 # "c-corr.T/N-->0"
	if(Error.typ == 3){
		Z.icor.sig2 <- sapply(Z.list, function(X) rowSums(resid[, 1:round(log(nc))]*X[, 1:round(log(nc))]), simplify = TRUE)
		D11 <- crossprod(Z.icor.sig2)/(nr*log(nc))
		asy.var <- inv.D0%*%D11%*%inv.D0/(nr*nc)
		}
#4 # "s-heteros.N/T-->0"
	if(Error.typ == 4){
		tim.sig2.hat <- diag(sqrt(diag(tcrossprod(resid))/nc)) # (nr x nr)
		Z.tim.sig2 <- sapply(Z.list, function(X) c(tim.sig2.hat%*%X), simplify = TRUE)
		D2 <- crossprod(Z.tim.sig2)/(nr*nc)
		asy.var <- inv.D0%*%D2%*%inv.D0/(nr*nc)
		}
#5 # "s-corr.N/T-->0"
	if(Error.typ == 5){
		wtw <- matrix(0, P, P)
		for(s in 2:Ltrank){#s=2
			w1 <- sapply(Z.list, function(Z) c((Z*resid)[1:(nr - s +1) , ]), simplify = TRUE)
			w2 <- sapply(Z.list, function(Z) c((Z*resid)[s:nr , ]), simplify = TRUE)
			wtw <- wtw + Weights[s]*crossprod(w1, w2)
			}
		D2 <- (wtw + t(wtw) + Weights[1]*crossprod(sapply(Z.list, function(Z) c(Z*resid), simplify = TRUE)))/(nr*nc)
		asy.var <- inv.D0%*%D2%*%inv.D0/(nr*nc)
		}
#6 # "c&s-heteros.T/N^2,N/T^2-->0"
	if(Error.typ == 6){
		Z.ti.sig2 <- sapply(Z.list, function(X) c(resid*X), simplify = TRUE)
		D3 <- crossprod(Z.ti.sig2)/(nr*nc)
		asy.var <- inv.D0%*%D3%*%inv.D0/(nr*nc)
		}
#7 # "c&s-heteros.T/N-->const"
	if(Error.typ == 7){
	  if(d > 0){
	  # cross section bias
		ind.sig2.hat <- diag(sqrt(diag(crossprod(resid))/nr)) # (nc x nc)
		inv.crossdiagA <- diag((crossdiagA)^{-1}, d)
		FUN1B <- function(X) sum(diag(t(X%*%M.A)%*%F%*%inv.crossdiagA%*%t(A)*ind.sig2.hat))/nr
		BiasBPar2 <- as.matrix(sapply(X.mat.list, FUN1B, simplify = TRUE))
		BiasB <- -inv.D0%*%BiasBPar2
	  # time serial bias
		tim.sig2.hat <- diag(sqrt(diag(tcrossprod(resid)))/nc) # (nr x nr)
		FUN2C <- function(X) sum(diag(t(X)%*%M.F%*%tim.sig2.hat%*%F%*%inv.crossdiagA%*%t(A)*diag(1, nc)))/nr
		BiasCpar2 <- as.matrix(sapply(X.mat.list, FUN2C, simplify = TRUE))
		BiasC <- -inv.D0%*%BiasCpar2

	  # bias corrected estimator 
		slope <- slope - BiasB/nc - BiasC/nr
	  }
	  # calculate D3
		Z.ti.sig2 <- sapply(Z.list, function(X) c(resid*X), simplify = TRUE)
		D3 <- crossprod(Z.ti.sig2)/(nr*nc)
		asy.var <- inv.D0%*%D3%*%inv.D0/(nr*nc)
		}

#8 #"c&s-corr.T/N --> const"
	if(Error.typ == 8){ 
	  # cross section bias
	  if(d > 0){
		inv.crossdiagA <- diag((crossdiagA)^{-1}, d)
		FUN1B <- function(X) sum((t(X%*%M.A)%*%F%*%inv.crossdiagA%*%t(A)%*%crossprod(resid)/nr)[1:round(log(nc)), 1:round(log(nc))])*nc/(nr*round(log(nc)))
		BiasBPar2 <- as.matrix(sapply(X.mat.list, FUN1B, simplify = TRUE))
		BiasB <- -inv.D0%*%BiasBPar2
	  # time serial bias
		#tim.sig2.hat <- diag(sqrt(diag(tcrossprod(resid)))/nc) # (nr x nr)
		wtw <- matrix(0, P, 1)
		for(s in 2:Ltrank){#s=2
			w1 <- sapply(X.mat.list, function(X) c((M.F%*%X*resid)[1:(nr - s +1) , ]), simplify = TRUE)
			w11 <- sapply(X.mat.list, function(X) c((M.F%*%X*resid)[s:nr , ]), simplify = TRUE)
			w2 <- matrix(F[s:nr , ]%*%inv.crossdiagA%*%t(A), ncol = 1)
			w21 <- matrix(F[1:(nr - s +1), ]%*%inv.crossdiagA%*%t(A), ncol = 1)
			wtw <- wtw + Weights[s]*crossprod(w1, w2) + Weights[s]*crossprod(w11, w21) 
			}
		XMOFAAA <- (wtw  + Weights[1]*crossprod(sapply(X.mat.list, function(X) c((M.F%*%X*resid)), simplify = TRUE), matrix(F%*%inv.crossdiagA%*%t(A), ncol = 1)))/(nr)
		BiasC <- -inv.D0%*%XMOFAAA

	  # bias corrected estimator 
		slope <- slope - BiasB/nc - BiasC/nr
	  }
	  # calculate DZ
		Z.icor.sig2 <- sapply(Z.list, function(X) rowSums(resid[, 1:round(log(nc))]*X[, 1:round(log(nc))]), simplify = TRUE)
		wtw <- matrix(0, P, P)
		for(s in 2:Ltrank){#s=2
			w1 <- Z.icor.sig2[1:(nr - s +1) , ]
			w2 <- Z.icor.sig2[s:nr , ]
			wtw <- wtw + Weights[s]*crossprod(w1, w2)
		}
		DZ <- (wtw + t(wtw) + Weights[1]*crossprod(Z.icor.sig2))/(nr*round(log(nc)))
		asy.var <- inv.D0%*%DZ%*%inv.D0/(nr*nc)
		}

	mpp <- sqrt(diag(asy.var))  
	remeind <- matrix(y - x%*%slope, nr, nc)
  
  ## Add intercept if it exists in the formula
 #if(is.intercept){
  #rownames(Intercept) <- "(Intercept)"
  #colnames(Intercept) <- " "
  #slope <- rbind(Intercept, slope)
  #asy.var.inter <- sig2.hat/(nc*nr) + 
  #        matrix(OvMeans[-1], 1, P, byrow = TRUE)%*%asy.var%*%matrix(OvMeans[-1], P, 1)
  #mpp11 <-  sqrt(c(asy.var.inter))
  #mpp <- c(mpp11, mpp)
  #}
  
	test <- slope/mpp 
	p.value <- (1 - pnorm(abs(test)))*2
	inf.result <- cbind(slope,
                            mpp,
                            test,
                            p.value)
	colnames(inf.result) <- c("Estimate", "Std.Err", "Z value", "Pr(>z)")
	#if(Error.typ == 7 |  Error.typ == 8) colnames(inf.result)[1] <- c("Bias Correct.Estimate")
	result <- list(ZZ = D0, inv.ZZ = inv.D0, inf.result=inf.result, sig2.hat=sig2.hat, BC.remeind = remeind) 	
  return(result)
}

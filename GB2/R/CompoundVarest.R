
# Computes the scores from the Gamma factors
scoreU.cgb2 <- function (fac, pl) 
{
    L <- length(pl)
    denom <- fac %*% pl
    num <- fac[, -L]
    midt <- num/as.vector(denom) - 1
    U <- t(pl[-L]*t(midt))
    return(U)
}

varscore.cgb2 <- function(U,w=rep(1,dim(U)[1]))
{
    if (dim(U)[1] != length(w) ){
		warning("error in dimensions: length of w= ",length(w),"; dim(U)= ",dim(U)[1])
		return()
	}
    Vsc <- t(U*w) %*% (U*w)
    return(Vsc)
}

desvar.cgb2 <- function(data=data, U=U, ids=NULL, probs=NULL, strata = NULL, variables = NULL, fpc=NULL,
	nest = FALSE, check.strata = !nest, weights=NULL, pps=FALSE, variance=c("HT","YG")) {
	datfull <- cbind(data,U)
	Names <- colnames(U)
	formul <-  as.formula(paste(" ~ ", paste(Names, collapse= "+")))
	dstr <- svydesign(data = datfull, ids=ids, probs=probs, strata = strata, variables = variables, fpc=fpc,
	nest = nest, check.strata = check.strata, weights=weights,pps=pps,variance=variance)
	v <-svytotal(formul, dstr, cov=TRUE)
	Vtheta <- vcov(v)
	return(list(svytotal=v,Vtheta=Vtheta))
}

hess.cgb2 <- function(U,pl,w=rep(1,dim(U)[1])){
	L <- length(pl) # 1 x L  	pl: vector of mixture probabilities
	dU <- dim(U)	 # n x (L-1) 	U: matrix of scores, output of scoreU.cgb2 (see eq. 15)
	L1 <- L-1
	Lw <- length(w)
	if  ((dU[2]!=L1)|(dU[1]!=Lw)){
		warning("error in dimensions: no of parameters= ",L1,"; length of w= ",Lw,"; dim(U)= ",dU[1],",",dU[2])
		return()
	}
	else{
		V1 <- -t(U)%*%(U*w)
		sumsc <- colSums(U*w) 
#		V2 <- V1-pl[-L]%*%t(sumsc) - sumsc%*%t(pl[-L]) + diag(sumsc)
		V2 <- V1-pl[-L]%*%t(sumsc) - sumsc%*%t(pl[-L])			# change 2014-04-21
		V2 <- V2 + ifelse(length(sumsc)==1, sumsc, diag(sumsc))	# change 2014-04-21
		colnames(V1)<- paste("v",1:L1,sep="")
		rownames(V1)<- paste("v",1:L1,sep="")
		dimnames(V2) <- dimnames(V1)
	eigv <- eigen(V2)[[1]]
	return(V2)
	}
}

vepar.cgb2 <- function(ml,Vsc, hess)
{
    estimate <- ml[[2]]$par
    V <- solve(hess)
    Vcov <- V %*% Vsc %*% V
    stderr <- sqrt(diag(Vcov))
    Vcor <- diag(1/stderr)%*%Vcov%*%diag(1/stderr)
    names(estimate) <- rownames(Vcov)
    dimnames(Vcor) <- dimnames(Vcov)
    return(list(type="parameter",estimate=estimate, stderr=stderr, Vcov=Vcov, Vcor=Vcor))
    
}

derivind.cgb2 <- function (shape1, scale, shape2, shape3, pl0, pl, prop=0.6, decomp="r") {
    par <-vofp.cgb2(pl)
    indic <- function(par) {
	pl <- pofv.cgb2(par)
	return(main.cgb2(prop, shape1, scale, shape2, shape3,pl0,pl,decomp=decomp))
    }
    estimate <- t(indic(par))
    names(estimate) <- c("median","mean","arpr","rmpg","qsr")
    rownames(estimate) <- ""
 
    MFDI <- jacobian(indic, par, method = "Richardson", method.args = list())
    rownames(MFDI) <- c("median","mean","arpr","rmpg","qsr")
    colnames(MFDI) <- paste("v",1:length(par),sep="")

    return(list(estimate=estimate, jacobian = MFDI))
}
veind.cgb2 <- function(Vpar,shape1, scale, shape2, shape3, pl0, pl, decomp="r") {
    esti <- derivind.cgb2(shape1, scale, shape2, shape3, pl0, pl, decomp=decomp) 
    MFDI <- esti[["jacobian"]]
    Vcov <- MFDI%*%Vpar[["Vcov"]]%*%t(MFDI)
    std <- sqrt(diag(Vcov))
    Vcor <- diag(1/std)%*%Vcov%*%diag(1/std)
    rownames(Vcor) <- colnames(Vcor) <- rownames(MFDI)
    return(list(type="indicator",estimate=esti[["estimate"]],stderr=std,Vcov=Vcov,Vcor=Vcor))
}
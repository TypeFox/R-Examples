# nfactor number of factor covariate
# ind.place: indicator of first factor covariate modalite in global covariate 
#N: number of covariates
# v.modal: a vector of number of parameters estimates by factor covariates
#Beta, vector of beta parameters estimated
#VarBeta, vector of variance beta parameters estimated
#modality vector of ddl

waldtest <- function(N,nfact,place,modality,b,Varb,Lfirts=NULL,Llast=NULL,Ntot=NULL){
	waldValue <- rep(0,nfact)
	for(i in 1:nfact){
   	if(place[i] == 1){
     taillecol <- N-modality[i]
			Lapres <- matrix(0,nrow=modality[i],ncol=taillecol)
			Lvar <- diag(rep(1,modality[i]))
			L <- cbind(Lvar,Lapres)
			if(!is.null(Lfirts)) L <- cbind(matrix(0,nrow=modality[i],ncol=Lfirts),L)
			if(!is.null(Llast)) L <- cbind(L,matrix(0,nrow=modality[i],ncol=Llast))
		}
		if(place[i] == 2){
      taillecol <- N-modality[i]
			Lapres <- matrix(0,nrow=modality[i],ncol=taillecol-1)
			Lvar <- diag(rep(1,modality[i]))
			L <- cbind(cbind(rep(0,modality[i]),Lvar),Lapres)
			if(!is.null(Lfirts)) L <- cbind(matrix(0,nrow=modality[i],ncol=Lfirts),L)
			if(!is.null(Llast)) L <- cbind(L,matrix(0,nrow=modality[i],ncol=Llast))
		}
		if(place[i] > 2){
    	taillecol <- place[i]-1 
    	Lavant <- matrix(0,nrow=modality[i],ncol=taillecol)	
			taillecol <- N-modality[i]-taillecol
			Lapres <- matrix(0,nrow=modality[i],ncol=taillecol)
      Lvar <- diag(rep(1,modality[i]))

      L <- cbind(Lavant,cbind(Lvar,Lapres))
			if(!is.null(Lfirts)) L <- cbind(matrix(0,nrow=modality[i],ncol=Lfirts),L)
			if(!is.null(Llast)) L <- cbind(L,matrix(0,nrow=modality[i],ncol=Llast))
		}
		if(is.null(Ntot))Ntot <- N
		out <- .Fortran("waldmultiv",
				as.double(b),
				as.integer(Ntot),
				as.integer(L),
				as.double(Varb),
				as.integer(modality[i]),
				Wald=as.double(0),
				PACKAGE = "frailtypack")
		waldValue[i] <- out$Wald
#		print(L)
	}
	return(waldValue)
}# end Wald

#################################################################
#######                    Fieldsim                      ########
#################################################################

## fieldsim.R  (2006-16)
##
##    Random field simulation by the FieldSim method
##
## Copyright 2006-16 Alexandre Brouste and Sophie Lambert-Lacroix


##    INPUT VARIABLES
#################################################################
##	process     : process (R process object)
##  Ne          : number of points to be simulated with the exact
##                method
##  nbNeighbor  : number of neighbors to be considered in the 
##                fieldsim method
#################################################################

##    OUTPUT VARIABLES
#################################################################
## function returns the value of the process in the corresponding
## slot of the process object with the FieldSim method
#################################################################

fieldsim <-function(process,Ne,nbNeighbor=4){
	
#-----------------Tests on input variables-----------------------
		
if(missing(process)){ 		
	cat("Error from fieldsim.R: parameter process is missing\n")
	return(NULL)
}	
	
if(!isS4(process)){ 
	cat("Error from fieldsim.R: parameter process is not of type process\n")
	return(NULL)
}else if(!(class(process)[1]=="process")){
	cat("Error from fieldsim.R: parameter process is not of type process\n")
	return(NULL)
}

manifold<-process@manifold
R<-process@covf	
parameter<-process@parameter
	
mesh<-manifold@atlas	
dimen<-dim(mesh)[1]
		
#Traitement des NaN
		
meshnan<-is.nan(mesh)
indexNan<-meshnan[1,]
	
for (i in 1:dimen){ 
	indexNan<-indexNan|meshnan[i,]
}
		
meshtmp<-mesh
mesh<-matrix(mesh[1:dimen,!indexNan],nrow=dimen)

Ntmp<-dim(meshtmp)[2]
N<-dim(mesh)[2]	
	

	
if(missing(Ne)){
		Ne<-N
}

if(!is.numeric(Ne)){
	cat("Error from fieldsim.R: Ne must be of numeric type\n")
	return(NULL)

}
	
if (Ne<=0|Ne>N){
	cat("Error from fieldsim.R: Ne have must be up to 0 and less than the size of the atlas\n")
	return(NULL)
}	

if(!is.numeric(nbNeighbor)){
	cat("Error from fieldsim.R: nbNeighbor must be of numeric type\n")
	return(NULL)
		
}
	
if(Ne<=nbNeighbor){
	cat("Error from fieldsim.R: Ne must be bigger than nbNeighbor\n")
	return(NULL)
}	

if((nbNeighbor<2)|(nbNeighbor>32)){
	cat("Error from fieldsim.R: nbNeighbor must belong to {1,...,32}\n")
	return(NULL)
}



	matcov<-matrix(0,Ne,Ne)
		
	for(i in 1:Ne){
		for (j in i:Ne){
			matcov[i, j]<- R(mesh[,i],mesh[,j])		
			matcov[j, i]<- matcov[i, j]
		}
	}

	eps=1e-10
	iso<-apply((matcov<(0+eps))&(matcov>(0-eps)),2,all)
	whereo<-which(iso)	

	if (any(iso)){
		matcovbis<-matcov[-whereo,-whereo]
		Nbis<-dim(matcovbis)[2]
	}else{
		matcovbis<-matcov
		Nbis<-Ne
	}
    
	L <-matrix(0,Nbis,Nbis)	
	L <- chol(matcovbis)
	Z <- rnorm(Nbis)
	restmp <- t(L) %*% Z
	
	if(any(iso)){
		res<-rep(0,length=Ne)
		res[-whereo]<-as.vector(restmp)
	}else{
		res<-as.vector(restmp)
	}	
	
	
	if (Ne<N){
		
		d<-manifold@distance
		dmat=nbNeighbor+1
				
		for (ts in (Ne+1):N){
			D<-0
			for (k in 1:(ts-1)){
				D[k]<- d(mesh[,ts],mesh[,k])
			}
			
			Voisins<-sort(D,index=TRUE)$ix[1:nbNeighbor]
			xx<-matrix(mesh[1:dimen, c(Voisins,ts)],nrow=dimen)	
					
			Rmat = matrix(0,dmat,dmat)

			for (i in 1:dmat){
				for (j in 1:dmat){ 
						Rmat[i,j]<- R(xx[,i],xx[,j])
				}
			}
			
			ismato<-apply((Rmat<(0+eps))&(Rmat>(0-eps)),2,all) 
			wheremato<-which(ismato)	
			
			if (any(wheremato==dmat)){
				xsim<-0
			}
			else{
			
				if (any(ismato)){
					Rmatbis<-Rmat[-wheremato,-wheremato]
					xvec <- res[Voisins[-wheremato]]
					dmatbis<-dim(Rmatbis)[2]
				}else{
					Rmatbis<-Rmat
					xvec <- res[Voisins]
					dmatbis<-dmat
				}
			
			
				Li=t(chol(Rmatbis))  
				D=diag(diag(Li))
				Ld<-dim(D)[2]	
				Lin=Li%*%solve(D)
				Linmoins1 = solve(Lin)
					
				xsim<-D[dmatbis,dmatbis]*rnorm(1)-sum(Linmoins1[dmatbis,1:(dmatbis-1)]*xvec)
			}
			
			res=c(res,xsim)	
			}
		}else{}
	
	
#restmp<-rep(NaN,Ntmp)
#restmp[!indexNan]<-res
#res<-restmp
	
	
    if (process@name=="bridge"){
		
		Gamma<-parameter$Gamma
		M<-dim(Gamma)[2]
		Rbis<-parameter$R
		Tp<-matrix(parameter$Tp,M,M)
		
		#x<-meshtmp
		x<-mesh
		
		taille<-dim(Gamma)[1]
		Gammavalue<-matrix(Gamma[taille,],M,1)
        
        
        Sigmatmp<- Tp %*% Gammavalue
        
		if (is.null(Gammavalue)){}
		else{
            
            for (l in 1:N){
                
                Qtmp<-matrix(0,1,M)
                for (j in 1:M){
                    Qtmp[1,j]<-Rbis(x[,l],Gamma[1:(taille-1),j])
                }
                
                res[l]<-res[l]+Qtmp%*%Sigmatmp
                
            }}}
	
	restmp<-rep(NaN,Ntmp)
	restmp[!indexNan]<-res
	res<-restmp
    
	
	nameProcess<-deparse(substitute(process))
	process@values<-res
	assign (nameProcess,process,envir = parent.frame())
	return(invisible(1))
}


#################################################################
#######                        Fin                         ######
#################################################################
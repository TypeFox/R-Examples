CreateQuadrature <- function(N,L,QuadPoly,ExpPoly,QuadType,ParamDistrib){

#	N <- 3  # number of random variables
# P <- 2  
#	L <- c(4,4) # quadrature levels in each dimention
#	M <- getM(N,P) # number of PCE termrs
#	ParamDistrib<- list(beta=c(0.0,0.0),alpha=c(0.0,0.0))	
#	Index = indexCardinal(N,P) # generate the basis of expansion
#	QuadType <- c("SPARSE") # type of quadrature
  M <- NA
	if (length(L) != N ){
		if (length(L) == 1 ){
			cat('Quadrature level assumed to be isotropic.\n')
      L = L*rep(1,N)
		} else {  
			cat('ERROR: Quadrature level not the same as number of random dimensiones\n')
		}	
	}
  
  # if (!is.null(ParamDistrib)){print(ParamDistrib)}
  
	if (QuadType == 'FULL'){
	  # cat("In Full loop \n")
    # determine the number of quadrature points
		Z = prod(L)
		
		# determine correct permulation for indices
		tmp = permute(1,1,N,L,Z,rep(0,N))
		NDIi = tmp$listNDIi
	  # cat("Permuted \n")
	  
		# populate the quadrature points and weights
		res = populateRootsWeights(N,L,Z,QuadPoly,ParamDistrib,NDIi)
	 	NDXi = res$listNDXi 
		NDWi = res$listNDWi
	  # cat('Past popRootsWeights \n')
		# populate matrix of phi_m(X_z), for m = 0 to M qnd z = 0 to Z-1
		# cat('L=',L)
		Index = indexCardinal(N,max(L)-1) 
		M = getM(N,max(L)-1)
		PolyNodes = populateNDPhZn(N,M,Z,ExpPoly,NDXi,Index,ParamDistrib)
		# cat('End first loop \n')
		
	} else if (QuadType == 'SPARSE'){
		L <- L[1]
		M = getM(N,L-1)
		
		# Checks automatically if 
		# ExpOdd + Fejer => nested Fejer used
		# ExpOdd + ClenshawCurtis => nested Clenshaw-Curtis used
		# *KP* => nested Gauss grid used
		Growth <- 'ExpOdd' # only nominal increase type coded 
		
		## Define the distribution of the 1D nodes    
		NominalSize <- GetNominalSize(L,Growth,QuadPoly) # 'Growth' is hardcoded as as there are currently no other implementation 
		## Define the representation of the 1D nodes
		NominalList <- GetNominalList(L,NominalSize,QuadPoly)
  	
		# Generate multivariate sparse grid according the Smolyak rules
		SparseGrid <- GenerateSparseGridSymmetry(N,L,NominalList,NominalSize,ExpPoly,ParamDistrib)
		Z <- SparseGrid$Size[L]

		# populate the quadrature points and weights
		NDXi <- t(SparseGrid$Node[1:SparseGrid$Size[L],1:N])
		NDWi <- SparseGrid$Weight[L,1:SparseGrid$Size[L]]

		# populate matrix of phi_m(X_z), for m = 0 to M qnd z = 0 to Z-1
		PolyNodes <- SparseGrid$PolyNodes
	}
	res = list(QuadSize = Z, QuadNodes = NDXi, QuadWeights = NDWi, PolyNodes = PolyNodes)
  # cat('Compile results \n')
  
  return(res)
}
#####################################################
## COVSTATIS processing
#####################################################

mpCOVSTATIS.core <- function(data, normalization = 'None', masses = NULL, table = NULL, make.table.nominal=TRUE)	
{ print('Preparing data')
  
  data = as.matrix(data)
  n.rows = dim(data)[1]
  n.cols = dim(data)[2]
  n.unique = length(unique(as.vector(data)))
	
   table <- mpTableCheck(data,table,make_table_nominal=TRUE)
   n.groups = dim(data)[2]/dim(data)[1]
   rawData = array(0,dim=c(rowSums(table)[1],rowSums(table)[1],dim(table)[1]))  
   for(i in 1:n.groups)
   {	to = i * n.rows
  	 	from = to - (n.rows-1)
  	 	rawData[,,i] <- data[,from:to]
    } 
     
  rawData.mat = c()
  for(i in 1:dim(rawData)[3])
  {   rawData.mat <- matrix(c(rawData.mat,rawData[,,i]),n.rows) 
  }
  
# masses
  if(is.null(masses)== TRUE)
   {   masses = matrix(1/dim(rawData)[1],dim(rawData)[1],1)
   }
		
  else if (is.null(masses)== FALSE)
  {    masses = matrix(diag(dim(rawData)[1]) * masses,dim(rawData)[1],1)
  }

# scalar product matrices
  scalarProductMatrices = array(0,dim=c(n.rows,n.rows,n.groups))
  phi = array(0,dim=c(n.rows,n.rows,n.groups))
  for(i in 1:n.groups)
  {   phi[,,i] = diag(dim(rawData[,,i])[1])-(matrix(1,dim(rawData[,,i])[1],1) %*% t(masses))
      scalarProductMatrices[,,i] = ((0.5) * phi[,,i] %*% rawData[,,i] %*% t(phi[,,i])) 
  }


# normalization
  norm = array(0,dim=c(n.rows,n.rows,n.groups))
  for(i in 1:n.groups)
  {  if(normalization != 'None' && normalization != 'MFA' && normalization != 'SumPCA' )
      {	  print(paste('WARNING: Normaliztion option not recognized. MFA was set as default'))
          norm[,,i]=scalarProductMatrices[,,i]/corePCA(scalarProductMatrices[,,i])$pdq$Dv[1]
      }
      else if(normalization == 'None')
      {	  norm[,,i] = scalarProductMatrices[,,i] 
      }
      else if(normalization == 'MFA')
      {	  norm[,,i]=scalarProductMatrices[,,i]/svd(scalarProductMatrices[,,i])$d[1]
      }
      else if(normalization == 'SumPCA')
      {	  norm[,,i]=scalarProductMatrices[,,i]/(sqrt(sum(scalarProductMatrices[,,i]*scalarProductMatrices[,,i])))
      }
    }

  norm2 <- array(norm, dim=c(n.rows*n.rows,n.groups))

#C Matrix
  CMatrix = t(norm2) %*% norm2

# RV Matrix
  Norm = sqrt(apply(norm2^2,2,sum))
  rvMatrix = CMatrix/(t(t(Norm)) %*% Norm)

  # # CMatrix
  # CMatrix = diag(1,n.groups)
  # for(i in 1:(n.groups))
  # {  for(j in i:n.groups)
  #    {   rv = rvCoeff(norm[,,i],norm[,,j],type=0)
  #        CMatrix[i,j] = rv
  #        CMatrix[j,i] = rv
  #     }
  #  }
 		
# eigen decomposition
   C.decomp = corePCA(CMatrix)
   decomp.C = C.decomp$pdq
 
# contribution
   ci = C.decomp$ci
   	if(is.null(rownames(table))){
	   rownames(ci) <- paste("Table", 1:dim(table)[1], sep = "")
	   table.names <- rownames(ci)
	}else{
	  	rownames(ci) <- rownames(table)
	  	table.names <- rownames(ci)
	}

# contribution
	cj = C.decomp$cj
	rownames(cj) <- rownames(ci)
    
# eigen vectors
   P = decomp.C$p
   rownames(P) <- rownames(ci)

# eigen values
   D= (decomp.C$Dv)
	
# cumulative eigen values
   D.cum = cumsum(D)

# factor scores
   G = decomp.C$p %*%  sqrt(decomp.C$Dd)
	rownames(G) <- rownames(P)

# percent of variance explained
   taus <- decomp.C$Dv/sum(decomp.C$Dv) * 100

# Alpha Weights
  alphaWeights = P[,1] / sum(P[,1])

##########################################
# Compromise
##########################################
	
# Compromise (S+)
  compromiseMatrix = matrix(0,n.rows, n.rows)
  for(i in 1:n.groups)
  {	compromiseMatrix = compromiseMatrix + alphaWeights[i] * norm[,,i]
  }
  
# analyzing the compromise
   PCA.compromise <- corePCA(compromiseMatrix)
   compromise.PCA <- PCA.compromise$pdq

# contribution
	compromise.ci <- PCA.compromise$ci
	rownames(compromise.ci)=rownames(data)

# contribution
	compromise.cj <- PCA.compromise$cj
	rownames(compromise.ci)=rownames(data)

# eigen vectors
   compromise.P = compromise.PCA$p
   rownames(compromise.P)=rownames(data)

# eigen values
   compromise.dd = (compromise.PCA$Dv)

# factor scores
   compromise.G = compromise.PCA$p %*% diag(sqrt(compromise.PCA$Dv))
   rownames(compromise.G)=rownames(data)
 	
# % of variance explained
   compromise.tau <- compromise.PCA$Dv/sum(compromise.PCA$Dv) * 100

##########################################
# Tables: Generalized PCA of X
##########################################	

# column names	
   table.colnames <- colnames(table)

# alpha weights
   table.alphaWeights <- alphaWeights

# weights and masses
   M =  rep(1/(dim(data)[1]),dim(data)[1])
   	
   w = c()
   for(i in 1:length(rowSums(table)))
   { w = c(w, rep(alphaWeights[i],rowSums(table)[i]))
   }

#general PDQ
	pdq.general = corePCA(rawData.mat,M=M,W=w)
	general.pdq = pdq.general$pdq

# contribution
	table.ci = pdq.general$ci

# contribution
	table.cj = pdq.general$cj
   
# Eigen vectors of the tables
   gpdq.vectors = general.pdq$p

# Eigen values of the tables 
   gpdq.eigenvalues = (general.pdq$Dd)^2
	
# Inertia
   gpdq.inertia = ((general.pdq$Dv) / sum(general.pdq$Dv))*100

# Loadings of the tables
   gpdq.loadings = general.pdq$q
   rownames(gpdq.loadings) = colnames(data)

# Factor scores of the tables
   gpdq.factorscores =  general.pdq$p %*%  (general.pdq$Dd)
   rownames(gpdq.factorscores)=rownames(data)
   
# Partial Factor Scores
   gpdq.partial = array(0,dim=c(dim(data)[1],dim(gpdq.loadings)[2],dim(table)[1]))
   to_partial = 0
   from_partial = 1
   for(i in 1:dim(table)[1])
   {   from = sum(table[i-1,]) + from_partial
       to = sum(table[i,]) + to_partial
       to_partial = to
       from_partial = from
       gpdq.partial[,,i] = rawData.mat[,from:to] %*% gpdq.loadings[from:to,]
    
   }
   
   gpdq.partialFS <- matrix(0,dim(data)[1]*dim(table)[1],dim(gpdq.loadings)[2])
   to.total = 0
   for(i in 1:dim(table)[1])
   {	from = to.total + 1
   		to = i*dim(gpdq.partial)[1]
   		to.total = to
   		gpdq.partialFS[from:to,]= gpdq.partial[,,i]
   	}
   	rownames(gpdq.partialFS) = paste(rep(table.names,each=dim(data)[1]),rep(rownames(data)))


######################################
# Results
######################################

res.distatis <- list(data = data, normalization= normalization, table = table, S=scalarProductMatrices, C = CMatrix, ci = ci, cj = cj, 
                eigs.vector= P, eigs = D, fi = G, tau = taus, alphaWeights = alphaWeights, rvMatrix = rvMatrix, 
			
                compromise = compromiseMatrix, compromise.ci = compromise.ci, compromise.cj = compromise.cj, compromise.eigs.vector = compromise.P, 
                compromise.eigs = compromise.dd, compromise.fi = compromise.G, compromise.tau = compromise.tau,

                 masses = M, table.partial.fi.array = gpdq.partial, table.cj = compromise.cj, table.ci = compromise.ci, 
			           table.eigs = compromise.dd, table.tau= compromise.tau, table.eigs.vector = compromise.P, 
			           table.fi = compromise.G, table.partial.fi= gpdq.partialFS, table.Q = gpdq.loadings)

return (res.distatis)
}





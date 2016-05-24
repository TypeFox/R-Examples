mpCANOSTATIS.core <- function (data, num.obs = num.obs, column.design, row.design , num.groups = num.groups, normalization = 'MFA', masses = NULL)
{	
	# Mahalanobis distances
	D <- array(0,dim=c(dim(row.design)[2],dim(row.design)[2],num.groups))
	for(i in 1:dim(column.design)[1])
	{	D[,,i] <- mpMahalanobis(data[,c(which(column.design[i,]==1))],row.design) 
	}
	
	D.mat = c()
  for(i in 1:dim(D)[3])
 	{ D.mat <- matrix(c(D.mat,D[,,i]),dim(row.design)[2]) 
  }
  	
  # masses
  if(is.null(masses)== TRUE)
  { masses = matrix(1/dim(D)[1],dim(D)[1],1)
  }
	else if (is.null(masses)== FALSE)
  { masses = matrix(diag(dim(D)[1]) * masses,dim(D)[1],1)
  }

  # scalar product matrices
  scalarProductMatrices = array(0,dim=c(dim(row.design)[2],dim(row.design)[2],num.groups))
  phi = array(0,dim=c(dim(row.design)[2],dim(row.design)[2],num.groups))
  for(i in 1:num.groups)
  { phi[,,i] = diag(dim(D[,,i])[1])-(matrix(1,dim(D[,,i])[1],1) %*% t(masses))
    scalarProductMatrices[,,i] = (-(0.5) * phi[,,i] %*% D[,,i] %*% t(phi[,,i])) 
  }
  
  # normalization
  norm = array(0,dim=c(dim(row.design)[2],dim(row.design)[2],num.groups))
  for(i in 1:num.groups)
  {   if(normalization != 'None' && normalization != 'MFA' && normalization != 'SumPCA' && normalization != '1Norm')
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

	norm2 <- array(norm, dim=c(dim(row.design)[2]*dim(row.design)[2],num.groups))

# CMatrix
	CMatrix = t(norm2) %*% norm2

# RV Matrix
  Norm = sqrt(apply(norm2^2,2,sum))
  rvMatrix = CMatrix/(t(t(Norm))%*%Norm)
		
# eigen decomposition
   C.decomp = corePCA(CMatrix)
   decomp.C = C.decomp$pdq
 
# contribution
   ci = C.decomp$ci

   rownames(ci) <- paste("Table", 1:num.groups, sep = "")
   table.names <- rownames(ci)

# contribution
	cj = C.decomp$cj
	rownames(cj) <- rownames(ci)
    
# eigen vectors
   P = decomp.C$p
   rownames(P) <- rownames(ci)

# eigen values
   D= (decomp.C$Dv)

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
  compromiseMatrix = apply(apply(norm,c(1,2),'*',t(alphaWeights)),c(2,3),sum)

# analyzing the compromise
   PCA.compromise <- corePCA(compromiseMatrix)
   compromise.PCA <- PCA.compromise$pdq

# contribution
	compromise.ci <- PCA.compromise$ci
	rownames(compromise.ci) <- colnames(row.design)

# contribution
	compromise.cj <- PCA.compromise$cj
	rownames(compromise.ci)=rownames(compromise.ci)

# eigen vectors
   compromise.P = compromise.PCA$p
   rownames(compromise.P)=rownames(compromise.ci)

# eigen values
   compromise.dd = (compromise.PCA$Dv)


# factor scores
   compromise.G = compromise.PCA$p %*% diag(sqrt(compromise.PCA$Dv))
   rownames(compromise.G)=rownames(compromise.ci)
   
# how each assessors interprets the space
	f <- compromiseMatrix %*% (compromise.PCA$p %*% diag(compromise.PCA$Dd^(-0.5)))
	
# % of variance explained
   compromise.taus <- compromise.PCA$Dv/sum(compromise.PCA$Dv) * 100

##########################################
# Tables: Generalized PCA of X
##########################################	

# alpha weights
   table.alphaWeights <- alphaWeights

# weights and masses 	
     M = diag(1/(dim(D.mat)[1]),dim(D.mat)[1],dim(D.mat)[1])    
   w = rep(alphaWeights,dim(row.design)[2])
   W = diag(w)

#general PDQ
	pdq.general = corePCA(D.mat,M=M,W=W)
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
   rownames(gpdq.loadings) <-  paste(rep(paste('Table',1:dim(column.design)[1],sep="."),each=dim(row.design)[2]),
                                rep(colnames(row.design),dim(column.design)[1]))

# Factor scores of the tables
   gpdq.factorscores =  general.pdq$p %*%  (general.pdq$Dd)
   rownames(gpdq.factorscores)=paste("Group", 1:dim(row.design)[2], sep = "")
   
# Partial Factor Scores
   gpdq.partial = array(0,dim=c(dim(row.design)[2],dim(row.design)[2],dim(column.design)[1]))
   for(i in 1:dim(column.design)[1])
   {   from = dim(row.design)[2]*(i-1) + 1
       to = dim(row.design)[2]*i
       gpdq.partial[,,i] = D.mat[,from:to] %*% gpdq.loadings[from:to,]
    
   }
   
   gpdq.partialFS <- array(gpdq.partial,dim=c(dim(row.design)[2]*dim(row.design)[2],dim(column.design)[1]))
	
	rownames(gpdq.partialFS)=(paste(rep(paste("Table",1:dim(row.design)[2],sep=""),each=dim(row.design)[2]),rep(colnames(row.design))))

######################################
# Results
######################################

res.canostatis <- list(mahalanobis=D, normalization= norm, S = scalarProductMatrices, column.design = column.design, row.design = row.design, 
      S=scalarProductMatrices, C = CMatrix, rvMatrix = rvMatrix, ci = ci, cj = cj, eigs.vectors = P, eigs = D, fi = G, tau = taus, alphaWeights = alphaWeights,
			
			compromise = compromiseMatrix, compromise.ci = compromise.ci, compromise.cj = compromise.cj, compromise.eigs.vector = compromise.P, compromise.eigs = compromise.dd, 
      compromise.fi = compromise.G, compromise.tau = compromise.taus,
			
			masses = M, table.partial.fi.array = gpdq.partial, table.cj = compromise.cj, table.ci = compromise.ci, table.eigs = compromise.dd, table.tau = compromise.taus, 
      table.eigs.vector = compromise.P, table.fi = compromise.G, table.partial.fi = gpdq.partialFS, table.Q = gpdq.loadings)

return (res.canostatis)
}
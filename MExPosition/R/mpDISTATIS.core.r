#####################################################
## DISTATIS processing
#####################################################

mpDISTATIS.core <- function(data, table, sorting = 'No', normalization = 'None', masses = NULL, make.table.nominal=TRUE)	
{ 
  print('Preparing data')
  data = as.matrix(data)
  n.rows = dim(data)[1]
  n.cols = dim(data)[2]
  n.unique = length(unique(as.vector(data)))

  #Sorting Task?
  if(sorting == 'No')
  { table <- mpTableCheck(data,table,make_table_nominal=make.table.nominal)
    n.groups = dim(data)[2]/dim(data)[1]
    D = array(0,dim=c(rowSums(table)[1],rowSums(table)[1],dim(table)[1]))
  
    for(i in 1:n.groups)
    {   to = i * n.rows
  	    from = to - (n.rows-1)
  	    D[,,i] <- data[,from:to]
    }
  }
 
  if(sorting == 'Yes')
  { print('Sorting Task')
  	table <- t(makeNominalData(table))
    colnames(table) = paste(rep(paste('Assessor',1:dim(data)[2],sep="."),each=dim(data)[1]),rep(rownames(data),dim(data)[1]))
    n.groups = dim(data)[2]
    L = array(0,dim=c(n.rows,n.unique,n.cols))
    for (i in 1:n.cols)
    {	   L[,1:length(unique(as.vector(data[,i]))),i] = makeNominalData(as.matrix(data[,i]))
    }

    R=array(0,dim=c(n.rows,n.rows,n.cols))
    for (i in 1:n.cols)
    {	   R[,,i] = L[,,i] %*% t(L[,,i])
    }
      
    D=array(0,dim=c(n.rows,n.rows,n.cols))

    for(i in 1:n.cols)
    { D[,,i]=matrix(1,n.rows,n.rows)-R[,,i]
    }	
  }	 

  D.mat = c()
  for(i in 1:dim(D)[3])
  {   D.mat <- matrix(c(D.mat,D[,,i]),n.rows) 
  }

  # masses
  if(is.null(masses)== TRUE)
   {  masses = diag(dim(D)[1])*(1/(dim(D)[1]));
   }
		
  if(is.null(masses)==FALSE)
   {  masses = diag(dim(D)[1]) * masses;
   }
 
# scalar product matrices
  scalarProductMatrices = array(0,dim=c(n.rows,n.rows,n.groups))
  phi = array(0,dim=c(n.rows,n.rows,n.groups))
  for(i in 1:n.groups)
  {   phi[,,i] = diag(dim(D[,,i])[1])-(matrix(1,dim(D[,,i])[1],1) %*% diag(masses))
      scalarProductMatrices[,,i] = (-(0.5) * phi[,,i] %*% D[,,i] %*% t(phi[,,i])) 
  }

# normalization
  norm = array(0,dim=c(n.rows,n.rows,n.groups))
  for(i in 1:n.groups)
  {  if(normalization != 'None' && normalization != 'MFA' && normalization != 'SumPCA')
      {	  print(paste('WARNING: Normaliztion option not recognized. MFA was set as default'))
          norm[,,i]=scalarProductMatrices[,,i]/corePCA(scalarProductMatrices[,,i])$pdq$Dv[1]
      }
      if(normalization == 'None')
      {	  norm[,,i] = scalarProductMatrices[,,i] 
      }
      if(normalization == 'MFA')
      {	  norm[,,i]=scalarProductMatrices[,,i]/svd(scalarProductMatrices[,,i])$d[1]
		 
      }
      if(normalization == 'SumPCA')
      {	  norm[,,i]=scalarProductMatrices[,,i]/(sqrt(sum(scalarProductMatrices[,,i]*scalarProductMatrices[,,i])))
      }
  }
	
  # CMatrix
  CMatrix = diag(1,n.groups)
  for(i in 1:(n.groups))
  {  for(j in i:n.groups)
     {   rv = rvCoeff(norm[,,i],norm[,,j],type=0)
         CMatrix[i,j] = rv
         CMatrix[j,i] = rv
      }
   }

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
   D= (decomp.C$Dv)^2

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
	
# Compromise
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
   compromise.dd = (compromise.PCA$Dv)^2

# factor scores
   compromise.G = compromise.PCA$p %*% (compromise.PCA$Dd)
   rownames(compromise.G)=rownames(data)

# % of variance explained
   compromise.tau <- compromise.PCA$Dv/sum(compromise.PCA$Dv) * 100

##########################################
# Tables: Generalized PCA of X
##########################################	

# alpha weights
   table.alphaWeights <- alphaWeights

# weights and masses
    M = rep(1/(dim(data)[1]),dim(data)[1])
    
   w = c()
   for(i in 1:length(rowSums(table)))
   { w = c(w, rep(alphaWeights[i],rowSums(table)[i]))
   }
	
#general PDQ
	pdq.general = corePCA(D.mat,M=M,W=w)
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
   rownames(gpdq.loadings) = colnames(table)

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
       gpdq.partial[,,i] = D.mat[,from:to] %*% gpdq.loadings[from:to,]
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

res.distatis <- list(data = data, normalization= normalization, sorting = sorting, table = table, S=scalarProductMatrices, C = CMatrix, ci = ci, cj = cj, 
                eigs.vector = P, eigs = D, fi = G, tau = taus,  

                compromise = compromiseMatrix, compromise.ci = compromise.ci, compromise.cj = compromise.cj, compromise.eigs.vector = compromise.P, 
                compromise.eigs = compromise.dd, compromise.fi = compromise.G, compromise.tau = compromise.tau,
			
			         alphaWeights = alphaWeights, masses = M, table.partial.fi.array = gpdq.partial , table.ci = compromise.ci, 
               table.cj = table.cj, table.eigs = compromise.dd, table.tau = compromise.tau, table.eigs.vector = compromise.P, table.Q = gpdq.loadings,
			         table.fi = compromise.G, table.partial.fi = gpdq.partialFS)

return (res.distatis)
}



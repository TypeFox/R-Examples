mpANISOSTATIS.core <- function (data, num.obs, column.design, num.groups, optimization.option = 'ANISOSTATIS_Type1')
{ 
##########################################
# Inner Product
##########################################
num.groups = dim(column.design)[1]

# Scalar Product Matrices (S)
   scalarProductMatrices = array(0,dim=c(num.obs,num.obs,num.groups))
   from = 1
   for(i in 1:num.groups)
   {   from = sum(column.design[i-1,])+from
       to = from + sum(column.design[i,])-1
       scalarProductMatrices[,,i] = data[,from:to] %*% t(data[,from:to])	
   }

scalarProductMatrices2 <- array(scalarProductMatrices, dim=c(num.obs*num.obs,num.groups))

#C Matrix
  CMatrix = t(scalarProductMatrices2) %*% scalarProductMatrices2

# RV Matrix
  Norm = sqrt(apply(scalarProductMatrices2^2,2,sum))
  rvMatrix = CMatrix/(t(t(Norm)) %*% Norm)

# eigen decomposition
  C.decomp = corePCA(CMatrix)
  decomp.C = corePCA(CMatrix)$pdq

# contribution
	ci = C.decomp$ci
	if(is.null(rownames(column.design)))
  { rownames(ci) <- paste("Table",1:dim(column.design)[1],sep = "")
		table.names <- rownames(ci)
	}
  else
  { rownames(ci) <- rownames(column.design)
		table.names <- rownames(ci)
	}

#contribution
	cj = C.decomp$cj
	rownames(cj) <- rownames(ci)
	
# eigen vectors
  P = decomp.C$p
  rownames(P) <- rownames(ci)

# eigen values
  D= decomp.C$Dv

# factor scores
   G = decomp.C$p %*%  sqrt(decomp.C$Dd)
   rownames(G) <- rownames(P)

# percent of variance explained
   tau <- decomp.C$Dv/sum(decomp.C$Dv) * 100

###########################################
# Compromise
###########################################

# Masses
    Masses <- diag(num.obs)*(1/nrow(data))

# Alpha Weights
    if(optimization.option == 'ANISOSTATIS_Type1')
    {	C <- (t(data) %*% Masses %*% data) * (t(data) %*% Masses %*% data)
        pcaC <- corePCA(C)$pdq
        P <- pcaC$p
        alphaWeights <- P[,1] / sum(P[,1])
    }
	
    if(optimization.option == 'ANISOSTATIS_Type2')
    {	 C <- (t(data) %*% Masses %*% data) * (t(data) %*% Masses %*% data)
        Z <- (C %*% t(column.design)) %*% t(C %*% t(column.design))
        pcaZ <- corePCA(Z)$pdq
        P <- pcaZ$p
        alphaWeights <- P[,1] / sum(P[,1]) 
     }
	
    if(optimization.option == 'ANISOSTATIS_Power1')
    {	alphaWeights = (CMatrix %*% matrix(1,num.groups,1)) / sum(CMatrix %*% matrix(1,num.groups,1))
     }
		

# Compute Compromise
   compromiseMatrix = matrix(0,num.obs,num.obs)
   for(i in 1:num.groups)
   {   compromiseMatrix = compromiseMatrix + alphaWeights[i] * scalarProductMatrices[,,i]
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
   compromise.dd = compromise.PCA$Dv

# factor scores
   compromise.G = compromise.PCA$p %*% sqrt(compromise.PCA$Dd)
   rownames(compromise.G)=rownames(data)

# % of variance explained
   compromise.tau <- compromise.PCA$Dv/sum(compromise.PCA$Dv) * 100

##########################################
# Tables: Generalized PCA of X
##########################################	

# alpha weights
   table.alphaWeights <- alphaWeights

# weights and masses
   M = diag(1/(dim(data)[1]),dim(data)[1],dim(data)[1])

   W =matrix(0,dim(data)[2],dim(data)[2])	
	
   to_total = 0
   from_total = 1
   for(i in 1:dim(column.design)[1])
   {	 from = sum(column.design[i-1,]) + from_total
	     to = sum(column.design[i,]) + to_total
	     to_total = to
	     from_total = from
	     W[from:to,from:to] = diag(alphaWeights[i],dim(data[,from:to])[2],dim(data[,from:to])[2])
   }

#general PDQ
	pdq.general = corePCA(data,M=M,W=W)
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
   gpdq.factorscores =  general.pdq$p %*%  sqrt(general.pdq$Dd)
   rownames(gpdq.factorscores)=rownames(data)
 
# Partial Factor Scores
   gpdq.partial = array(0,dim=c(dim(data)[1],dim(gpdq.loadings)[2],num.groups))
   to_partial = 0
   from_partial = 1
   for(i in 1:dim(column.design)[1])
   {   from = sum(column.design[i-1,]) + from_partial
       to = sum(column.design[i,]) + to_partial
       to_partial = to
       from_partial = from
       gpdq.partial[,,i] = data[,from:to] %*% gpdq.loadings[from:to,]
   }
   
   gpdq.partialFS <- matrix(0,dim(data)[1]*num.groups,dim(gpdq.loadings)[2])
   to.total = 0
   for(i in 1:num.groups)
   {	from = to.total + 1
   		to = i*dim(gpdq.partial)[1]
   		to.total = to
   		gpdq.partialFS[from:to,]= gpdq.partial[,,i]
   	}

##########################################
# Results
##########################################	

res.anisostatis.core <- list(S=scalarProductMatrices, RVMatrix = rvMatrix, C = CMatrix, ci = ci, cj = cj, eigs.vector = P, eigs = D, 
                        fi = G, alphaWeights = alphaWeights, tau = tau,
          
                        compromise = compromiseMatrix, compromise.ci = compromise.ci, compromise.cj = compromise.cj, 
                        compromise.eigs.vector = compromise.P, compromise.eigs = compromise.dd, compromise.fi = compromise.G, 
                        compromise.tau = compromise.tau,
      
                        masses = M, 
                        table.partial.fi.array = gpdq.partial,table.cj = table.cj, table.ci = table.ci, 
                        table.eigs = gpdq.eigenvalues, table.tau = gpdq.inertia, table.eigs.vector = gpdq.vectors, 
                        table.loadings = gpdq.loadings,  table.fi = gpdq.factorscores,  
                        table.partial.fi = gpdq.partialFS)
    
return (res.anisostatis.core)
}

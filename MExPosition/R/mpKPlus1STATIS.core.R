mpKPlus1STATIS.core <- function(data, plus1data, num.obs, column.design, num.groups, optimization.option = 'STATIS')
{	
# Scalar Product Matrices (S)
   scalarProductMatrices = array(0,dim=c(num.obs,num.obs,num.groups))
   from = 1
   for(i in 1:num.groups)
   
   {   from = sum(column.design[i-1,])+from
       to = from + sum(column.design[i,])-1
       scalarProductMatrices[,,i] = data[,from:to] %*% t(data[,from:to])	
   } 
   
#################################
## Finding alpha *
################################   

#scalar product matrices star
	scalarProductMatrices.star = array(0,dim=c(dim(plus1data)[2],dim(plus1data)[2],num.groups))
	for(i in 1:num.groups)
	{	scalarProductMatrices.star[,,i] <- t(plus1data) %*% scalarProductMatrices[,,i] %*% plus1data
	}
	
	CubeSP2.star <- array(scalarProductMatrices.star,dim=c(num.obs*num.obs,num.groups))
    
# C Matrix star
    CMatrix.star = t(CubeSP2.star) %*% CubeSP2.star

  # RV Matrix
   Norm = sqrt(apply(CubeSP2.star^2,2,sum))
   rvMatrix.star = CMatrix.star/(t(t(Norm)) %*% Norm)
    
# eigen decomposition star
    C.decomp = corePCA(CMatrix.star)
    decomp.C.star = C.decomp$pdq

# contribution star
   ci.star = C.decomp$ci
   if(is.null(rownames(column.design)))
   {  rownames(ci.star) <- paste("Table", 1:dim(column.design)[1], sep = "")
      table.names <- rownames(ci.star)
   }
   else
   {  rownames(ci.star) <- rownames(column.design)
      table.names <- rownames(ci.star)
   }

# contribution
  cj.star = C.decomp$cj
  rownames(cj.star) <- rownames(ci.star)
    
# eigen vectors star
   P.star = decomp.C.star$p
  rownames(P.star) <- rownames(ci.star)

# eigen values star
   D.star= (decomp.C.star$Dv)

# factor scores star
   G.star = decomp.C.star$p %*%  diag(sqrt(D.star))
   rownames(G.star) <- rownames(P.star)

# percent of variance explained
   tau.star <- D.star/sum(D.star) * 100

# Alpha Weights
   if (optimization.option == "None")
   {  alphaWeights.star = matrix(1/num.groups,1,num.groups)
   }
  
   if (optimization.option == "Multitable")
   {  alphaWeights.star = matrix(1,1,num.groups)
   }
  
   if(optimization.option == 'RV_Matrix')
   {    alphaWeights.star = (corePCA(rvMatrix.star)$pdq$p)[,1]/sum((corePCA(rvMatrix.star)$pdq$p)[,1])
   }
  
   if(optimization.option == "STATIS")
   {  alphaWeights.star = P.star[,1] / sum(P.star[,1]) 
   }
  
   if(optimization.option == 'STATIS_Power1')
   {  alphaWeights.star = (CMatrix.star %*% matrix(1,num.groups,1)) / sum(CMatrix.star %*% matrix(1,num.groups,1))
   }   	

##########################################
#Compromise
##########################################

#compromise (S+)
	compromiseMatrix = apply(apply(scalarProductMatrices,c(1,2),'*',t(alphaWeights.star)),c(2,3),sum)

#analyze the compromise
	PCA.compromise <- corePCA(compromiseMatrix)
	compromise.PCA <- PCA.compromise$pdq

#contribution
	compromise.ci <- PCA.compromise$ci
	rownames(compromise.ci) = rownames(data)
	
	compromise.cj <- PCA.compromise$cj
	rownames(compromise.cj) = rownames(data)

#eigen vectors
	compromise.P = compromise.PCA$p
	rownames(compromise.P) = rownames(data)

# eigen values
   compromise.dd = compromise.PCA$Dv

# factor scores
   compromise.G = compromise.PCA$p %*% diag(sqrt(compromise.PCA$Dv))
   rownames(compromise.G)=rownames(data)
	
# % of variance explained
   compromise.tau <- compromise.PCA$Dv/sum(compromise.PCA$Dv) * 100

##########################################
# Tables: Generalized PCA of data
##########################################	

# column names	
   table.colnames <- colnames(table)

# alpha weights
   table.alphaWeights <- alphaWeights.star

# weights and masses
   M = rep(1/(dim(data)[1]),dim(data)[1])
	
   w = c()
   for(i in 1:length(rowSums(column.design)))
   { w = c(w, rep(alphaWeights.star[i],rowSums(column.design)[i]))
   }

#general PDQ
	pdq.general = corePCA(data,M=M,W=w)
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
   	rownames(gpdq.partialFS) = paste(rep(table.names,each=dim(data)[1]),rep(rownames(data)))

##########################################
# Results
##########################################  

res.KPlus1Statis.core <- list(S=scalarProductMatrices, S.star = scalarProductMatrices.star, C = CMatrix.star, ci = ci.star, cj = cj.star, eigs.vector = P.star, eigs = D.star, 
                        fi.star = G.star, alphaWeights.star = alphaWeights.star, tau.star = tau.star, rvMatrix.star = rvMatrix.star,
          
                        compromise = compromiseMatrix, compromise.ci = compromise.ci, compromise.cj = compromise.cj, 
                        compromise.eigs.vector = compromise.P, compromise.eigs = compromise.dd, compromise.fi = compromise.G, 
                        compromise.taus = compromise.tau,
      
                        masses = M, table.partial.fi.array = gpdq.partial,table.cj = table.cj, table.ci = table.ci, 
                        table.eigs = gpdq.eigenvalues, table.tau = gpdq.inertia, table.eigs.vector = gpdq.vectors, 
                        table.loadings = gpdq.loadings,  table.fi = gpdq.factorscores,  
                        table.partial.fi = gpdq.partialFS)
    
return (res.KPlus1Statis.core)
}
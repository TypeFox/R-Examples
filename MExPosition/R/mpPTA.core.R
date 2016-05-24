mpPTA.core <- function(data, num.obs, column.design, num.groups, optimization.option = 'STATIS')
{	
	# getting into array
   data.array = array(0,dim=c(num.obs,sum(column.design[1,]),num.groups))
   from = 1
   for(i in 1:num.groups) 
   {   from = sum(column.design[i-1,])+from
       to = from + sum(column.design[i,])-1
       data.array[,,i] = data[,from:to] 	
   }
   
   # Cross Product
   CubeCP <- array(0, dim=c(num.obs,num.obs,num.groups))
	for(i in 1:(num.groups))
	{ CubeCP[,,i] <- data.array[,,i] %*% t(data.array[,,i])
	}
	
   CubeCP2 <- array(CubeCP,dim=c(num.obs*num.obs,num.groups))
   # C Matrix
   CMatrix = t(CubeCP2) %*% CubeCP2
	
	#RV Matrix	 
   Norm = sqrt(apply(CubeCP2^2,2,sum))
   rvMatrix = CMatrix / (t(t(Norm)) %*% Norm)

   #checking if all entries are positive
   if(sum(CMatrix<0)>= 1)
   { stop('PTA requires all values in C matrix to be positive') 	
   }
   
	# eigen decomposition
   C.decomp = corePCA(CMatrix)
   decomp.C = C.decomp$pdq
   
   # contribution
   ci = C.decomp$ci
   #rownames(ci)=paste('assessor',1:num.groups)
   if(is.null(rownames(table)))
   {  rownames(ci) <- paste("Table", 1:dim(column.design)[1], sep = "")
	   table.names <- rownames(ci)
	}
   else
   {  rownames(ci) <- rownames(column.design)
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

   # factor scores
   G = decomp.C$p %*%  diag(sqrt(D))
   rownames(G) <- rownames(P)

   # percent of variance explained
   tau <- D/sum(D) * 100
	
   # Alpha Weights
   if (optimization.option == "None")
   {	alphaWeights = matrix(1/num.groups,1,num.groups)
   }
	
   if (optimization.option == "Multitable")
   {	alphaWeights = matrix(1,1,num.groups)
   }
	
   if(optimization.option == 'RV_Matrix')
   {    alphaWeights = (corePCA(rvMatrix)$pdq$p)[,1]/sum((corePCA(rvMatrix)$pdq$p)[,1])
   }
	
   if(optimization.option == "STATIS")
   {	alphaWeights = P[,1] / sum(P[,1])	
   }
	
   if(optimization.option == 'STATIS_Power1')
   {	alphaWeights = (CMatrix %*% matrix(1,num.groups,1)) / sum(CMatrix %*% matrix(1,num.groups,1))
   }
	#make sure all positive

##########################################
#Compromise
##########################################

   #compromise (S+)
	compromiseMatrix <- apply(apply(CubeCP,c(1,2),'*',t(alphaWeights)),c(2,3),sum)
	print(compromiseMatrix)

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
   rownames(compromise.G)=rownames(column.design)
   	
   # how each assessors interprets the space
	f <- compromiseMatrix %*% (compromise.PCA$p %*% diag(compromise.PCA$Dd^(-0.5)))

   # % of variance explained
   compromise.tau <- compromise.PCA$Dv/sum(compromise.PCA$Dv) * 100

##########################################
# Tables: Generalized PCA of X
##########################################	

   # column names	
   table.colnames <- colnames(column.design)

   # alpha weights
   table.alphaWeights <- alphaWeights

   # weights and masses
   M = diag(1/(dim(data)[1]),dim(data)[1],dim(data)[1])
	
   w = c()
   for(i in 1:length(rowSums(column.design)))
   { w = c(w, rep(alphaWeights[i],rowSums(column.design)[i]))
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
   rownames(gpdq.factorscores)= rownames(data)

   
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

res.pta.core <- list(S=CubeCP, RVMatrix = rvMatrix, C = CMatrix, ci = ci, cj = cj, eigs.vector = P, eigs = D, 
                fi = G, alphaWeights = alphaWeights, tau = tau,
          
               compromise = compromiseMatrix, compromise.ci = compromise.ci, compromise.cj = compromise.cj, 
               compromise.eigs.vector = compromise.P, compromise.eigs = compromise.dd, compromise.fi = compromise.G, 
               compromise.tau = compromise.tau,
      
               masses = M, table.partial.fi.array = gpdq.partial,table.cj = table.cj, table.ci = table.ci, 
               table.eigs = gpdq.eigenvalues, table.tau = gpdq.inertia, table.eigs.vector = gpdq.vectors, 
               table.loadings = gpdq.loadings,  table.fi = gpdq.factorscores,  
               table.partial.fi = gpdq.partialFS)
    
return (res.pta.core)
}






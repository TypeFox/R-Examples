mpDOACT.STATIS.core <- function(dataset1, column.design.1, dataset2, column.design.2)
{ 
  num.groups <- dim(column.design.1)[1]

   # cross product of dataset 1
   CubeSP.1= array(0,dim=c(dim(dataset1)[1],dim(dataset1)[1],dim(column.design.1)[1]))
   from = 1
   for(i in 1:dim(column.design.1)[1]) 
   {   from = sum(column.design.1[i-1,])+from
       to = from + sum(column.design.1[i,])-1
       CubeSP.1[,,i] = dataset1[,from:to] %*% t(dataset1[,from:to])
   }	
   
   # cross product of dataset 2
   CubeSP.2= array(0,dim=c(dim(dataset2)[1],dim(dataset2)[1],dim(column.design.2)[1]))
   from = 1
   for(i in 1:dim(column.design.2)[1]) 
   {   from = sum(column.design.2[i-1,])+from
       to = from + sum(column.design.2[i,])-1
       CubeSP.2[,,i] = dataset2[,from:to] %*% t(dataset2[,from:to])
   }
   
  CubeSP2.1 <- array(CubeSP.1,dim=c(dim(dataset1)[1]*dim(dataset1)[1],dim(column.design.1)[1]))
  CubeSP2.2 <- array(CubeSP.2,dim=c(dim(dataset2)[1]*dim(dataset2)[1],dim(column.design.2)[1]))
 

### NEED TO ADD RV HERE


  # C Matrix using cross products of dataset 1 and dataset 2
  CMatrix = t(CubeSP2.1) %*% CubeSP2.2

  C.decomp = corePCA(CMatrix)
  decomp.C = C.decomp$pdq

  # contribution
   ci = C.decomp$ci
   #rownames(ci)=paste('assessor',1:num.groups)
   if(is.null(rownames(column.design.1)))
   {  rownames(ci) <- paste("Table", 1:dim(column.design.1)[1], sep = "")
      table.names <- rownames(ci)
  }
    else
   {  rownames(ci) <- rownames(column.design.1)
       table.names <- rownames(ci)
    }

# contribution
  cj = C.decomp$cj
  rownames(cj) <- rownames(ci)
    
# eigen vectors
   P = decomp.C$p
   rownames(P) <- rownames(ci)

   Q = decomp.C$q
   rownames(Q) <- rownames(ci)

# eigen values
   D= (decomp.C$Dv)

# factor scores
   G = decomp.C$p %*%  diag(sqrt(D))
   rownames(G) <- rownames(P)

# percent of variance explained
   tau <- D/sum(D) * 100

  # alpha and beta weights 
  alphaWeights = P[,1] / sum(P[,1])
  betaWeights = Q[,1] / sum(Q[,1])


##########################################
#Compromise
##########################################

#compromise (S+)
	compromiseMatrix.1 <- apply(apply(CubeSP.1,c(1,2),'*',t(alphaWeights)),c(2,3),sum)
	compromiseMatrix.2 <- apply(apply(CubeSP.2,c(1,2),'*',t(betaWeights)),c(2,3),sum)

#analyze the compromise
	PCA.compromise.1 <- corePCA(compromiseMatrix.1)
	compromise.PCA.1 <- PCA.compromise.1$pdq
	
	PCA.compromise.2 <- corePCA(compromiseMatrix.2)
	compromise.PCA.2 <- PCA.compromise.2$pdq

#contribution
	compromise.ci.1 <- PCA.compromise.1$ci
	rownames(compromise.ci.1) = rownames(dataset1)
	
	compromise.cj.1 <- PCA.compromise.1$cj
	rownames(compromise.cj.1) = rownames(dataset1)
	
	compromise.ci.2 <- PCA.compromise.2$ci
	rownames(compromise.ci.1) = rownames(dataset1)
	
	compromise.cj.2 <- PCA.compromise.2$cj
	rownames(compromise.cj.2) = rownames(dataset1)
	
#eigen vectors
	compromise.P.1 = compromise.PCA.1$p
	rownames(compromise.P.1) = rownames(dataset1)
	
	compromise.P.2 = compromise.PCA.2$p
	rownames(compromise.P.2) = rownames(dataset1)
	
# eigen values
   compromise.dd.1 = (compromise.PCA.1$Dv)
   compromise.dd.2 = (compromise.PCA.2$Dv)

# factor scores
   compromise.G.1 = compromise.PCA.1$p %*% diag(sqrt(compromise.PCA.1$Dv))
   rownames(compromise.G.1)=rownames(dataset1)
   
   compromise.G.2 = compromise.PCA.2$p %*% diag(sqrt(compromise.PCA.2$Dv))
   rownames(compromise.G.2)=rownames(dataset1)
	
# % of variance explained
   compromise.tau.1 <- compromise.PCA.1$Dv/sum(compromise.PCA.1$Dv) * 100
   compromise.tau.2 <- compromise.PCA.2$Dv/sum(compromise.PCA.2$Dv) * 100

##########################################
# Tables: Generalized PCA of data
##########################################	

# alpha weights
   table.1.alphaWeights <- alphaWeights
   table.2.betaWeights <- alphaWeights

# weights and masses
   M.1 = rep(1/(dim(dataset1)[1]),dim(dataset1)[1])
   M.2 = rep(1/(dim(dataset2)[1]),dim(dataset2)[1]) 

   w.1 = c()
   for(i in 1:length(rowSums(column.design.1)))
   { w.1 = c(w.1, rep(alphaWeights[i],rowSums(column.design.1)[i]))
   }

   w.2 = c()
   for(i in 1:length(rowSums(column.design.2)))
   { w.2 = c(w.2, rep(betaWeights[i],rowSums(column.design.2)[i]))
   }

#general PDQ
	pdq.general.1 = corePCA(dataset1,M=M.1,W=w.1)
	general.pdq.1 = pdq.general.1$pdq
	
	pdq.general.2 = corePCA(dataset2,M=M.2,W=w.2)
	general.pdq.2 = pdq.general.2$pdq


# contribution
	table.1.ci = pdq.general.1$ci
	table.2.ci = pdq.general.2$ci

# contribution
	table.1.cj = pdq.general.1$cj
	table.2.cj = pdq.general.2$cj
   
# Eigen vectors of the tables
   gpdq.vectors.1 = general.pdq.1$p
   gpdq.vectors.2 = general.pdq.2$p

# Eigen values of the tables 
   gpdq.eigenvalues.1 = (general.pdq.1$Dd)^2
   gpdq.eigenvalues.2 = (general.pdq.2$Dd)^2
	
# Inertia
   gpdq.inertia.1 = ((general.pdq.1$Dv) / sum(general.pdq.1$Dv))*100
   gpdq.inertia.2 = ((general.pdq.2$Dv) / sum(general.pdq.2$Dv))*100

# Loadings of the tables
   gpdq.loadings.1 = general.pdq.1$q
   rownames(gpdq.loadings.1) = colnames(dataset1)
   
   gpdq.loadings.2 = general.pdq.2$q
   rownames(gpdq.loadings.2) = colnames(dataset2)

# Factor scores of the tables
   gpdq.factorscores.1 =  general.pdq.1$p %*%  (general.pdq.1$Dd)
   rownames(gpdq.factorscores.1)=rownames(dataset1)
   
   gpdq.factorscores.2 =  general.pdq.2$p %*%  (general.pdq.2$Dd)
   rownames(gpdq.factorscores.2)=rownames(dataset2)
   
# Partial Factor Scores
   gpdq.partial.1 = array(0,dim=c(dim(dataset1)[1],dim(gpdq.loadings.1)[2],dim(column.design.1)[1]))
   to_partial = 0
   from_partial = 1
   for(i in 1:dim(column.design.1)[1])
   {   from = sum(column.design.1[i-1,]) + from_partial
       to = sum(column.design.1[i,]) + to_partial
       to_partial = to
       from_partial = from
       gpdq.partial.1[,,i] = dataset1[,from:to] %*% gpdq.loadings.1[from:to,]
   }
   
   gpdq.partial.2 = array(0,dim=c(dim(dataset2)[1],dim(gpdq.loadings.2)[2],dim(column.design.2)[1]))
   to_partial = 0
   from_partial = 1
   for(i in 1:dim(column.design.2)[1])
   {   from = sum(column.design.2[i-1,]) + from_partial
       to = sum(column.design.2[i,]) + to_partial
       to_partial = to
       from_partial = from
       gpdq.partial.2[,,i] = dataset2[,from:to] %*% gpdq.loadings.2[from:to,]
   }
   
   gpdq.partialFS.1 <- matrix(0,dim(dataset1)[1]*dim(column.design.1)[1],dim(gpdq.loadings.1)[2])
   to.total = 0
   for(i in 1:dim(column.design.1)[1])
   {	from = to.total + 1
   		to = i*dim(gpdq.partial.1)[1]
   		to.total = to
   		gpdq.partialFS.1[from:to,]= gpdq.partial.1[,,i]
   	}
   	# rownames(gpdq.partialFS.1) = paste(rep(table.names,each=dim(data)[1]),rep(rownames(data)))

   gpdq.partialFS.2 <- matrix(0,dim(dataset2)[1]*dim(column.design.2)[1],dim(gpdq.loadings.2)[2])
   to.total = 0
   for(i in 1:dim(column.design.2)[1])
   {	from = to.total + 1
   		to = i*dim(gpdq.partial.2)[1]
   		to.total = to
   		gpdq.partialFS.2[from:to,]= gpdq.partial.2[,,i]
   	}


##########################################
# Results
##########################################	

res.doact.statis.core <- list(S.1 = CubeSP.1, S.2 = CubeSP.2, C = CMatrix, ci = ci, cj = cj, eigs.vector = P, eigs = D, 
                        fi = G, alphaWeights = alphaWeights, tau = tau,

                        alphaWeights = alphaWeights, betaWeights = betaWeights,
						
						            compromiseMatrix.1 = compromiseMatrix.1, compromise.ci.1 = compromise.ci.1, compromise.cj.1 = compromise.cj.1, compromise.P.1 = compromise.P.1,
                        compromise.eigs.value.1 = compromise.dd.1, compromise.fi.1 = compromise.G.1, compromise.tau.1 = compromise.tau.1, 
                        compromiseMatrix.2 = compromiseMatrix.2, compromise.ci.2 = compromise.ci.2, compromise.cj.2 = compromise.cj.2, compromise.P.2 = compromise.P.2,  
						            compromise.eigs.value.2 = compromise.dd.2, compromise.fi.2 = compromise.G.2, compromise.tau.2 = compromise.tau.2,
						
                        masses.1 = M.1, 
                        table.partial.fi.array.1 = gpdq.partial.1, table.cj.1 = table.1.cj, table.ci.1 = table.1.ci, 
                        table.eigs.1 = gpdq.eigenvalues.1, table.tau.1 = gpdq.inertia.1, table.eigs.vector.1 = gpdq.vectors.1, 
                        table.loadings.1 = gpdq.loadings.1,  table.fi = gpdq.factorscores.1,  
                        table.partial.fi.1 = gpdq.partialFS.1,

                        masses.2 = M.2, 
                        table.partial.fi.array.2 = gpdq.partial.2, table.cj.2 = table.2.cj, table.ci.2 = table.2.ci, 
                        table.eigs.2 = gpdq.eigenvalues.2, table.tau.2 = gpdq.inertia.2, table.eigs.vector.2 = gpdq.vectors.2, 
                        table.loadings.2 = gpdq.loadings.2,  table.fi.2 = gpdq.factorscores.2,  
                        table.partial.fi.2 = gpdq.partialFS.2)

return (res.doact.statis.core)

}
reduce <-
function(beta, indices, assured.intercept, accuracy)
{
C <- diag(1, length(beta))
dimnames(C) <- list(rep(1:ncol(indices),times=indices["index1",]), rownames(beta))
A <- C

indexs <- as.matrix(indices[-c(1,3),])
index1 <- indices["index1",]
index2 <- indices["index2",]
index2b <- indices["index2b",]
indexCo <- rep(1:ncol(indices), indices["index1",]) # zu welchem coefficient welche Kovariable?!

i <- 1
while (i <= ncol(C) ) { 
  if ((colSums(indexs))[eval(parse(text=rownames(C)[which(C[,i]==1)]))]!=0 )  # exclude only penalized coefficients        
  { 
      j <- indexCo[which(C[,i]==1)] # bei welcher variable ist i grade?!
      married <- c( (sum(index1[1:j]) - index1[j] + 1) : sum(index1[1:j]) ) # welche koefs gehoeren dazu?

      # exclude zero coefficients
      if (sum(beta*C[,i])==0 && ifelse(index2b[j]==0 && index2[j]!=0, sum(as.matrix(C[,-i])[married,])!=0, TRUE))  # keep assured.intercept, even if zero
      { 
          namen <- colnames(C)[-i]
          C <- as.matrix(C[,-i])
          A <- as.matrix(A[,-i])
          colnames(C) <- colnames(A) <- namen
          
      # exclude equal coefficients     
      } else {
          equal <- which((beta - rep(sum(beta*C[,i]), length(beta)))==0)
          # married <- which(rownames(C)==eval(parse(text=rownames(C)[which(C[,i]==1)])) )
          null <- equal[which(equal %in% married)]
          # A <- C
          A[null,i] <- 1
          raus <- which(colSums(A * A[,i])==1)
          if (length(null)>1 && ifelse(assured.intercept, 
               sum(as.matrix(C[,-raus])[c(1:indices["index1",1]),])!=0, TRUE) && 
               any(indices[c("index2", "index3", "index4", "index5", "index7", "index8"),j]!=0)
             )
             {
             namen <- colnames(C)[-raus]
             C <- as.matrix(C[,-raus])
             A <- as.matrix(A[,-raus])
             colnames(C) <- colnames(A) <- namen
             }
          if (all(indices[c("index2", "index3", "index4", "index5", "index7", "index8"),j]==0)) A[,i] <- C[,i]
          i <- i+1
      }
      
  } else {i <- i+1}
}

return(list(C=A, A=C, beta=t(C) %*% beta))

}

# in the function, it holds:
# C to reduce beta
# A to reduce x

# in the output, it holds:
# A to reduce beta
# C to reduce x

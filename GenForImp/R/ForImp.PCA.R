ForImp.PCA <- function(mat, stand=FALSE, cor=FALSE, r=2, probs = seq(0, 1, 0.1), q="10%", tol=1e-6){
     # mat = a matrix or a dataframe
     # stand = should mat be standardized? Default is FALSE
     # r = order of the Minkowski distance (r >= 1). In particular, r=Inf is Chebyshev distance
     # cor = should principal components be extracted from correlation matrix? Default is FALSE
     # probs = sequence of probabilities for computing donor quantiles. The default is: probs=(0, .1, .2, ..., 1)
     # q = a choosen donor quantile. The default is: q="10%", i.e. the first 10% of donors.
     # tol = tolerance factor controlling for floating points. Default is: tol=1e-6
origMatDat <- MatDat <- as.matrix(mat)
# number of units
n <- nrow(origMatDat)
# number of variables
p <- ncol(origMatDat)
#
# row and column names
if(is.null(rownames(MatDat)))
      rownames(MatDat) <- 1:n
if(is.null(colnames(MatDat))) 
      colnames(MatDat) <- paste("x", 1:p, sep=".")
case.name <- rownames(MatDat)
var.name <- colnames(MatDat)
#
# Standardization
if(stand){
        # mean vector and variance-covariance matrix
        mean.vec <- apply(origMatDat, 2, function(x) mean(x, na.rm=T) )
        sd.vec <- apply(origMatDat, 2, function(x) sd(x, na.rm=TRUE)*sqrt( (n-1)/n ) )
        #
MatDat <- scale(origMatDat, center=mean.vec, scale=sd.vec )
        names(mean.vec) <- names(sd.vec) <- var.name
}
#
#
# Remark: the maximum number of within-rows missing values might not coincide with the number of submatrices #
# Map of missing values
ind.NA <- is.na(MatDat)
count.miss <- as.numeric(names( table(apply(ind.NA, 1, function(x) table(x)["TRUE"])) ))
max.na <- max(count.miss)
# number of submatrices with missing values (the complete matrix is excluded) #
nr.Xk <- length(count.miss)
#
list.ind.na <- vector("list", nr.Xk + 1)
# Indices of units in the complete submatrix (complete units)
list.ind.na[[1]] <- which(apply(ind.NA, 1, function(x) table(x)["FALSE"]== p ))
# Indices of units in the incomplete submatrices (incomplete units)
for(i in 2:(nr.Xk + 1)){
      list.ind.na[[i]] <- which(apply(ind.NA, 1, function(x) table(x)["TRUE"]== count.miss[i-1] ))
}
#
# PARTITION OF MAT IN SUBMATRICES
list.sub.mtr.k <- lapply(1:(nr.Xk + 1), function(i){ 
                           mat.i <- MatDat[list.ind.na[[i]],]
                           if(is.null(dim(mat.i)))
                              dim(mat.i) <- c(1, p)
                              dimnames(mat.i) <- list(list.ind.na[[i]], var.name)
                            mat.i
                            })
#
# Complete matrix
X0 <- list.sub.mtr.k[[1]]
nX0 <- nrow(X0)
#
list.donor.k <- vector("list", 3)
#
#--- 
#--- ROUTINE STARTS ---#
#--- 
for(k in 1:nr.Xk){
   Xk <- list.sub.mtr.k[[k+1]]
   # number of units in Xk
   nXk <- nrow(Xk)
#
# PCA on the complete matrix   
#
if(nX0 < p){ # if the rank of the cor (or cov) matrix is not full #
   if(cor) corX0 <- cor(X0) else corX0 <- cov(X0)*((nX0 - 1)/nX0)
   corX0.eig <- eigen(corX0)
   pca.eigen <- corX0.eig$values
   # controlling for numerical negative eigenvalues
   pca.eigen[pca.eigen < 0] <- 0
   pca.coeff <- corX0.eig$vectors
   names(pca.eigen) <- paste("Comp", 1:p, sep=".")
   colnames(pca.coeff) <- paste("Comp", 1:p, sep=".")
   rownames(pca.coeff) <- var.name 
   }
#
# By default function 'princomp' uses divisor N for the covariance matrix (if cor=F). 
else{ # the rank of the cor (or cov) matrix is full #
   pca.compl <- princomp(X0, cor=cor, scores=T)
   pca.coeff <- unclass(pca.compl$loadings)
   pca.eigen <- (pca.compl$sdev)^2
   }
#
## weights are given by the square root of the normalized variance of each component #
   weight.eig <- sqrt( pca.eigen /sum(pca.eigen) )
#

#--- Set-up of k-combinations of the p variables #
   ind.var.k <- combn(var.name, count.miss[k] )
   ncol.ind.var.k <- ncol(ind.var.k)
#  
list.ind.var.k <- lapply(1:ncol.ind.var.k, function(j){
            Xk.j <- Xk[,ind.var.k[,j]]
            ind.unit <- is.na(Xk.j)
            if(nXk == 1){ 
                ind <- all(ind.unit)
                names(ind) <- rownames(Xk)
                 }
           else{ 
             if(is.null(dim(ind.unit))){ind <- ind.unit}
                else{ ind <- apply(ind.unit,1,all) }
               }
           ind 
          })
#
#--- Retaining only k-combinations of the p variables with missing values (so-called "iota set")#
#--- Submatrices with more than one row
#
if(nXk != 1){
   list.ret.ind.var.k <- lapply(1:ncol.ind.var.k, function(j){
      t.false.j <- table(list.ind.var.k[[j]])["FALSE"]
   if(!is.na(t.false.j)){
             if( t.false.j == nXk ) NULL
             else{ list.ind.var.k[[j]] }
          }
     else list.ind.var.k[[j]]
             })
#
   dummy.ret.ind.var.k <- sapply(1:ncol.ind.var.k, function(j){
             ifelse(is.null(list.ret.ind.var.k[[j]]), 0, 1)         
     })
   attr(dummy.ret.ind.var.k, "names") <- NULL
}
#
#--- Submatrix with one row 
else{
     list.ret.ind.var.k <- list.ind.var.k
     dummy.ret.ind.var.k <- sapply(1:ncol.ind.var.k, function(j){
                              ifelse(list.ind.var.k[[j]] != "FALSE", 1, 0)         
                             })
     attr(dummy.ret.ind.var.k, "names") <- NULL
     }
#
   mat.ret.ind.var.k <- matrix(unlist(list.ret.ind.var.k), ncol=length(list.ind.na[[k+1]]), byrow=T)
   colnames(mat.ret.ind.var.k) <- list.ind.na[[k+1]]

#
ret.ind.var.k <- matrix(ind.var.k[,dummy.ret.ind.var.k==1], ncol=sum(dummy.ret.ind.var.k))
ncol.ret.ind.var.k <- ncol(ret.ind.var.k)

#--- Pseudo PC scores on the complete submatrix #
list.mtr.PPC.compl <- lapply(1:ncol.ret.ind.var.k, function(j){ 
                           ind <-  setdiff(var.name, ret.ind.var.k[,j]) 
                           as.matrix(X0[,ind]) %*% pca.coeff[ind,] 
                           })
#
#--- Pseudo PC scores on the incomplete submatrix #
#---
#--- Submatrices with one row 
#---
#
if(nXk == 1){
        list.mtr.reord.k  <- Xk
        list.mtr.reord.k[,ret.ind.var.k] <- 0
        list.mtr.PPC.incompl.k <- list.mtr.reord.k %*% pca.coeff
#
#--- Comparing PPC scores btw. complete and incomplete matrices #
#
# Construction of the distance matrix 
  repl.unit.list.k <- matrix(rep( as.matrix(list.mtr.PPC.incompl.k), nrow(list.mtr.PPC.compl[[1]]) ), 
       nrow(list.mtr.PPC.compl[[1]]), p, byrow=T)
## Weighted Minkowski distances of order r ##
  if(r == Inf){
  # Chebyshev distance #
  dist.mtr.list.k <- apply(abs(repl.unit.list.k - list.mtr.PPC.compl[[1]]), 1, max)
  }
  else{
  dist.mtr.list.k <- apply(abs(repl.unit.list.k - list.mtr.PPC.compl[[1]])^r * weight.eig^r, 1, sum)^(1/r) 
  }
#
#--- Detection of donors #
#
# Quantile table #
     list.donor.k[[1]] <- quantile(dist.mtr.list.k, probs = probs )
# Selection of donors #
     list.donor.k[[2]] <- dist.mtr.list.k[(dist.mtr.list.k <= list.donor.k[[1]][q]*(1+tol) )]     
# Donors' names #
     list.donor.k[[3]] <- names(list.donor.k[[2]])
#
#--- Imputation #
#
   sel.mat <- X0[list.donor.k[[3]], ret.ind.var.k]
#
#  In case of a single donor #
#
   if(length(list.donor.k[[3]]) == 1) imput.j <- sel.mat
#  more than one donor #
   else{
   #--- in case of perfect coincidence btw. PPC scores of complete and incomplete units ---#
        if(any((list.donor.k[[2]])==0)){
          list.donor.k[[2]][which((list.donor.k[[2]])==0)] <- 0.0000001
          }
#
     #-- one variable --#
   if(is.null(dim(sel.mat))) imput.j <- sum(sel.mat * 1/list.donor.k[[2]]) / sum(1/list.donor.k[[2]])
     else{
     #-- more than one variable --#
   imput.j <- apply(sel.mat * 1/list.donor.k[[2]],2,sum) / sum(1/list.donor.k[[2]]) }
   }#
   dim(imput.j) <- c(1, length( ret.ind.var.k))
   dimnames(imput.j) <- list(list.ind.na[[k+1]], ret.ind.var.k)
        if(stand){
          imput.j <- imput.j * sd.vec[ret.ind.var.k ] + mean.vec[ret.ind.var.k]
    }
#
# Arranging matrices #
   mtr.imput <- Xk
   mtr.imput[ , ret.ind.var.k]  <- imput.j
   mtr.imput
#
} #--- end of submatrices with one row 
#---
#--- Submatrices with more than one row
#---
else{
   list.mtr.reord.k <- lapply(1:ncol.ret.ind.var.k, function(j){
#
           ind <- mat.ret.ind.var.k[j,]
     matr.j <- Xk[ind,]
           if(is.null(dim(matr.j))){ 
    dim(matr.j) <- c(1, p) 
          dimnames(matr.j) <- list(rownames(Xk)[ind], var.name)
             }
     matr.j[ , ret.ind.var.k[,j]] <- 0
     matr.j
          })
   list.mtr.PPC.incompl.k <- lapply(1:length(list.mtr.reord.k), function(j) list.mtr.reord.k[[j]] %*% pca.coeff )
#
#--- Comparing pseudo pc scores btw. complete and incomplete matrices #
# Construction of the distance matrix 
   repl.unit.list.k <- lapply(1:length(list.mtr.PPC.incompl.k), function(j){
                         lapply(1: nrow(list.mtr.PPC.incompl.k[[j]]), 
                             function(i) matrix(rep( as.matrix(list.mtr.PPC.incompl.k[[j]][i,]), nrow(list.mtr.PPC.compl[[j]]) ), 
                             nrow(list.mtr.PPC.compl[[j]]), p, byrow=T))
}) 
  #
   dist.mtr.list.k <- lapply(1:length(list.mtr.PPC.incompl.k), function(j){
#
         dist.i <- t( sapply(1:length(repl.unit.list.k[[j]]), function(i) 
               #
               if(r == Inf){ 
               ## Chebyshev distance (r=Inf) ##
            (apply(abs(repl.unit.list.k[[j]][[i]] - list.mtr.PPC.compl[[j]]), 1, max)) }
                #
               else{#
                ## Weighted Minkowski distances of order r ##
            (apply(abs(repl.unit.list.k[[j]][[i]] - list.mtr.PPC.compl[[j]])^r * weight.eig^r, 1, sum))^(1/r) } 
                #
                )) #
          rownames(dist.i) <- rownames(list.mtr.PPC.incompl.k[[j]]) 
      dist.i
})
#
#--- Detection of donors #
#
# Quantile table #
  list.donor.k[[1]] <- lapply(1:length(dist.mtr.list.k), function(j){
            tab.j <- t(sapply(1:nrow(dist.mtr.list.k[[j]]), 
                     function(i) quantile(dist.mtr.list.k[[j]][i,], probs = probs ))) 
            rownames(tab.j) <- rownames(dist.mtr.list.k[[j]])
            tab.j
            })
# Selection of donors #
  list.donor.k[[2]] <- lapply(1:length(dist.mtr.list.k), function(j){
      lapply(1:nrow(dist.mtr.list.k[[j]]), function(i) 
            dist.mtr.list.k[[j]][i, ( dist.mtr.list.k[[j]][i,] <= list.donor.k[[1]][[j]][i, q]*(1+tol) )] )
            })
# Donors' names #
  list.donor.k[[3]] <- lapply(1:length(dist.mtr.list.k), function(j){         
            nam.don.j <- lapply(1:nrow(dist.mtr.list.k[[j]]), function(i){ 
            #---- one donor  ----#
                   if(length(list.donor.k[[2]][[j]][[i]])==1){
                         nam.dons <- names(dist.mtr.list.k[[j]][i, ])
                         nam.dons[ dist.mtr.list.k[[j]][i,] <= list.donor.k[[1]][[j]][i, q]*(1+tol) ]                        
                    }
            #--- more than one donor ---#
                  else{              
                      names(list.donor.k[[2]][[j]][[i]])  
                      }
                  })
            nam.don.j
            })
#
#--- Imputation #
  list.imput.k <- lapply(1:length(dist.mtr.list.k), function(j){
            ret.ind.var.k.j <- ret.ind.var.k[,j]
            row.nam.j   <- rownames(dist.mtr.list.k[[j]])
            imput.j <- t(sapply(1:nrow(dist.mtr.list.k[[j]]), function(i) {
                       mat.ij <- X0[list.donor.k[[3]][[j]][[i]], ret.ind.var.k.j]
                       if(is.null(dim(mat.ij))){ 
                           dim(mat.ij) <- c(length(list.donor.k[[3]][[j]][[i]]), length(ret.ind.var.k.j))
                           dimnames(mat.ij) <- list(list.donor.k[[3]][[j]][[i]], ret.ind.var.k.j)
                        } #
                     #--- in case of perfect coincidence btw. PPC scores of complete and incomplete units ---#
                       if(any((list.donor.k[[2]][[j]][[i]])==0)){
                           list.donor.k[[2]][[j]][[i]][which((list.donor.k[[2]][[j]][[i]])==0)] <- 0.0000001
                       }
                       apply(mat.ij * 1/list.donor.k[[2]][[j]][[i]],2,sum) / sum(1/list.donor.k[[2]][[j]][[i]]) 
                       } ) ) #-- end sapply
             if( dim(imput.j)[1]==1 & length(unique(colnames(imput.j))) == 1 ){
                  imput.j <- t(imput.j)
              }
             dimnames(imput.j) <- list(row.nam.j, ret.ind.var.k.j) 
             imput.j
      })
#
# Arranging matrices #
mtr.imput <- Xk
#
for(j in 1:length(list.imput.k)){ 
        val.imput.k.j <- list.imput.k[[j]]
        ret.ind.var.k.j <- ret.ind.var.k[,j]
        mtr.imput[rownames(dist.mtr.list.k[[j]]), ret.ind.var.k.j]  <- val.imput.k.j
 }
#
} #--- end of submatrices with more than one row 
#
#----- Updating the complete part of data with imputed values #
X0 <- rbind(X0, mtr.imput)
X0 <- X0[order(as.numeric(rownames(X0))), ]
nX0 <- nrow(X0)
  } 
#--- 
##--- ROUTINE ENDS ---#
#--- 
#
if(stand){
         X0 <- X0 %*% diag(sd.vec) + matrix(rep(mean.vec,n), nrow=n, byrow=T)
}
#
dimnames(X0) <- list(case.name, var.name)
X0
}

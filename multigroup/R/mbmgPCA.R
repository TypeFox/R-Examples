#' @title  multiblock and multigroup Principal Component Analysis
#' 
#' @description 
#'  multiblock and multigroup PCA (mbmgPCA)
#' 
#' @param Data a numeric (quantitative) matrix or data frame
#' @param Group a vector of factors associated with group structure
#' @param nBlock a vector of number of variables in each block
#' @param Block.name vector of name of blocks 
#' @param ncomp number of components, if NULL number of 
#' components is equal to min(rank(Data), M-1)
#' @param niter number of iteration, if NULL number of 
#' iteration is equal to 10
#' @param ScaleGroup scaling variables in each group and block, by defalt is FALSE
#' @param ScaleDataA scaling variables in each block after group preprocessing, by defalt is FALSE
#' @param ScaleDataB scaling variables in each block befor group preprocessing, by defalt is FALSE
#' @param norm normalize each block, by defalt is FALSE
#' @return list with the following results:
#' @return \item{K.Data}{Block data}
#' @return \item{concat.Data}{Concatenated data}
#' @return \item{concat.block.Data}{Block concatenated data}
#' @return \item{res.iter}{Result of iteration} 
#' @return \item{CRIT.h}{Maximization criterion for each diemnsion}
#' @return \item{CRIT}{Maximization criterion}
#' @return \item{crit.group}{Maximization criterion associated with each group}
#' @return \item{crit.block}{Maximization criterion associated with each block}
#' @return \item{omega}{Weight of each block in construction of common scores}
#' @return \item{block.common.loading}{Common loadings for each block}
#' @return \item{block.group.loadings}{Partial loadings for each block and group}
#' @return \item{similarity}{Similarity among common and partial loadings for each block}
#' @return \item{global.scores}{Global scores among blocks}
#' @return \item{block.scores}{Scores for each block}
#' @return \item{block.group.scores}{Scores for each block and group}
#' @return \item{block.scores}{Scores for each block}
#' @return \item{global.expvar}{Global explained variance}
#' @return \item{cum.exp.var.block.group}{Cumulative explained variance for each block and group}
#' @seealso \code{\link{mgPCA}} 
#' @export
#' @references A. Eslami, E. M. Qannari, A. Kohler and S. Bougeard, Under Review. Multivariate data analysis of multi-groups datasets. Application to sensory analysis,
#'  \emph{Chemolab}, 25, 108-123.
#' 
#'     
#' @examples
#' data(wine)
#' Select=c(which(wine[,2]=="Env1"),which(wine[,2]=="Env2"),which(wine[,2]=="Reference"))
#' WineData = wine[Select,-c(1,2)]
#' Group <- as.factor(c(rep("Env1",7), rep("Env2",5), rep("Reference",7)))
#' nBlock <- c(5, 3, 10, 9)
#' BlockNames    <- c("Olfaction at rest", "Vision", "Olfaction  after shaking", "Taste")
#' res = mbmgPCA(Data = WineData, Group, nBlock , ncomp=5)
mbmgPCA <- function(Data, Group, nBlock, Block.name=NULL, ncomp=NULL, niter=NULL, ScaleGroup=FALSE, ScaleDataA=FALSE, 
                    ScaleDataB=FALSE, norm=FALSE){
  
  #============================================================================
  #                              Checking inputs
  #============================================================================
  Group = as.factor(Group)
  
  
  if(!is.data.frame(Data) & !is.matrix(Data))
    stop("\nOops the class of Data must be data.frame or matrix")
  
  # if(has_missing(Data))
  if(sum(is.na(Data)) != 0)
    stop("Oops there are missing values")
  
  
  if(nrow(Data) != length(Group))
    stop("\nOops the length of group is different of number of rows in data")
  
  
  if(length(levels(Group)) == 1)
    stop("\nGroups must have more than one level") 
  
  
  if (class(Data) == 'data.frame') {
    Data = as.matrix(Data)
  }
  
  
  if(is.null(niter)) {niter=10}  
  
  if(is.null(colnames(Data))) {
    colnames(Data) = paste('V', 1:ncol(Data), sep='')
  }
  
  
  
  rownames(Data) = Group                 #----rownames of data=groups
  M = length(levels(Group))              #----number of groups: M
  n = as.vector(table(Group))            #----number of individuals in each group
  N = sum(n)                             #----number of individuals
  K = length(nBlock)                     #----number of blocks
  if(is.null(Block.name)) {Block.name=paste("Block", 1:K, sep="")}
  P = nBlock
  PP = sum(P)
  if(is.null(ncomp)) {min(qr(Data)$rank, M-1)}  # or H=qr(g_data_group)$rank
  
  #---------------------------------------------------------------------------
  #                       2. prepare data
  #---------------------------------------------------------------------------
  
  res.prepare = prepare_data(Data, Group, nBlock, Block.name, ScaleGroup, ScaleDataA, ScaleDataB, norm)
  K.Data = res.prepare$K.Data                     #----------- split each block to M groups of individuals
  concat.Data = res.prepare$concat.Data           #----------- concatenate dataset by row as groups for each block
  concat.block.Data = res.prepare$concat.block.Data
  
  
  #============================================================================
  #        	                       Outputs
  #============================================================================
  res <- list(
    K.Data = K.Data,
    concat.Data = concat.Data,
    concat.block.Data = concat.block.Data
  )
  
  res$res.iter= rep(0,ncomp)
  
  
  # Criterion
  res$CRIT.h = list()
  res$CRIT = matrix(0, ncol=ncomp, nrow=1)
  colnames(res$CRIT) = paste("Dim", 1:ncomp, sep="")
  res$crit.group = matrix(0,ncol=M,nrow=ncomp)
  rownames(res$crit.group) = paste("Dim", 1:ncomp, sep="")
  colnames(res$crit.group) = paste("Group", 1:M, sep="")
  res$crit.block = matrix(0, ncol=K, nrow=ncomp)
  rownames(res$crit.block) = paste("Dim", 1:ncomp, sep="")
  colnames(res$crit.block) = paste("Block", 1:K, sep="")
  
  # omega
  res$omega = matrix(0,nrow=K, ncol=ncomp)
  colnames(res$omega) = paste("Dim", 1:ncomp, sep="")
  rownames(res$omega) = paste("Block", 1:K, sep="")
  
  
  # loadingds  
  res$block.common.loading = vector("list",K)
  for(k in 1:K){
    res$block.common.loading[[k]] = matrix(0, ncol=ncomp, nrow=P[k])
    rownames(res$block.common.loading[[k]]) = colnames(concat.Data[[k]])
    colnames(res$block.common.loading[[k]]) = paste("Dim", 1:ncomp, sep="")
  }
  
  res$block.group.loadings = vector("list",K)
  for(k in 1:K){
    res$block.group.loadings[[k]] = vector("list",M)
    for(m in 1:M){
      res$block.group.loadings[[k]][[m]] = matrix(0, ncol=ncomp, nrow=P[k])
    }
  }
  
  
  #scores
  res$block.scores = vector("list",K)
  for(k in 1:K){
    res$block.scores[[k]] = matrix(0, ncol=ncomp, nrow=N)
  }
  
  res$global.scores = matrix(0, ncol=ncomp, nrow=N)
  rownames(res$global.scores) = Group
  colnames(res$global.scores) = paste("Dim", 1:ncomp, sep="")
  

  res$block.group.scores = vector("list",K)
  for(k in 1:K){
    res$block.group.scores[[k]] = vector("list",M)
    for(m in 1:M){
      res$block.group.scores[[k]][[m]] = matrix(0, ncol=ncomp, nrow=n[m])
    }
  }
  
  res$block.scores = vector("list",K)
  for(k in 1:K){
    res$block.scores[[k]] = matrix(0, ncol=ncomp, nrow=N)
    rownames(res$block.scores[[k]]) = Group
    colnames(res$block.scores[[k]]) = paste("Dim", 1:ncomp, sep="")
  }
  

  # explined variance
  #============================================================================
  #      		            Iterative algorithm of mbmgPCA
  #============================================================================
  eps = 1e-7       # for iteration loops
  for(h in 1:ncomp){
    
    
    #------------------------------ initial value of tm for m=(1,...,M)
    initial.value = vector("list", niter)
    max.crit = rep(0,niter) 
    for(i in 1:niter){
      
      tm = vector("list", M)
      
      for(m in 1: M){
        tmm = rnorm(n[m]) 
        tm[[m]] = normv(tmm)
      }
      initial.value[[i]]$tm = tm
      max.crit[i] = select_initial_Iter(tm=initial.value[[i]]$tm, K.Data=K.Data, concat.Data=concat.Data, K=K, M=M, P=P, n=n, N=N)
    }
    selected = which(max.crit == max(max.crit))
    
    
    #------------------------------
    CRIT = 0
    iter = 0
    x    = 1.0
    res.it = list()
    res.it$tm = initial.value[[selected]]$tm
  
    
    save.global.t =matrix(0, ncol=1, nrow=N)
    while (x>eps){
      iter = iter + 1  
      res.it = iteraite_fun_Iter(tm=res.it$tm, K.Data=K.Data, concat.Data=concat.Data, K=K, M=M, P=P, n=n, N=N)
      save.global.t = cbind(save.global.t,res.it$global.t)
      #========================================================================
      #                          optimization criterion
      #========================================================================      
      #Global optimization criterion        
      crit = 0
      for(k in 1: K){
        for(m in 1: M){
          crit = crit +  t(res.it$tm[[m]]) %*%  K.Data[[k]][[m]]  %*%  res.it$load.block[[k]]  %*% t(res.it$load.block[[k]])  %*% t(K.Data[[k]][[m]])  %*% res.it$tm[[m]]   
        }
      }  
      CRIT[iter]=crit
      
      #block optimization criterion
      crit.block = c(rep(0,K))
      for(k in 1: K){
        crit.b=0
        for(m in 1: M){
          crit.b = crit.b+  t(res.it$tm[[m]]) %*%  K.Data[[k]][[m]] %*%  res.it$load.block[[k]] %*% t(res.it$load.block[[k]])  %*% t(K.Data[[k]][[m]])  %*% res.it$tm[[m]]   
        }
        crit.block[k] = crit.b / crit
      } 
      
      #group optimization criterion
      crit.group = c(rep(0,M))
      for(m in 1: M){
        critm = 0
        for(k in 1: K){
          critm = critm +  t(res.it$tm[[m]]) %*%  K.Data[[k]][[m]] %*%  res.it$load.block[[k]] %*% t(res.it$load.block[[k]])  %*% t(K.Data[[k]][[m]])  %*% res.it$tm[[m]]   
        }
        crit.group[m] = critm / crit
      }
      
      #---------------------
      if (iter>1){
        x = CRIT[iter] - CRIT[(iter-1)]
      }
    
        
    } # end of iteration
    
  
    #========================================================================
    #                                 saving results 
    #========================================================================
    res$omega[,h]= res.it$omega
    
    res$CRIT.h[[h]] = CRIT
    res$CRIT[h] = crit
    res$crit.block[h,] = crit.block
    res$crit.group[h,] = crit.group
    
    
    #loadings
    for(k in 1:K){
      res$block.common.loading[[k]][,h] = res.it$load.block[[k]]
    }
    
    for(k in 1:K){
      for(m in 1:M){
        res$block.group.loadings[[k]][[m]][,h] = res.it$block.group.loadings[[k]][[m]]
      }
    }
    
    
    #scores
    res$global.scores[,h] = res.it$global.t
    for(k in 1:K){
      res$block.scores[[k]][,h] = res.it$Tblocks[,k]
    }
    
    for(k in 1:K){
      for(m in 1:M){
        res$block.group.scores[[k]][[m]][,h] = res.it$block.group.scores[[k]][[m]]
      }
    }
    
    #========================================================================
    #                                 Deflation 
    #========================================================================
    for(k in 1:K){
      concat.Data[[k]] = deflation(concat.Data[[k]], as.matrix(res$global.scores[,h]))
      K.Data[[k]] = split( concat.Data[[k]], Group)
      
      for(m in 1:M){  
        K.Data[[k]][[m]] = matrix(K.Data[[k]][[m]], ncol=P[k])
        colnames(K.Data[[k]][[m]]) = colnames(res$K.Data[[k]][[m]])
      }
      
    }
    
  }  # END OF DIMENSION
  
  #========================================================================
  #                            explained variance
  #========================================================================
  
  # explained variance for each block and group X(b)_m
  res$cum.exp.var.block.group = vector("list",K)
  for(k in 1:K){ 
    res$cum.exp.var.block.group[[k]] = matrix(0, nrow=(M+1), ncol=ncomp)
  }
  
  if(M>1){
    Gnames = paste("Group", 1:M, sep="")
    for(k in 1:K){ 
      for(m in 1:M){ 
        for(h in 1:ncomp){ 
          proj    = project(res$block.group.scores[[k]][[m]][,1:h])
          Xkm.hat = proj %*%  res$K.Data[[k]][[m]]
          explvar = Trace(Xkm.hat) / Trace(K.Data[[k]][[m]])
          res$cum.exp.var.block.group[[k]][m,h]=explvar
        }
      }
      colnames(res$cum.exp.var.block.group[[k]]) = paste("Dim", 1:ncomp, sep="")
      rownames(res$cum.exp.var.block.group[[k]]) = c(Gnames, "Average")
      res$cum.exp.var.block.group[[k]][(M+1),] = apply(res$cum.exp.var.block.group[[k]][1:M,],2,mean)
    }
  }
  
  if(M==1){
    Gnames = paste("Group", 1:M, sep="")
    for(k in 1:K){ 
      for(m in 1:M){ 
        for(h in 1:ncomp){ 
          proj    = project(res$block.group.scores[[k]][[m]][,1:h])
          Xkm.hat = proj %*%  res$K.Data[[k]][[m]]
          explvar = Trace(Xkm.hat) / Trace(K.Data[[k]][[m]])
          res$cum.exp.var.block.group[[k]][m,h]=explvar
        }
      }
      colnames(res$cum.exp.var.block.group[[k]]) = paste("Dim", 1:ncomp, sep="")
      rownames(res$cum.exp.var.block.group[[k]]) = c(Gnames, "Average")
      res$cum.exp.var.block.group[[k]][(M+1),] = res$cum.exp.var.block.group[[k]][1:M,]
    }
  } 
  
  
  # explained variance for each block
  res$cum.expvar.block = matrix(0, nrow=K, ncol=ncomp)
  rownames(res$cum.expvar.block)=paste("Block", 1:K, sep=" ")
  colnames(res$cum.expvar.block)=paste("Dim", 1:ncomp, sep="")  
  for(k in 1:K){ 
    for(h in 1:ncomp){ 
      proj    = project(res$block.scores[[k]][,1:h])
      Xk.hat = proj %*%  res$concat.Data[[k]]
      explvar = sum(diag(Xk.hat %*% t(Xk.hat))) / Trace (res$concat.Data[[k]])
      res$cum.expvar.block[k,h] = explvar
    }
  }
  
  # global explained variance 
  res$global.expvar = matrix(0, nrow=1, ncol=ncomp)
  for(h in 1:ncomp){
    aa = sum(diag(t(res$concat.block.Data) %*% res$concat.block.Data))
    bb = c(t(res$global.scores[,h]) %*% res$global.scores[,h] )
    cc = t(res$global.scores[,h]) %*% res$concat.block.Data %*% t(res$concat.block.Data) %*% res$global.scores[,h]
    res$global.expvar[,h] = 100* cc/ (aa*bb)
    #res$global.expvar[,h] = cc/ (bb)
  }
  colnames(res$global.expvar) = paste("Dim", 1:ncomp, sep="")
  
  
  # similarity between block-group loadings and block common loadings
  res$similarity = list()
  for(k in 1:K){
    res$similarity[[k]] = matrix(0, nrow=M, ncol=ncomp)
    colnames(res$similarity[[k]]) = paste("Dim", 1:ncomp, sep="")
    rownames(res$similarity[[k]]) = levels(Group)
    
    for(m in 1:M){
      c = 0
      for(h in 1:ncomp){
        c = c +  abs(sum (res$block.group.loadings[[k]][[m]][,h] * res$block.common.loading[[k]][,h]))
        res$similarity[[k]][m,h] = c/h
      }
    }
  }
  
  #-------------------------------------
  # add class
  class(res) = "mbmgPCA"
  return(res)
}


#' @S3method print mbmgPCA
print.mbmgPCA <- function(x, ...)
{
  cat("\nmultiblock and multigroup PCA\n")
  cat(rep("-",43), sep="")
  cat("\n$omega                        ", "Weight of blocks")
  cat("\n$global.scores                ", "global.scores")
  cat("\n$block.common.loading         ", "block.common.loading")
  
  invisible(x)
}

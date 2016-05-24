#' @title Iteration function
#' 
#' @description 
#' To calculate outputs in each iteration
#' 
#' @param The parameters in the main functions
#' @return The results for each iteration
#' @export
#' @keywords internal
iteraite_fun_Iter <- function(tm, K.Data, concat.Data, K, M, P, n, N){
  #------------------- for saving results  
  block.group.scores = vector("list",K)  # t(b)m
  for(k in 1:K){
    block.group.scores[[k]] = vector("list",M)
    for(m in 1:M){
      block.group.scores[[k]][[m]] = matrix(0, ncol=1, nrow=n[m])
    }
  }
  block.group.loadings = vector("list", K)   # a(b)m
  for(k in 1:K){
    block.group.loadings[[k]] = vector("list",M)
    for(m in 1:M){
      block.group.loadings[[k]][[m]] = matrix(0, ncol=1, nrow=P[k])
    }
  }
  
#----------------------------------------------------------------
  # a(k) : common loadings for each block
  load.block = list()
  for(k in 1: K){
    summ = 0
    for(m in 1: M){
      cc = t(K.Data[[k]][[m]]) %*%  tm[[m]]  %*% t(tm[[m]])  %*%  K.Data[[k]][[m]]
      summ = summ + cc
    }
    load.block[[k]] = as.matrix(eigen(summ)$vectors[,1]) 
    rownames(load.block[[k]]) = colnames(K.Data[[k]])
  }
  
  # a(k)m: loadings for each block and each group
  for(k in 1: K){
    for(m in 1: M){
      block.group.loadings[[k]][[m]]= t(K.Data[[k]][[m]]) %*% tm[[m]]
      block.group.loadings[[k]][[m]] = normv(block.group.loadings[[k]][[m]])
    }
  }
  


  #---------------------------- block score
  Tblocks=matrix(0,nrow=N, ncol=K)
  for(k in 1:K){
    Tblocks[,k] = concat.Data[[k]] %*% load.block[[k]]
  }
 

  
  #block group score t(k)m
  for(k in 1: K){
    for(m in 1: M){
      block.group.scores[[k]][[m]] = K.Data[[k]][[m]] %*% load.block[[k]]
    }
  }
  
  # Tm=[t(1)m|...|t(K)m]
  TTm = list()
  for(m in 1: M){
    TTm[[m]] = matrix(0, ncol=K, nrow=n[m])
  }
  
  for(m in 1: M){
    for(k in 1: K){         
      TTm[[m]][,k] = K.Data[[k]][[m]] %*% load.block[[k]]
    }
  }
  

  #---------------------------- omega: eigen vector
  sum.omeg=0
  for(m in 1: M){
    smm = t(TTm[[m]]) %*% TTm[[m]] %*%  t(TTm[[m]]) %*% TTm[[m]]         
    sum.omeg = sum.omeg + smm
  }
  
  omega = eigen(sum.omeg)$vectors[,1]
  
  
  #global score: tm 
  global.t = 0
  for(m in 1: M){
    tm[[m]]  = TTm[[m]] %*% omega
    global.t = rbind(global.t,tm[[m]])
  }
  
  global.t = global.t[-1,]
  

res = list()

res$load.block = load.block
res$block.group.scores = block.group.scores
res$block.group.loadings = block.group.loadings
res$Tblocks = Tblocks


res$omega = omega
res$global.t = global.t
res$TTm = TTm
res$tm = tm

return(res)
}
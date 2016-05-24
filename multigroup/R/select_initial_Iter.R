#' @title Selection of initial values
#' 
#' @description 
#' To Select the initial values
#' 
#' @param X a numeric matrix
#' @param parametrs from iteraite.fun.Iter function
#' @return maximum of the optimal criterion 
#' @export
#' @keywords internal
select_initial_Iter = function(tm=tm, K.Data=K.Data, concat.Data=concat.Data, K=K, M=M, P=P, n=n, N=N){

#------------------------------
CRIT = 0
iter = 0
x    = 1.0
eps = 1e-7   
res.it = list()
res.it$tm = tm


while (x > eps){
  iter = iter + 1
  
  
  res.it = iteraite_fun_Iter(tm=res.it$tm, K.Data=K.Data, concat.Data=concat.Data, K=K, M=M, P=P, n=n, N=N)
  
  
  
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
  CRIT[iter] = crit
  #---------------------
  if (iter>1){
    x = CRIT[iter] - CRIT[(iter-1)]
  }

  } # end of iteration

return(max(CRIT))

}
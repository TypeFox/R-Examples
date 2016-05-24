MultiOptim <- function(retries, optim.fun, ...){
  best <- list(eval=Inf)
  for(i in 1:retries){
    res <- optim.fun(...)
    if(res$eval < best$eval){
      best <- res
    }
  }
  
  res
}

# Scalarizing function optimization ---------------------------------------

OptimizeScalMultiCMAES <- function(goal, n, optim.params){
  MultiOptim(optim.params$retries, OptimizeScalCMAES, 
                 goal,n, optim.params)
}

OptimizeScalCMAES <- function(goal,n, optim.params){  
  res <- cmaes::cma_es(runif(n,min=optim.params$par.min,max=optim.params$par.max),goal,lower=rep(optim.params$par.min,n),upper=rep(optim.params$par.max,n), control=optim.params$inner)                  
  
  list(eval=res$value, w.par=res$par[[1]], d=res$par[-1])
}


OptimizeScalMultiNM <- function(goal, n, optim.params){
  MultiOptim(optim.params$retries, OptimizeScalNM, 
             goal,n, optim.params$inner)
}

OptimizeScalNM <- function(goal,n, optim.params){
  res <- optim(runif(n,min=optim.params$par.min,max=optim.params$par.max),goal, control=optim.params$inner)  
  
  list(eval=res$value, w.par=res$par[[1]], d=res$par[-1])
}

OptimizeScalMultiDE <- function(goal, n, optim.params){
  MultiOptim(optim.params$retries, OptimizeScalDE, 
                 goal,n, optim.params)
}

OptimizeScalDE <- function(goal,n, optim.params){  
  res <- DEoptim::DEoptim(goal,lower=rep(optim.params$par.min,n),upper=rep(optim.params$par.max,n), optim.params$inner)                  
  
  list(eval=res$optim$bestval, w.par=res$optim$bestmem[[1]], d=res$optim$bestmem[-1])
}



# Activation function optimization ----------------------------------------

OptimizeActivMultiCMAES <- function(goal, optim.params){
  MultiOptim(optim.params$retries, OptimizeActivCMAES, goal,optim.params)
}

OptimizeActivCMAES<- function(goal, optim.params){  
  res <- cmaes::cma_es(runif(2,min=optim.params$par.min,max=optim.params$par.max),goal,lower=rep(optim.params$par.min,2),upper=rep(optim.params$par.max,2), control=optim.params$inner)                  
  
  list( eval=res$value, 
        a=res$par[[1]],
        b=res$par[[2]]
  )
}


OptimizeActivMultiNM <- function(goal, optim.params){
  MultiOptim(optim.params$retries, OptimizeActivNM, goal,optim.params)
}


OptimizeActivNM<- function(goal, optim.params){
  res <- optim(runif(2,min=optim.params$par.min,max=optim.params$par.max),goal, control=optim.params$inner)  
  
  list( eval=res$value, 
        a=res$par[[1]],
        b=res$par[[2]]
  )
}

OptimizeActivDE<- function(goal, optim.params){  
  res <- DEoptim::DEoptim(goal,lower=rep(optim.params$par.min,2),upper=rep(optim.params$par.max,2),optim.params$inner)                 
  
  list( eval=res$optim$bestval, 
        a=res$optim$bestmem[[1]],
        b=res$optim$bestmem[[2]]
  )
}

OptimizeActivMultiDE <- function(goal, optim.params){
  MultiOptim(optim.params$retries, OptimizeActivDE, goal,optim.params)
}
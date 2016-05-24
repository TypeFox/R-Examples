newton <-
function(x0, lb, ub, f, g, h, alpha = 0.25, beta = 0.5, max.iter = 100, tol = 1e-2){
  count = 0  
  x = x0
  
  repeat {
    count = count + 1
    
    ## newton's step
    delta = - g(x)/h(x)
    
    ## line search to pick the stepsize 
    size = 1
    while ( (x + size * delta <=0) || (log(f(x + size * delta)) > log(f(x) + alpha * size * g(x) * delta)) ) {
      size = beta * size
    }
    ## Update
    x.new = x + size * delta
    
    if ( count >= max.iter || abs(x-x.new)< tol || x.new > ub || x.new < lb ){
      
      if(count == max.iter) warning("Maximum number of iterations reached!")
      break
    }
    x = x.new
    
  }	
  
  return(list(solution = x, iter = count, stepSize = size))
}

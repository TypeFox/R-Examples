design.reg.polynom <- function(...){
  design.regression.polynom(...)
}

design.regression.polynom <- function(a, b, k, n){
  p <- n%/%(k+1)
  rem <- n%%(k+1)
    loc <- (a+b)/2+(b-a)/2*
    solve(deriv(legendre.polynomials(k)[[k+1]])*polynomial(c(1,0,-1)))
  repl <- rep(p,k+1)+c(rep(1,rem),rep(0,k+1-rem))
  d <- list(model="polynomial", order=k, locations=loc, replications=repl, interval=c(a,b))
  class(d) <- "design.regression"
  d
}

print.design.regression <- function(x, epl=6, ...){
  cat(paste("Design for a ", x$model, "model\n"))
  if(!is.null(x$order))
    cat(paste("Order:",x$order,"\n"))
  cat(paste("Interval: [",x$interval[1],",",x$interval[2],"]\n"))
  cat("\n")
  
  
  lines <- length(x$locations)%/%epl+1
  entries <- length(x$locations)
  for(i in 1:lines){
    if(lines==1)
      cat("( ")
    else if(i==1)
      cat("( ")
    else
      cat("  ")

    cat(sprintf("%1.3e",x$locations[((i-1)*epl+1):(min(entries,((i-1)*epl+1)+epl))]))
    
    if(lines==1)
      cat(" )\n( ")
    else if(i==1)
      cat("\n( ")
    else if(i==lines)
      cat(" )\n  ")
    else
      cat("\n  ")
    
    cat(sprintf("%9i",x$replications[((i-1)*epl+1):(min(entries,((i-1)*epl+1)+epl))]))
    if(lines==1)
      cat(" )\n")
    else if(i==1)
      cat("\n  ")
    else if(i==lines)
      cat(" )\n")
    else
      cat("\n  ")
    cat("\n")
  }
}


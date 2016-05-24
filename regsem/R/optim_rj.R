
optim_rj <- function(start,func,grad,hess,tol,max.iter){
  count = 0
  ret <- list()
  #max.iter = 1500
  convergence = 1
  vals <- rep(NA,max.iter)
  vals[1] <- 100
  new.pars <- matrix(NA,max.iter+1,length(start))
  new.pars[1,] <- start


  while(count < max.iter){
    count=count+1


    if(is.null(hess)==TRUE){
      if(count < 150){
        new.pars[count+1,] = new.pars[count,] - 0.5 * grad(new.pars[count,])
      }else if(count >= 150){
        new.pars[count+1,] = new.pars[count,] - 0.1 * grad(new.pars[count,])
      }else if(count >= 500){
        new.pars[count+1,] = new.pars[count,] - 0.01 * grad(new.pars[count,])
      }
    }else{
      hh = try(solve(hess(new.pars[count,])))

      if(inherits(hh, "try-error")){
        delta = - hh %*% (grad(new.pars[count,]))
        new.pars[count+1,] = new.pars[count,] + delta
      }else{
        if(count < 100){
          new.pars[count+1,] = new.pars[count,] - 0.5 * grad(new.pars[count,])
        }else if(count >= 200){
          new.pars[count+1,] = new.pars[count,] - 0.1 * grad(new.pars[count,])
        }else if(count >= 500){
          new.pars[count+1,] = new.pars[count,] - 0.01 * grad(new.pars[count,])
        }
      }
    }


    vals[count+1] = func(new.pars[count+1,])


    st.crit = try(abs(vals[count+1] - vals[count])<tol)
    if(inherits(st.crit, "try-error")){
      convergence=99
    }else if(is.na(st.crit)==TRUE){
      convergence=99
    }else{
      if(st.crit==TRUE){
        convergence = 0
        break
      }
    }
  }
  ret$iterations <- count
  ret$value <- vals[count+1]
  ret$pars <- new.pars[count+1,]
  ret$convergence <- convergence
  ret

}

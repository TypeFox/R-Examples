random.test.data <- function(nrows=10, ncols=10, vars.no=3, dep.var.dis="normal", xycoords=TRUE){
  
  obs.no <- nrows * ncols
  
  if (dep.var.dis=="normal") {
    dep <- runif(obs.no)
  } 
  if (dep.var.dis == "poisson"){
    dep <- rpois(obs.no, lambda = 7)
  }
  
  if (dep.var.dis == "zip"){
    dep <- rpois(obs.no, lambda = 7)
    idx <- sample(obs.no, 0.2*obs.no)
    dep[ idx ] <- 0
  }
  
  if (xycoords == TRUE){
    X <- rep(1:nrows, each=ncols) 
    Y <- rep(1:ncols, nrows)
  }
  
  vars<-matrix(data=NA, nrow=obs.no, ncol=vars.no-1)
  
  for (i in 1:vars.no-1) {
    vars[,i] <- runif(obs.no)
  }
  
  if (xycoords == TRUE){
      random.df <- data.frame(dep=dep, vars, X=X, Y=Y)
  } else {
      random.df <- data.frame(dep=dep, vars)
  }
  return(random.df)
}
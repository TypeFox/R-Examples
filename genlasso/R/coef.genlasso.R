coef.genlasso <- function(object, lambda, nlam, df, type=c("primal", "dual", "both"), ...) {
  type = type[[1]]
  if (!(type %in% c("primal", "dual", "both"))) {
    stop("Invalid type, must be \"primal\", \"dual\", or \"both\".")
  }
  
  if (missing(lambda) && missing(nlam) && missing(df)) {
    # Nothing was specified, use default lambdas
    if (type == "primal") return(list(beta=object$beta,lambda=object$lambda))
    if (type == "dual") return(list(u=object$u,lambda=object$lambda))
    if (type == "both") return(list(beta=object$beta,u=object$u,lambda=object$lambda))
  }
  else if (!missing(lambda) || !missing(nlam)) {
    # Check that they don't specify a lambda too small
    if (!missing(lambda)) {
      if (!is.numeric(lambda) || any(lambda < 0)) {
        stop("All specified lambda values must be nonnegative.")
      }
      if (min(lambda) < min(object$lambda) && !object$completepath) {
        stop(paste(sprintf("The path was truncated at %0.3f,",min(object$lambda)),
                   " and you asked for the solution at a smaller lambda.",sep=""))
      }
    }

    # Create nlam evenly spaced lambdas on the log scale
    else {
      m1 = max(object$lambda)
      m2 = min(object$lambda)
      if (m1==m2) {
        lambda = object$lambda
      }
      else {
        if (m2==0) m2 = min(object$lambda[object$lambda>m2])
        if (object$completepath) {
          # Over the first half of the path
          big = log(m1)
          sma = (log(m1)+log(m2))/2
          lambda = exp(seq(big,sma,length=nlam))
        }
        else {
          # Over the entire computed path
          big = log(m1)
          sma = log(m2)
          lambda = exp(seq(big,sma,length=nlam))
          # We do this because of numerical imprecision issues...
          lambda[nlam] = m2
        }
      }
    }
    
    # Sort the lambdas in decreasing order
    o = order(lambda,decreasing=TRUE)
    lambda = lambda[o]

    knots = object$lambda
    betas = object$beta
    us = object$u
    dfs = object$df
    
    # If the path is complete, append the solutions
    # at lambda=0
    if (object$completepath) {
      knots = c(knots,0)
      betas = cbind(betas,object$bls)
      us = cbind(us,rep(0,nrow(us)))
      dfs = cbind(dfs,dfs[length(dfs)]+1)
    }
    
    k = length(lambda)
    mat = matrix(rep(knots,each=k),nrow=k)
    b = lambda >= mat
    blo = max.col(b,ties.method="first")
    bhi = pmax(blo-1,1)

    i = bhi==blo
    p = numeric(k)
    p[i] = 0
    p[!i] = ((lambda-knots[blo])/(knots[bhi]-knots[blo]))[!i]
  
    if (type != "dual") {
      beta = t((1-p)*t(betas[,blo,drop=FALSE]) + p*t(betas[,bhi,drop=FALSE]))
      colnames(beta) = as.character(round(lambda,3))
    }
    if (type != "primal") {
      u = t((1-p)*t(us[,blo,drop=FALSE]) + p*t(us[,bhi,drop=FALSE]))
      colnames(u) = as.character(round(lambda,3))
    }
    df = dfs[blo]
    
    # Return in original order
    o = order(o)
    if (type == "primal") return(list(beta=beta[,o,drop=FALSE],lambda=lambda[o],df=df[o]))
    if (type == "dual") return(list(u=u[,o,drop=FALSE],lambda=lambda[o],df=df[o]))
    if (type == "both") return(list(beta=beta[,o,drop=FALSE],u=u[,o,drop=FALSE],lambda=lambda[o],df=df[o]))
  }
  else {
    # For each specified df value, find the largest df along the path
    # that is <= the specified value (often, this will be =). If there
    # is more than one such spot on the path, take the last spot
    a = matrix(object$df,length(df),length(object$df),byrow=TRUE)
    b = a <= df
    i = which(rowSums(b)!=0)
    b[b!=0] = a[b]
    b = b == apply(b,1,max)
    j = max.col(b,ties.method="last")
    j = j[i]
    dfs = a[i+(j-1)*nrow(a)]

    if (type != "dual") {
      beta = object$beta[,j,drop=FALSE]
      colnames(beta) = as.character(dfs)
    }
    if (type != "primal") {
      u = object$u[,j,drop=FALSE]
      colnames(u) = as.character(dfs)
    }
    lambdas = object$lambda[j]
    
    if (type == "primal") return(list(beta=beta,lambda=lambdas,df=dfs))
    if (type == "dual") return(list(u=u,lambda=lambdas,df=dfs))
    if (type == "both") return(list(beta=beta,u=u,lambda=lambdas,df=dfs))
  }
}


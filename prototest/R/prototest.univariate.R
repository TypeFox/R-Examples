prototest.univariate <-
function(x, y, type=c("ALR", "ELR", "MS", "F"), selected.col=NULL, lambda, mu=NULL, sigma=1, hr.iter=50000, hr.burn.in=5000, verbose=FALSE, tol=10^-8){
  if (!is.matrix (x)){stop("x should be a matrix")}
  
  ### sort out the data
  y = as.matrix (y)
  x = apply (x, 2, function(col){col-mean(col)})
  x = apply (x, 2, function(col){col/sqrt(sum(col^2))})


  #####################
  #### INPUT CHECK ####
  #####################
  
  # input dimensions
  if (length(y) != nrow(x)){stop("x and y should have the same number of rows")}
  n = length(y)
  p = ncol (x)
  y.bar = mean(y)
  
  
  # type
  type = type[1]
  if (!(type %in% c("ALR", "ELR", "MS", "F"))){stop("\"type\" should be one of \"ALR\", \"ELR\", \"MS\" or \"F\"")}
  
  # selected col
  if (!is.null(selected.col)){
    selected.col = sort (unique(selected.col))
    if (length(selected.col) == 0){stop("must select non-zero number of columns")}
    if (length(selected.col) > p | min(selected.col) <= 0 | max(selected.col) > p){stop("invalid columns selected")}
    if ((length(selected.col) + !is.null(mu)) > n){stop("too many columns selected")}
    if (type == "MS" & length(selected.col) != 1){stop("must select only one column for the marginal screening prototype method")}
  }
  
  
  
  #################
  #### TESTING ####
  #################
  
  y.ref = NULL
  if (is.null(selected.col)){ ##### SELECTIVE #####
                              
                              
    #### select columns #####
    if (verbose){print("Selecting columns and generating constraints...", quote=FALSE)}
    if (type == "MS"){ # select only one column -- maximum marginal correlation prototype
      Ab.obj = mc.selection.A.b(y=y, X=x, mu=mu)
    }else{ # select a bunch using the lasso
      Ab.obj = enet.selection.A.b(y=y, X=x, lambda=lambda, mu=mu)
    }
    selected.col = Ab.obj$which.col
    
    if (length(selected.col) == 0){ ## exit, with warning
      warning ("no columns selected, returning p-val = 1")
      return (list (ts=0, p.val=1, ts.ref=NULL))
    }
    
    
    
    #### extract constraints
    A = Ab.obj$A
    b = Ab.obj$b
    
    
    
    
    #### generate hit-and-run distribution?
    if (hr.iter > 0){ ##### HIT-AND-RUN ######
      
      ## generate hit-and-run samples 
      if (verbose){print ("Hit-and-run: generating samples...", quote=FALSE)}
      
      A.row.sum = apply (A, 1, sum)
      if (is.null(mu)){ # mu not specified, so condition on the sample average of the ys
        
        # condition on the value of y bar
        delta = A.row.sum*mean(y)
      
        # adjust constraint matrices
        A.tilde = sigma*(A - A.row.sum%*%t(rep(1, ncol(A)))/n)
        b.tilde = b - delta
        y.init = (y - y.bar)/sigma
      
        # generate samples
        y.hr = sigma*rcpp_generate_hit_and_run_samples (num_samples=hr.iter, burn_in=hr.burn.in, init_y = y.init, A=A.tilde, b=b.tilde)
        y.hr = apply (y.hr, 2, function(col){col-mean(col)})
        y.hr = y.hr + y.bar
      }else{
        
        # adjust constraint matrices
        A.tilde = sigma*A
        b.tilde = b - mu*A.row.sum
        y.init = (y - mu)/sigma
        
        # generate samples
        y.hr = sigma*rcpp_generate_hit_and_run_samples (num_samples=hr.iter, burn_in=hr.burn.in, init_y = y.init, A=A.tilde, b=b.tilde) + mu
      }
      
      ## compute test statistic replications and p-value
      if (verbose){print ("Hit-and-run: computing p-value...", quote=FALSE)}
      test.stat.obj = compute.test.statistic (x=x[, selected.col, drop=FALSE], y=cbind(y, y.hr), type=type, mu=mu, sigma=sigma, verbose=FALSE, tol=tol)
      ts = test.stat.obj$ts
      test.stat = ts[1]
      test.stat.ref = ts[-1]
      p.val = mean (test.stat <= test.stat.ref)
      y.ref = y.hr
      
    }else{ ##### NON-SAMPLING #####
      # should we exit? 
      if (is.null(mu) & (type == "ELR" | type == "ALR")){stop("non-sampling tests not implemented for ELR and ALR with unspecified mu: use hit-and-run reference distribution")}
      
      # compute the test statistic and p-value
      ts.p.obj = compute.selective.ts.and.p.val (y=y, x=x[,selected.col,drop=FALSE], type=type, A=A, b=b, mu=mu, sigma=sigma)
      test.stat = ts.p.obj$ts
      p.val = ts.p.obj$p.val
      ts.ref.dist = NULL
      
    }
  }else{ ##### NON-SELECTIVE -- just do the tests #####
    df1 = df2 = NULL
    M = length (selected.col)
    
    #### test statistic
    test.stat = compute.test.statistic (x=x[, selected.col, drop=FALSE], y=y, type=type, mu=mu, sigma=sigma, verbose=verbose, tol=tol)
    
    
    #### degrees of freedom, if required
    if (type == "F"){
      df1 = test.stat$df1
      df2 = test.stat$df2
    }
    else if (type == "MS"){
      df1 = test.stat$df
    }
    test.stat = test.stat$ts
    
    
    #### reference distrbution and p-value
    p.val = compute.non.selective.p.val(ts=test.stat, type=type, df1=df1, df2=df2, mu=mu, sigma=sigma, M=M)
  }
  
  out = list (ts=test.stat, p.val=p.val, selected.col=selected.col, y.ref=y.ref)
  class (out) = "prototest"
  return(out)
}

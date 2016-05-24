cv.glmlep <-
function(x, y, family = c("gaussian", "binomial"), lambda=NULL,
           lambda.min=ifelse(n<p,0.05,0.001), nlambda=100, lambda2=0,
           kappa=ifelse(n<p,0.1,0.05), pen.fac=rep(1,p), tol=1e-6, max.ite = 1e3, foldid, nfolds=5, cv.seed=100){
    
    family = match.arg(family)
    this.call = match.call()

    penalty.factor = pen.fac

    n <- dim(x)[1]
    p <- dim(x)[2] 
    if (length(y) != n) 
      stop(paste("number of observations in y (", length(y), ") not equal to the number of rows of x (", 
                 n, ")", sep = ""))
    
    if (missing(foldid)){
      set.seed(cv.seed)
      foldid = sample(rep(seq(nfolds), length = n))
    }
    
    ## generate lambda sequence
    if(!is.null(lambda)){
      if (any(lambda < 0)) 
        stop("lambdas should be non-negative")
      nlambda <- length(lambda)
      ls <-  lambda
    }  else {
      if (lambda.min >= 1) 
        stop("lambda.min should be less than 1")
      ls <- SetLambda(x, y, lambda.min, nlambda, penalty.factor)
    }
    
    beta0<-rep(0,p+1)
    ## give estimates for different lambda values.
    beta = switch(family, 
                  gaussian = sapply(1:length(ls), function(i) .C("gaulep",as.double(x),as.double(y),as.double(ls[i]),
                                                                 as.double(kappa),as.double(tol),as.integer(max.ite),
                                                                 as.integer(n),as.integer(p),as.double(beta0),DUP=FALSE)[[9]]),
                  binomial = sapply(1:length(ls), function(i) .C("binlep",as.double(x),as.double(y),as.double(ls[i]),
                                                                 as.double(kappa),as.double(tol),as.integer(max.ite),
                                                                 as.integer(n),as.integer(p),as.double(beta0),DUP=FALSE)[[9]]))
    
    
    ## using cross validation to choose lambda
    
    loss <- matrix(,nrow=nlambda,ncol=nfolds)
    for (j in 1:nfolds){
      ind <- foldid==j # test data index
      train.x <- x[-ind,]
      train.y <- matrix(y[-ind],ncol=1)
      test.x <- x[ind,]
      test.y <- matrix(y[ind],ncol=1)
      
      nj <- dim(train.x)[1]
      ## fit the model and estimate some type of criterion on the test data.
      betaj = switch(family, 
                    gaussian = sapply(1:nlambda, function(i) .C("gaulep",as.double(train.x),as.double(train.y),as.double(ls[i]),
                                                                   as.double(kappa),as.double(tol),as.integer(max.ite),
                                                                   as.integer(nj),as.integer(p),as.double(beta0),DUP=FALSE)[[9]]),
                    binomial = sapply(1:nlambda, function(i) .C("binlep",as.double(train.x),as.double(train.y),as.double(ls[i]),
                                                                   as.double(kappa),as.double(tol),as.integer(max.ite),
                                                                   as.integer(nj),as.integer(p),as.double(beta0),DUP=FALSE)[[9]]))
      loss[,j] <- sapply(1:nlambda, function(i)loglike(test.x, test.y, beta=betaj[,i], family))

    }
    mloss <- apply(loss,1,sum)
    ind_l <- which.min(mloss)
       
    df <- apply(beta[-1,], 2, function(x) sum(x!=0)) 
    
    beta_min <- beta[,ind_l]
    
    fit <- list(beta=beta, lambda=ls, df=df, loss=loss, lambda.min=ls[ind_l], beta.min=beta_min)
    fit$call <- this.call
    class(fit) = c(paste(substr(family,1,3),"lep",sep=""), "glmlep")
    return(fit)
  }

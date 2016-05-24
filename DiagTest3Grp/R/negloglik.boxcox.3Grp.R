##########################################################################################
####This function calculates the negative joint loglikelihood of 3 groups of diagnostic marker data, modified based on the function for 2 groups
####This function is called by "Youden3Grp.boxcox.fit" as optimization goal function
##########################################################################################
negloglik.boxcox.3Grp <- function (lambda.val, data1, xmat1=as.matrix(rep(1,length(data1))),data2,xmat2=as.matrix(rep(1,length(data2))),data3,xmat3=as.matrix(rep(1,length(data3))), lik.method = "ML")
  {

    .negloglik.boxcox <- function (lambda.val, data, xmat, lik.method = "ML") 
      {
      ####This function is called by boxcox.fit() to get the negative logliklihood at a fixed lambda for a group of diagnostic measurements
       ###For now, we only use the method of "ML". The use of the method of "RML" for parameter estimation (a type of ML but the nuisance para has no effect on the likelihood)
        ### allows unbiased estimate of var/covar parameters
        if (length(lambda.val) == 2) {
          data <- data + lambda.val[2]
          lambda <- lambda.val[1]
        }
        else lambda <- lambda.val
        lambda <- unname(lambda)    
        n <- length(data)
        beta.size <- ncol(xmat)
        if (isTRUE(all.equal(unname(lambda), 0))) 
          yt <- log(data)
        else yt <- ((data^lambda) - 1)/lambda
        beta <- solve(crossprod(xmat), crossprod(xmat, yt))
        ss <- sum((drop(yt) - drop(xmat %*% beta))^2)/n
        
        if (lik.method == "ML") 
         #neglik <- (n/2) * log(ss) - ((lambda - 1) * sum(log(data)))
          neglik <- (n/2) * (log(ss)+log(2*pi)+1) - ((lambda - 1) * sum(log(data)))
        if (lik.method == "RML") {
          xx <- crossprod(xmat)
          if (length(as.vector(xx)) == 1) 
            choldet <- 0.5 * log(xx)
          else choldet <- sum(log(diag(chol(xx))))
          neglik <- ((n - beta.size)/2) * log(ss) + choldet - ((lambda - 1) * sum(log(data)))
        }
        if (mode(neglik) != "numeric") 
          neglik <- Inf
        return(drop(neglik))
      }


    neglik1 <- .negloglik.boxcox (lambda.val, data1, xmat1, lik.method)
    neglik2 <- .negloglik.boxcox (lambda.val, data2, xmat2, lik.method)
    neglik3 <- .negloglik.boxcox (lambda.val, data3, xmat3, lik.method)
    
    neglik <- neglik1+neglik2+neglik3
    return(neglik)
  }

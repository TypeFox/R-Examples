lambda.reg.bayes <- function(object, columns, ret.mcmc = TRUE){

  if (class(object) != "eiRegBayes")
    stop("'object' must be output from 'ei.reg.bayes'")
  if (missing(columns) | length(columns) < 2)
    stop("'columns' requires at least two column names")

  lambda.out <- array(NA, dim=c(length(rownames(object$draws)),
                       length(columns), dim(object$draws)[3]))
  rownames(lambda.out) <- rownames(object$draws)
  colnames(lambda.out) <- columns

  for(i in columns){
    lambda.out[,i,] <-
      object$draws[,i,]/apply(object$draws[,columns,],c(1,3),sum)
  }

  if (ret.mcmc){
    lambda.out <- t(matrix(lambda.out,
                           nrow(lambda.out)*ncol(lambda.out),
                           dim(lambda.out)[3]))
    colnames(lambda.out) <- apply(expand.grid(rownames(object$draws),
                                              columns)[,1:2], 1,
                                  paste, collapse=".") 
    lambda.out <- coda::mcmc(lambda.out)
  }
  
  class(lambda.out) <- c("lambdaRegBayes", class(lambda.out))
  lambda.out
}

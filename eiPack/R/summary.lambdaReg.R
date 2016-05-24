summary.lambdaReg <- function(object, ...){ 
  nidx <- apply(expand.grid(dimnames(object$lambda.out)), 1, paste, collapse = ".")
  object$table <- cbind(c(object$lambda.out), c(object$se), 
                       c(object$lambda.out)/c(object$se))
  colnames(object$table) <- c("Estimate", "Std. Error", "t-stat")
  rown <- array(NA, nrow(object$table))
  for(i in 1:ncol(object$lambda.out)){
    rown[((i-1)*nrow(object$lambda.out)+1):(i*nrow(object$lambda.out))] <-
      paste("Proportion", rownames(object$lambda.out), "in",
            colnames(object$lambda.out)[i])
  }
  rownames(object$table) <- rown
  object$lambda.out <- object$table
  object$table <- NULL
  object$lambda.out
}

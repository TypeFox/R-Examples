
summary.pscore <- function(object,
                           ...){

  res <- list(sum=summary(object$model.pscore))

  class(res) <- "summary.pscore"
  
  res


}

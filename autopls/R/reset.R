reset <- function (object, verbose = TRUE)
{

  iter <- object$metapls$autopls.iter
  lv <- object$metapls$autopls.lv
  
  ## Reconstruct first part
  new.object <- object$iterations [[iter]]
  att <- attributes (new.object$model)
  X <- object$metapls$X
  subX <- X [,new.object$predictors]
  Y <- object$metapls$Y
  new.object$model <- list (Y, subX)
  attributes (new.object$model) <- att
  names (new.object$model) <- c('Y', 'Xsub')
  
  ## Construct $metapls object
  part2 <- object$metapls
  part2$current.iter = iter
  part2$current.lv = lv
  
  ## Get iterations
  part3 <- object$iterations
  
  ##Put things together
  result <- new.object
  result [[length (result) + 1]] <- part2
  names (result) [length (result)] <- 'metapls'
  result [[length (result) + 1]] <- part3
  names (result) [length (result)] <- 'iterations'
  
  class (result) <- 'autopls'

  ## Reporting
  if (verbose) print (result)

  ## Output
  invisible (result)
}




  
  



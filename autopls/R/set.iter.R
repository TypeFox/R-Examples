set.iter <- function (object, iteration, verbose = TRUE)
{
  if (iteration > length (object$metapls$lv.history)) 
    stop (paste ('Maximum:', length (object$metapls$lv.history)))

  ## Model corresponding to the selected iteration
  new.object <- object$iterations [[iteration]]

  ## Reconstruct $model part
  att <- attributes (new.object$model)  
  X <- object$metapls$X
  subX <- X [,new.object$predictors]  
  Y <- object$metapls$Y    
  new.object$model <- list (Y, subX)
  attributes (new.object$model) <- att
  names (new.object$model) <- c('Y', 'subX')
  
  ## Set lv to the originally selected number
  current.lv <- object$metapls$lv.history [iteration]  
  names (current.lv) <- NULL
  
  part2 <- object$metapls
  part2$current.iter = as.integer (iteration)
  part2$current.lv = current.lv

  ## Get iteration data
  part3 <- object$iterations

  ## Put things together
  result <- new.object
  result [[length (result) + 1]] <- part2
  names (result) [length (result)] <- 'metapls'
  result [[length (result) + 1]] <- part3
  names (result) [length (result)] <- 'iterations'

  ## Output
  class (result) <- 'autopls'
  if (verbose) print (result)
  
  invisible (result)
}

npmixapply <- function(object, FUN, ...)  {
  ### object is of class "npmix"
  ### FUN is some vectorized function
  
  ### e.g., for posterior mean, FUN should be FUN = function(x) { x }
  
  switch(object$family$family,
         gaussian={
             tmp <- PostProbNorm(x = object$data[,1],std_err = object$data[,2],
                               support = object$support,mix.prop = object$mix.prop)
             PP <- tmp$postprobs
             ans <- PP%*%FUN(object$support, ...)
         },
         poisson={
             tmp <- PostProbPois(x = object$data[,1],eta = object$data[,2],
                               support = object$support,mix.prop = object$mix.prop)
             PP <- tmp$postprobs
             ans <- PP%*%FUN(object$support, ...)
         },
         binomial={
             tmp <- PostProbBinomial(x = object$data[,1],ntrials = object$data[,2],
                               support = object$support,mix.prop = object$mix.prop)
             PP <- tmp$postprobs
             ans <- PP%*%FUN(object$support, ...)
         },
         tdist={
             tmp <- PostProbT(x = object$data[,1],std_err = object$data[,2],
                              df = object$family$df, support = object$support,
                              mix.prop = object$mix.prop)
             PP <- tmp$postprobs
             ans <- PP%*%FUN(object$support, ...)
         }
  )
  return(ans)
}

## ####################################################################
## A function for checking the balance in covariate distribution
## between treatment groups if class of object is 'matched.data.frames'
## ####################################################################

ps.balance.matched.data.frames <- function(object,
                                           ...) 
{

  data <- rbind(as.data.frame(object$data[1]),
                as.data.frame(object$data[2]))

  data.matched <- rbind(as.data.frame(object$data.matched[1]),
                        as.data.frame(object$data.matched[2]))

  match.index <- c(object$match.index[[1]],
                   object$match.index[[2]])
   
  data <- list(data              = data,
               data.matched      = data.matched,
               match.index       = match.index,
               matched.by        = object$matched.by,
               name.match.index  = object$name.match.index,
               match.parameters  = object$match.parameters)
  
  ps.balance.matched.data.frame(object = data,
                                ...)
  
}

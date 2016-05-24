
ps.estimate.matched.data.frames <- function(object,
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
   
  ps.estimate.matched.data.frame(object = data,
                                 ...)
  
}



#ps.estimate.matched.data.frames <- function(object,
#                                            resp,
#                                            treat         = NULL,
#                                            stratum.index = NULL,
#                                            match.index   = NULL,
#                                            adj           = NULL,
#                                            weights       = "rr",
#                                            family        = "gaussian",
#                                            regr          = NULL,
#                                            ...) 
#{

#  data <- rbind(as.data.frame(object$data[1]),
#                as.data.frame(object$data[2]))

#  data.matched <- rbind(as.data.frame(object$data.matched[1]),
#                        as.data.frame(object$data.matched[2]))

#  match.index1 <- c(object$match.index[[1]],
#                    object$match.index[[2]])

#  treat1 <- c(object$treat[[1]],
#              object$treat[[2]])
   
#  data <- list(data             = data,
#               data.matched     = data.matched,
#               match.index      = match.index1,
#               name.match.index = object$name.match.index,
#               match.parameters = object$match.parameters,
#               treat            = treat1,
#               name.treat       = object$name.treat)

# # ps.estimate.matched.data.frame(object = data,
##                                 resp,
##                                 treat         = NULL,
##                                 stratum.index = NULL,
##                                 match.index   = NULL,
##                                 adj           = NULL,
##                                 weights       = "rr",
##                                 family        = "gaussian",
##                                 regr          = NULL,
##                                 ...)

  
#}

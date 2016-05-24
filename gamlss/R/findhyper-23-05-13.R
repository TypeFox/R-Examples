# last modification 23-05-13 DS
find.hyper <-function ( model = NULL, 
                   parameters = NULL, 
                        other = NULL, 
                            k = 2,  
                        steps = c(0.1), 
                        lower = -Inf, 
                        upper = Inf, 
                       method = "L-BFGS-B", 
                               ... )# 
{   
# ------------------------------------------------------------------------------
# local function
#-------------------------------------------------------------------------------
  fn <- function(p, pen=k)
    { 
    if (!is.null(other)) 
    {
      newXvar <- if (isDataIn) with(TheDataFrame, eval(other))
      else eval(other)
      assign( as.character(other[[2]]) , newXvar, envir=envFindHyp)
    }
       mod.1 <- if (isDataIn)  with(TheDataFrame, eval(model))
             else eval(model)
        call <<- mod.1$call
     call$mu.start <<- fitted(mod.1,"mu")
     if ("sigma"%in%mod.1$parameters)  call$sigma.start <<- fitted(mod.1,"sigma")
     if (   "nu"%in%mod.1$parameters)  call$nu.start    <<- fitted(mod.1,"nu")
     if (  "tau"%in%mod.1$parameters)  call$tau.start   <<- fitted(mod.1,"tau")
     cat("par",p,"crit=",IC(mod.1,pen),"with pen=",pen,"\n")
      IC(mod.1,pen)
    }
#-------------------------------------------------------------------------------
 if(is.null(model)) stop("you have not defined the model")
 if(is.null(parameters)) stop("you have not define the starting values of the parameters")
  lp <- length(parameters)
 if (lp==length(steps)) ndeps <- steps else ndeps <-rep(steps[1],lp)
# env <- attach(NULL, name="Save_Params_Env")
# assign("p", parameters, envir=env)
#if (!is.null(model$data)) 
#  {
#      attach(eval(model$data), name="The_Data_Env")
#      on.exit(detach("The_Data_Env"))
#      assign("p", parameters, envir=as.environment("The_Data_Env")) 
#   }
#else
#   {
#      env <- attach(NULL, name="The_Data_Env")
#      assign("p", parameters, envir=env)   
#     on.exit(detach("The_Data_Env")) 
#   }
#  assign("p", parameters, envir=globalenv())
#  rm(p 
#  assign("objectAll", object, envir=env)
#  
#  if (!is.null()) 
#  { 
#       attach(eval(model$data))
#    #attach(eval(substitute(model$data))
#    on.exit(detach(eval(model$data)))
#  }
    envFindHyp <- environment()
      isDataIn <- !is.null(model$data)
if (isDataIn) TheDataFrame <- eval(model$data)
            p <- parameters
  if (!is.null(other)) 
    {
    newXvar  <- if (isDataIn) with(TheDataFrame, eval(other))
                 else          eval(other)
    assign( as.character(other[[2]]) , newXvar, envir=envFindHyp)
     }
          mod.1 <- if (isDataIn) with(TheDataFrame, eval(model))
                   else eval(model)
           call <- mod.1$call
  call$mu.start <- fitted(mod.1,"mu")
 if ("sigma"%in%mod.1$parameters)  call$sigma.start <- fitted(mod.1,"sigma")
 if (   "nu"%in%mod.1$parameters)  call$nu.start    <- fitted(mod.1,"nu")
 if (  "tau"%in%mod.1$parameters)  call$tau.start   <- fitted(mod.1,"tau")
  o2 <- optim(p, fn, lower = lower, upper = upper, 
            method = method, control = list(ndeps = ndeps , ...))
 # on.exit(rm(p,envir=.GlobalEnv),add=TRUE)
  o2
}

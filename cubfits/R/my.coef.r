### Provide number of coefficients for the given model.

### Get the constant according to the options.
get.my.ncoef <- function(model, assign.Env = TRUE){
  if(model[1] == "rocnsef"){
    ret <- 3
  } else if(model[1] == "roc"){
    ret <- 2
  } else if(model[1] == "nsef"){
    ret <- 2
  } else{
    stop("model is not found.")
  }

  if(assign.Env){
    assign("my.ncoef", ret, envir = .cubfitsEnv)
  }
  ret
} # End of get.my.ncoef().

get.my.coefnames <- function(model, assign.Env = TRUE, as.delta.eta = T){
  pName <- "t"
  if(as.delta.eta){pName <- "eta"}
  
  if(model[1] == "rocnsef"){
    ret <- c("log.mu", paste("Delta", pName, sep="."), "omega")
  } else if(model[1] == "roc"){
    ret <- c("log.mu", paste("Delta", pName, sep="."))
  } else if(model[1] == "nsef"){
    ret <- c("log.mu", "omega")
  } else{
    stop("model is not found.")
  }

  if(assign.Env){
    assign("my.coefnames", ret, envir = .cubfitsEnv)
  }
  ret
} # End of get.my.coefnames().

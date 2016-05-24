# restructure the coefficients (outer loadings and path coefficients)
coef.sempls <- function(object, ...){
  blocks <- object$model$blocks
  model <- object$model
  latent <- model$latent           # names of the latent variables
  strucmod <- model$strucmod
  pC <- object$path_coefficients
  crossL <- object$cross_loadings
  W <- object$outer_weights

  arrows <- NULL
  coef_names <- NULL
  estimates <- NULL

  ## measurement model
  # reflectice: mode A
  fooA <- function(latent, blocks){
    paste(rep(latent , length(blocks[[latent]])),
          " -> ", blocks[[latent]], sep="")
  }
  # formative: mode B
  fooB <- function(latent, blocks){
    paste(blocks[[latent]], 
          " -> ", rep(latent , length(blocks[[latent]])), sep="")
  }

  ## iterate over all blocks
  for(i in 1:length(blocks)){
    if(attr(blocks[[i]], "mode")=="A"){
      arrows <- append(arrows, fooA(names(blocks)[i], blocks))
      # outer loadings for Mode 'A' (reflective)
      estimates <- append(estimates, crossL[blocks[[i]], names(blocks)[i]])
      coef_names <- append(coef_names, paste("lam_", i, "_", 1:length(blocks[[i]]), sep=""))
    }
    if(attr(blocks[[i]], "mode")=="B"){
      arrows <- append(arrows, fooB(names(blocks)[i], blocks))
      # outer weights for Mode 'B' (formative)
      estimates <- append(estimates, W[blocks[[i]], names(blocks)[i]])
      coef_names <- append(coef_names, paste("gam_", i, "_", 1:length(blocks[[i]]), sep=""))
    }
  }
  

  ## structural model
  foo <- function(strucmod){
    paste(strucmod[,1], " -> ", strucmod[,2], sep="")
  }
  
  arrows <- append(arrows, foo(strucmod))
  indx <- which(pC != 0, arr.ind=TRUE)
  coef_names <- append(coef_names, paste("beta_", indx[,1], "_", indx[,2], sep=""))
  estimates <- append(estimates, pC[pC != 0])
  
  coefficients <- data.frame(Path=arrows, Estimate=estimates)
  rownames(coefficients) <- coef_names
  return(coefficients)
}

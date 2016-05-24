kernCreate <-
function(x, kernType, kernOptions=NULL) {
  if ( is.list(x) ) {
    dim <- array()
    for ( i in 1:length(x) ) {
      dim[i] <- dim(as.matrix(x[[i]]))[2]
      if ( (dim[i] == 1) & (dim(as.matrix(x[[i]]))[1] == 1) )
        dim[i] <- x[[i]]
    }
  } else {
    dim <- dim(as.matrix(x))[2]
    if ( (dim == 1) & (dim(as.matrix(x))[1] == 1) )
      dim <- x
  }

  if ( is.list(kernType) && kernType$type == "parametric" ) {
    kernOptions <- kernType$options
    kernType <- kernType$realType
  }

  if ( is.list(kernType) && ("options" %in% names(kernType)) ) {
    kernOptions <- kernType$options
  }
  
  if ( is.list(kernType) && ("complete" %in% names(kernType)) ) {
    if ( kernType$complete == 1 ) {
      kern <- kernType
    }
    
  } else if ( is.list(kernType) ) {
    
    kern <- list(inputDimension=dim, type=kernType$type)

    if (!is.null(kernOptions))
      kern$options <- kernOptions
    
    start <- 1    
    
    if ( kern$type == "multi" ) {
      for ( i in start:length(kernType$comp) ) {
        if ( is.list(kernType$comp) ) {
          iType <- kernType$comp[[i]]
        } else {
          iType <- kernType$comp[i]
        }
        
        if ( is.list(x) ) {
          kern$comp[[i-start+1]] <- kernCreate(x[[i-start+1]], iType)
          kern$diagBlockDim[i-start+1] <- dim(as.array(x[[i-start+1]]))[1]
        } else {
          kern$comp[[i-start+1]] <- kernCreate(x, iType)
        }

        kern$comp[[i-start+1]]$index = array()
      }
      
    } else if ( kern$type %in% c("cmpnd", "tensor", "translate",
                                 "selproj") )  {
      for ( i in start:length(kernType$comp) ) {
        if ( is.list(kernType$comp) ) {
          iType <- kernType$comp[[i]]
        } else {
          iType <- kernType$comp[i]
        }
        
        if (kern$type == "selproj") {
          if ( (dim(as.matrix(x))[2] == 1) && (dim(as.matrix(x))[1] == 1) )
            x_proj <- x-1
          else
            x_proj <- x[,-1]
            
          kern$comp[[i-start+1]] <- kernCreate(x_proj, iType)
        } else {
          kern$comp[[i-start+1]] <- kernCreate(x, iType)
        }
        kern$comp[[i-start+1]]$index = array()
      }
      
    } else if ( kern$type == "exp" ) {
      ## need double check
      if ( start == length(kernType$comp) ) {
        kern$argument <- kernCreate(x, kernType$comp[start])
      } else {
        kern$argument <- kernCreate(x, kernType$comp[start:length(kernType$comp)])
      }
    }

    kern <- kernParamInit(kern)

  } else {
    kern <- list(type=kernType, inputDimension=dim)

    if (!is.null(kernOptions))
      kern$options <- kernOptions

    kern <- kernParamInit(kern)
  }

  kern$Kstore <- matrix()
  kern$diagK <- matrix()      

#   if (!is.null(kernOptions) && 'priors' %in% names(kernOptions)) {
#     kern$priors <- list()
#     for (k in seq_along(kernOptions$prior))
#       kern$priors[[k]] <- priorCreate(kernOptions$prior[[k]])
#   }
  
  return (kern)
  
}

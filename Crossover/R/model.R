# ppp=proportionality parameter
linkMatrix <- function(model, v, ppp=0.5, placebos=1) {
  model <- getModelNr(model)
  if(missing(v)) stop("Please specify number of treatments")
  mI <- diag(v)
  m1 <- matrix(1,v,1)
  m0 <- matrix(0, v, v)
  if(model=="Standard additive model" || model==1) { #TODO Convert the string comparisons to comments.
    return(rbind(cbind(mI, m0), 
                 cbind(kronecker(m1,mI),kronecker(mI,m1))))
  }
  if(model=="Second-order carry-over effects" || model==8) {
    return(rbind(cbind(mI,m0,m0),
                 cbind(kronecker(m1,mI),kronecker(mI,m1),matrix(0,v*v,v)),
                 cbind(kronecker(kronecker(m1,m1),mI), kronecker(kronecker(m1,mI),m1),kronecker(kronecker(mI,m1),m1))))
  }
  if(model=="Full set of interactions" || model==7) {
    M <- matrix(0, sum(1:v)*2, v*v)    
    for (j in (v+1):(sum(1:v)*2)) {
      jv <- (j-1)%/%v
      M[j, v*(j-1-v*jv)+jv] <- 1
    }
    return(cbind(linkMatrix("Standard additive model", v), M))
  }
  if(model=="Self-adjacency model" || model==2) {
    M <- cbind(linkMatrix("Standard additive model", v), matrix(0, sum(1:v)*2, v))
    for (j in (v+1):(sum(1:v)*2)) {
      jv <- (j-1)%/%v      
      if (jv==j-v*jv) {
        M[j,v+jv] <- 0
        M[j,2*v+jv] <- 1
      }
    }
    return(M)
  }
  if(model=="Placebo model" || model==4) {
    M <- matrix(0, sum(1:v)*2, 2*v)     
    for (j in 1:(sum(1:v)*2)) {
      jv <- (j-1)%/%v
      M[j,j-v*jv] <- 1
      if (j>v*(placebos+1)) {
        M[j,v+jv] <- 1
      }
    }
    return(M)
  }
  if(model=="No carry-over into self model" || model==5) {
    M <- linkMatrix("Standard additive model", v)
    for (j in (v+1):(sum(1:v)*2)) {
      jv <- (j-1)%/%v
      if (jv==j-v*jv) {
        M[j,v+jv] <- 0
      }
    }
    return(M)
  }
  if(model=="Treatment decay model" || model==6) {
    M <- matrix(0, sum(1:v)*2, 2*v)     
    for (j in 1:(sum(1:v)*2)) {
      jv <- (j-1)%/%v
      M[j,j-v*jv] <- 1
      if (j>v && jv==j-v*jv) {
        M[j,v+jv] <- -1
      }
    }
    return(M)
  }
  if(model=="Proportionality model" || model==3) {
    M <- matrix(0, sum(1:v)*2, v) 
    M[1:v,1:v] <- diag(v)
    for (j in (v+1):(sum(1:v)*2)) {
      jv <- (j-1)%/%v
      if (jv==j-v*jv) {
        M[j, jv] <- 1+ppp
      } else {
        M[j, j-v*jv] <- 1
        M[j, jv] <- ppp
      }
    }
    return(M)
  }
  if (model=="No carry-over effects" || model==9) {
    #M <- matrix(0, sum(1:v)*2, v) 
    #M[1:v,1:v] <- diag(v)
    #for (j in (v+1):(sum(1:v)*2)) {
    #  jv <- (j-1)%/%v      
    #  M[j, j-v*jv] <- 1      
    #}
    #return(M)
    return(diag(v))
  }  
  stop(paste("Sorry model \"",model,"\" is not known.", sep=""))
}

models <- c("Standard additive model",
            "Self-adjacency model",
            "Proportionality model",
            "Placebo model",
            "No carry-over into self model",
            "Treatment decay model",
            "Full set of interactions",
            "Second-order carry-over effects",
            "No carry-over effects")
#"No carry-over effects")

#' Get the number or character string specifying the model
#' 
#' Get the number or character string specifying the model
#'
#' @param model Number or character string specifying the model
#' @param type Eiher \code{"numeric"} or \code{"character"}. If numeric the number of the model will be returned.
#' Otherwise the character string description of the model.
#' @return Either number or character string specifying the model.
#' @examples
#' Crossover:::getModelNr("Self-adjacency model")==Crossover:::getModelNr(2)
#' "Self-adjacency model"==Crossover:::getModelNr(2, type="character")
#' Crossover:::getModelNr("Self-adjacency model")==2
getModelNr <- function(model, type="numeric") {
  if (type!="numeric") {
    if (type=="character") {
      model <- models[getModelNr(model)]
      return(model)
    } else {
      stop("Parameter type must be either \"numeric\" or \"character\".")
    }
  }
  if (is.numeric(model)) {
    if (model %in% 1:9) {
      return(model)
    } else {
      if (model==0) return(9) 
      stop("Model must be number between 1 and 9.")
    }
  }
  modelNr <- which(models==model)
  if (length(modelNr)==0) stop(paste("Unknown model \"", model ,"\".", sep=""))
  return(modelNr)
}

#' Create a row column design 
#'
#' @param X cross-over design
#' @param v number of treatments
#' @param model String or number describing the model. See \code{\link{getModelNr}}.
#' @return A row-column design (as matrix - but not the design matrix).
#' @seealso \code{\link{rcdMatrix}} gives the row-column design matrix.
#' @examples
#' # TODO
rcd <- function(X, v, model) {
  model <- getModelNr(model)
  return(.Call( "rcd2R", X, v, model, PACKAGE = "Crossover" ))
}

#' Create the design matrix for a given row column design 
#'
#' @param X row-column design
#' @param v number of treatments
#' @param model String or number describing the model. See \code{\link{getModelNr}}.
#' @return The design matrix for a row-column design.
#' @seealso \code{\link{rcd}} gives the row-column design to a given crossover design.
#' @examples
#' # TODO
rcdMatrix <- function(X, v, model) {
  if (length(levels(as.factor(X)))<=v && model !=9) warning("It looks like you called rcdMatrix with a crossover design,\nbut you should provide the row-column design.")
  model <- getModelNr(model)
  return(.Call( "rcdMatrix2R", X, v, model, PACKAGE = "Crossover" ))
}

infMatrix <- function(X, v, model, xtx=FALSE) {
  model <- getModelNr(model)  
  return(.Call( "infMatrix2R", X, v, model, xtx, PACKAGE = "Crossover" ))
}

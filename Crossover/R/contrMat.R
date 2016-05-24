#' Create the design matrix, variance-covariance matrix, the variance of each
#' pairwise comparison and the efficicency of each pairwise comparison for a
#' cross-over design
#' 
#' Function to read in a cross-over design and create the design matrix X, 
#' the variance of each pairwise comparison and the efficicency of each pairwise comparison.
#' 
#' See the vignette of this package for further details.
#' 
#' @param type Type of contrast. A character vector containing the following: "Dunnett", "Tukey", "none".
#' If the length is 1, this contrast is only applied for the treatment effects and for carry-over effects a "Tukey" contrast is used.
#' Otherwise the specified contrasts are used, see also the examples.
#' @param v Number of treatments
#' @param model Model - one of the following: 1) "Standard additive model",
#' 2) "Second-order carry-over effects", 3) "Full set of interactions",
#' 4) "Self-adjacency model", 5) "Placebo model", 6) "No carry-over into self
#' model", 7) "Treatment decay model", 8) "Proportionality model", 9) "No carry-over effects". 
#' Can be specified as number or as character string.
#' @param eff.factor Weight applied to the different sub contrast matrices. A warning is given if it does not sum up to one. See examples.
#' @return A contrast matrix
#' @author Kornelius Rohmeyer \email{rohmeyer@@small-projects.de}
#' @keywords misc
#' @examples
#' 
#' contrMat2("Tukey", v=3, model=1)
#' contrMat2("Dunnett", v=3, model=1)
#' contrMat2(c("Dunnett", "Dunnett"), v=3, model=1)
#' contrMat2(c("Dunnett", "none"), v=3, model=1)
#' contrMat2(c("Dunnett", "none", "none"), v=3, model=8)
#' contrMat2("Dunnett", v=3, model=1, eff.factor=c(0.9, 0.1))
#' contrMat2("Dunnett", v=3, model=8, eff.factor=c(0.5, 0.3, 0.2))
#' 
#' @export contrMat2
contrMat2 <- function(type, v, model, eff.factor=rep(1, length(parameterCount(model, v)))) {
  model <- getModelNr(model)
  if (length(eff.factor)!=length(parameterCount(model, v))) {
    stop(paste("For model ", model, " the paramter eff.factor must have length ", length(parameterCount(model, v)), ".", sep=""))
  }
  if (!isTRUE(all.equal(1, sum(eff.factor), check.attributes = FALSE))) {
    warning("Parameter eff.factor does not sum up to 1.")
  }
  if (length(type)!=length(parameterCount(model, v)) && length(type)==1) {
    type <- c(type, rep("Tukey", length(parameterCount(model, v))-1))
  }
  if (length(type)!=length(parameterCount(model, v))) {
    stop(paste("For model ", model, " the paramter type must have length ", length(parameterCount(model, v)), ".", sep=""))
  }
  if (all(type %in% c("Dunnett", "Tukey", "none"))) {
    Csub <- contrMat(n=rep(1, v), type=type[1])
    class(Csub) <- "matrix"
    C <- appendZeroColumns(Csub, model=model, v)
    if (length(eff.factor)<=1 || all(eff.factor[-1]==0) || model %in% c(3,9)) {
      return(C)
    }
    if (length(type)>1) {      
      if (type[2]=="none") {
        Csub2 <- matrix(0, 0, v)
      } else {
        Csub2 <- contrMat(n=rep(1, v), type=type[2])
        class(Csub2) <- "matrix"
      }      
      m <- matrix(0, dim(Csub2)[1], dim(Csub)[2])
    }
    if (model %in% c(1, 4, 5, 6)) { # v+v parameters
      C <- rbind(C*eff.factor[1], cbind(m, Csub2)*eff.factor[2])
    } else if (model %in% c(2, 8)) { # v+v+v parameters
      
      if (type[3]=="none") {
        Csub3 <- matrix(0, 0, v)
      } else {
        Csub3 <- contrMat(n=rep(1, v), type=type[3])
        class(Csub3) <- "matrix"
      }
      m <- matrix(0, dim(Csub3)[1], dim(Csub)[2])
      
      C <- rbind(C*eff.factor[1], cbind(m, Csub2, matrix(0,dim(Csub2)[1],v))*eff.factor[2])      
      C <- rbind(C, cbind(m, matrix(0,dim(Csub2)[1],v), Csub3)*eff.factor[3])
    } else if (model %in% c(7)) { # Full set of interactions v+v+v^2
      warning("Full set of interactions are not yet implemented in contrMat2.")
      #C <- rbind(C)
      # TODO
    }
    return(C)    
  }
  stop("Unrecognized argument for 'type'.")
}

# TODO Merge getPairwiseContrasts and contrMat2.
# getPairwiseContrasts(model=2, v=5)
# getPairwiseContrasts(7, 3)
getPairwiseContrasts <- function(model, v) {
  pc <- parameterCount(model, v)
  contrasts <- list()
  p.prev <- 0
  p.follow <- sum(pc)
  for (p in pc) {
    Csub <- contrMat(n=rep(1, p), type="Tukey")
    class(Csub) <- "matrix"
    p.follow <- p.follow - p
    C <- as.matrix(cbind(matrix(0,dim(Csub)[1], p.prev), Csub,matrix(0,dim(Csub)[1], p.follow)))
    p.prev <- p.prev + p    
    contrasts <- c(contrasts, list(C))
  } 
  return(contrasts)
}


nrOfParameters <- function(model, v) {
  model <- getModelNr(model)
  if (model %in% c(3,9)) return(v)
  if (model %in% c(1, 4, 5, 6)) return(2*v)
  if (model %in% c(2, 8)) return(3*v)
  if (model==7) return(v+v+v*v) 
}


parameterCount <- function(model, v) {
  model <- getModelNr(model)
  if (model %in% c(2,8)) {
    return(c(v, v, v))
  } else if (model %in% c(3,9)) {
    return(c(v))
  } else if (model == 7) {
    return(c(v, v, v*v))
  } else if (model %in% c(1,4,5,6) ) {
    return(c(v, v))
  }  
}

corMat <- function(correlation, s, p, rho, q=0) {
  if (correlation=="equicorrelated") {
    V <- diag(p)
    for (i in 1:p) {
      for (j in 1:p) {
        if (i!=j) {
          V[i,j] <- rho
        }
      }
    }
    #V <- q*diag(p)
  } else if (correlation=="autoregressive") {
    V <- diag(p)
    for (i in 1:p) {
      for (j in 1:p) {
        V[i,j] <- rho^abs(i-j)
      }
    }
  }
  # Our design matrix is indexed p=1,1,1,2,2,2,3,3,3; s=1,2,3,1,2,3,1,2,3 therefore we have to exchange the arguments:
  if (q==0) {
    sigmaI <- kronecker(solve(V), diag(s)) #kronecker(diag(s), f(V))
  } else {
    sigmaI <- solve(q*diag(p*s)+(1-q)*kronecker(V, diag(s)))
  }
  return(sigmaI)
}
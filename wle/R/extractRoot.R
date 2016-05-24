#############################################################
#                                                           #
#	extractRoot.wle.glm function                              #
#	Author: Claudio Agostinelli                               #
#	E-mail: claudio@unive.it                                  #
#	Date: October, 09, 2012                                   #
#	Version: 0.1-2                                            #
#                                                           #
#	Copyright (C) 2012 Claudio Agostinelli                    #
#                                                           #
#############################################################

###############################################
extractRoot <- function(object, root=1, ...) UseMethod("extractRoot")

###############################################
extractRoot.wle.clm <- function(object, root=1, ...) {
  root <- round(as.numeric(root))
  if (length(root)!=1)
    stop("'root' arguments must be a vector of length 1")
  if (root > object$tot.sol | root <= 0)
    stop("'root' must be between 1 and object$tot.sol")

  object <- c(object, eval(parse(text=paste('object$root', root,"$model.store", sep=''))))
  #object[1:object$tot.sol] <- NULL
  class(object) <- c("wle.clm.root", "clm", "glm", "lm")
  return(object)
}

extractRoot.wle.glm <- function(object, root=1, ...) {
  root <- round(as.numeric(root))
  if (length(root)!=1)
    stop("'root' arguments must be a vector of length 1")
  if (root > object$tot.sol | root <= 0)
    stop("'root' must be between 1 and object$tot.sol")

  object <- c(object, eval(parse(text=paste('object$root', root, sep=''))))
  object[1:object$tot.sol] <- NULL
  class(object) <- c("wle.glm.root", "glm", "lm")            
  return(object)
}

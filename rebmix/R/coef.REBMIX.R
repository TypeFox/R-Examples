setMethod("coef",
          signature(object = "REBMIX"),
function(object, pos, ...)
{
  if (missing(object)) {
    stop(sQuote("object"), " object of class REBMIX is requested!", call. = FALSE)
  }

  if (!is.wholenumber(pos)) {
    stop(sQuote("pos"), " integer is requested!", call. = FALSE)
  }
  
  length(pos) <- 1

  if ((pos < 1) || (pos > nrow(object@summary))) {
    stop(sQuote("pos"), " must be greater than 0 and less or equal than ", nrow(object@summary), "!", call. = FALSE)
  }
  
  w <- matrix(object@w[[pos]], nrow = 1)
  
  c <- ncol(w)  
  
  rownames(w) <- "w"
  colnames(w) <- paste("comp", if (object@summary[pos, "c"] > 1) 1:object@summary[pos, "c"] else "", sep = "") 
  
  print(w, quote = FALSE, ...)
  
  d <- length(object@Variables)   
  
  names <- names(object@Theta[[pos]]) 

  Names <- names[grep("theta1", names, fixed = TRUE)]
  
  theta1 <- matrix(unlist(object@Theta[[pos]][Names]), ncol = d, byrow = TRUE)
  
  rownames(theta1) <- Names
  colnames(theta1) <- paste(1:d, sep = "")
  
  print(theta1, quote = FALSE, ...) 

  Names <- names[grep("theta2", names, fixed = TRUE)]
  
  theta2 <- matrix(unlist(object@Theta[[pos]][Names]), ncol = d, byrow = TRUE)
  
  rownames(theta2) <- Names
  colnames(theta2) <- paste(1:d, sep = "")
  
  print(theta2, quote = FALSE, ...)
  
  rm(list = ls())    
}) ## coef

setMethod("coef",
          signature(object = "REBMVNORM"),
function(object, pos, ...)
{
  if (missing(object)) {
    stop(sQuote("object"), " object of class REBMVNORM is requested!", call. = FALSE)
  }

  if (!is.wholenumber(pos)) {
    stop(sQuote("pos"), " integer is requested!", call. = FALSE)
  }
  
  length(pos) <- 1

  if ((pos < 1) || (pos > nrow(object@summary))) {
    stop(sQuote("pos"), " must be greater than 0 and less or equal than ", nrow(object@summary), "!", call. = FALSE)
  }
  
  w <- matrix(object@w[[pos]], nrow = 1)
  
  c <- ncol(w)  
  
  rownames(w) <- "w"
  colnames(w) <- paste("comp", if (object@summary[pos, "c"] > 1) 1:object@summary[pos, "c"] else "", sep = "") 
  
  print(w, quote = FALSE, ...)
  
  d <- length(object@Variables)   
  
  names <- names(object@Theta[[pos]]) 

  Names <- names[grep("theta1", names, fixed = TRUE)]
  
  theta1 <- matrix(unlist(object@Theta[[pos]][Names]), ncol = d, byrow = TRUE)
  
  rownames(theta1) <- Names
  colnames(theta1) <- paste(1:d, sep = "")
  
  print(theta1, quote = FALSE, ...) 

  Names <- names[grep("theta2", names, fixed = TRUE)]
  
  theta2 <- matrix(unlist(object@Theta[[pos]][Names]), ncol = d * d, byrow = TRUE)
  
  rownames(theta2) <- Names
  colnames(theta2) <- paste(rep(1:d, each = d), rep(1:d, d), sep = "-") 
  
  print(theta2, quote = FALSE, ...)
  
  rm(list = ls())   
}) ## coef

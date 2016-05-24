setMethod("summary", 
          signature(object = "REBMIX"),
function(object, ...)
{
  if (missing(object)) {
    stop(sQuote("object"), " object of class REBMIX is requested!", call. = FALSE)
  }
  
  p <- match(c("Dataset", "Preprocessing", "Criterion", "c", "v/k", "IC", "logL", "M"), names(object@summary), nomatch = 0)

  DF <- object@summary[p]

  is.num <- sapply(DF, is.number); DF[is.num] <- lapply(DF[is.num], as.number)

  print(DF, quote = FALSE, ...)

  cat(paste("Maximum logL = ", DF[object@pos, "logL"], " at pos = ", object@pos, ".\n", sep = "", collapse = ""))
  
  rm(list = ls()) 
}) ## summary

setMethod("summary", 
          signature(object = "REBMVNORM"),
function(object, ...)
{
  if (missing(object)) {
    stop(sQuote("object"), " object of class REBMVNORM is requested!", call. = FALSE)
  }
  
  p <- match(c("Dataset", "Preprocessing", "Criterion", "c", "v/k", "IC", "logL", "M"), names(object@summary), nomatch = 0)
  
  DF <- object@summary[p]

  is.num <- sapply(DF, is.number); DF[is.num] <- lapply(DF[is.num], as.number)

  print(DF, quote = FALSE, ...)  

  cat(paste("Maximum logL = ", DF[object@pos, "logL"], " at pos = ", object@pos, ".\n", sep = "", collapse = ""))
  
  rm(list = ls()) 
}) ## summary

setMethod("summary", 
          signature(object = "REBMIX.boot"),
function(object, ...)
{
  if (missing(object)) {
    stop(sQuote("object"), " object of class REBMIX.boot is requested!", call. = FALSE)
  }
  
  w.cv <- matrix(object@w.cv, nrow = 1) 
  
  c <- ncol(w.cv)   
  
  rownames(w.cv) <- "w.cv"
  colnames(w.cv) <- paste("comp", if (object@c.mode > 1) 1:object@c.mode else "", sep = "") 
  
  print(w.cv, quote = FALSE, ...)

  d <- length(object@x@Variables)  
  
  names <- names(object@Theta.cv) 

  Names <- names[grep("theta1", names, fixed = TRUE)]
  
  theta1.cv <- matrix(unlist(object@Theta.cv[Names]), ncol = d, byrow = TRUE)
  
  rownames(theta1.cv) <- Names
  colnames(theta1.cv) <- paste(1:d, sep = "")
  
  print(theta1.cv, quote = FALSE, ...) 
  
  Names <- names[grep("theta2", names, fixed = TRUE)]
  
  theta2.cv <- matrix(unlist(object@Theta.cv[Names]), ncol = d, byrow = TRUE)
  
  rownames(theta2.cv) <- Names
  colnames(theta2.cv) <- paste(1:d, sep = "") 
  
  print(theta2.cv, quote = FALSE, ...)

  cat(paste("Mode probability = ", as.number(object@c.prob), " at c = ", object@c.mode, " components.\n", sep = "", collapse = ""))
  
  rm(list = ls()) 
}) ## summary

setMethod("summary", 
          signature(object = "REBMVNORM.boot"),
function(object, ...)
{
  if (missing(object)) {
    stop(sQuote("object"), " object of class REBMVNORM.boot is requested!", call. = FALSE)
  }
  
  w.cv <- matrix(object@w.cv, nrow = 1) 
  
  c <- ncol(w.cv)   
  
  rownames(w.cv) <- "w.cv"
  colnames(w.cv) <- paste("comp", if (object@c.mode > 1) 1:object@c.mode else "", sep = "") 
  
  print(w.cv, quote = FALSE, ...)

  d <- length(object@x@Variables)  
  
  names <- names(object@Theta.cv) 

  Names <- names[grep("theta1", names, fixed = TRUE)]
  
  theta1.cv <- matrix(unlist(object@Theta.cv[Names]), ncol = d, byrow = TRUE)
  
  rownames(theta1.cv) <- Names
  colnames(theta1.cv) <- paste(1:d, sep = "")
  
  print(theta1.cv, quote = FALSE, ...) 
  
  Names <- names[grep("theta2", names, fixed = TRUE)]
  
  theta2.cv <- matrix(unlist(object@Theta.cv[Names]), ncol = d * d, byrow = TRUE)
  
  rownames(theta2.cv) <- Names
  colnames(theta2.cv) <- paste(rep(1:d, each = d), rep(1:d, d), sep = "-") 
  
  print(theta2.cv, quote = FALSE, ...)

  cat(paste("Mode probability = ", as.number(object@c.prob), " at c = ", object@c.mode, " components.\n", sep = "", collapse = ""))
  
  rm(list = ls()) 
}) ## summary

setMethod("summary", 
          signature(object = "RCLRMIX"),
function(object, ...)
{
  if (missing(object)) {
    stop(sQuote("object"), " object of class RCLRMIX is requested!", call. = FALSE)
  }
  
  print(object@Zp, quote = FALSE, ...)
   
  rm(list = ls()) 
}) ## summary

setMethod("summary", 
          signature(object = "RCLRMVNORM"),
function(object, ...)
{
  if (missing(object)) {
    stop(sQuote("object"), " object of class RCLRMVNORM is requested!", call. = FALSE)
  }
  
  print(object@Zp, quote = FALSE, ...)
  
  rm(list = ls()) 
}) ## summary

setMethod("summary", 
          signature(object = "RCLSMIX"),
function(object, ...)
{
  if (missing(object)) {
    stop(sQuote("object"), " object of class RCLSMIX is requested!", call. = FALSE)
  }
  
  CM <- as.data.frame(object@CM)
  
  colnames(CM) <- c("Test", "Predictive", "Frequency")
 
  print(CM, quote = FALSE, ...)
   
  cat(paste("Error = ", as.number(object@Error), ".\n", sep = "", collapse = ""))
  
  rm(list = ls()) 
}) ## summary

setMethod("summary", 
          signature(object = "RCLSMVNORM"),
function(object, ...)
{
  if (missing(object)) {
    stop(sQuote("object"), " object of class RCLSMVNORM is requested!", call. = FALSE)
  }
  
  CM <- as.data.frame(object@CM)
  
  colnames(CM) <- c("Test", "Predictive", "Frequency")
 
  print(CM, quote = FALSE, ...)
   
  cat(paste("Error = ", as.number(object@Error), ".\n", sep = "", collapse = ""))
  
  rm(list = ls()) 
}) ## summary

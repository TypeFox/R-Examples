setMethod("RCLRMIX",
          signature(model = "RCLRMIX"),
function(model, ...)
{
  Names <- names(model@x@Theta[[model@pos]])
    
  pdf <- unlist(model@x@Theta[[model@pos]][grep("pdf", Names)])
    
  theta1 <- unlist(model@x@Theta[[model@pos]][grep("theta1", Names)])
      
  theta1[is.na(theta1)] <- 0

  theta2 <- unlist(model@x@Theta[[model@pos]][grep("theta2", Names)])
      
  theta2[is.na(theta2)] <- 0

  c <- length(model@x@w[[model@pos]])

  w <- model@x@w[[model@pos]]
      
  d <- length(pdf) / c
  
  dataset <- as.matrix(model@x@Dataset[[model@pos]])
  
  n <- nrow(dataset)

  output <- .C("RCLRMIX",
    n = n,
    X = as.double(dataset),
    d = as.integer(d),
    c = as.integer(unlist(c)),
    w = as.double(unlist(w)),
    pdf = as.character(unlist(pdf)),
    theta1 = as.double(unlist(theta1)),
    theta2 = as.double(unlist(theta2)),
    Z = integer(n),
    error = integer(1),
    PACKAGE = "rebmix")

  if (output$error == 1) {
    stop("in RCLRMIX!", call. = FALSE); return(NA)
  }

  model@Zp <- as.factor(output$Z)

  rm(list = ls()[!(ls() %in% c("model"))])

  return(model)
}) ## RCLRMIX

setMethod("RCLRMIX",
          signature(model = "RCLRMVNORM"),
function(model, ...)
{
  Names <- names(model@x@Theta[[model@pos]])
    
  pdf <- unlist(model@x@Theta[[model@pos]][grep("pdf", Names)])
    
  theta1 <- unlist(model@x@Theta[[model@pos]][grep("theta1", Names)])
      
  theta1[is.na(theta1)] <- 0

  theta2 <- unlist(model@x@Theta[[model@pos]][grep("theta2", Names)])
      
  theta2[is.na(theta2)] <- 0

  c <- length(model@x@w[[model@pos]])

  w <- model@x@w[[model@pos]]
      
  d <- length(pdf) / c
  
  dataset <- as.matrix(model@x@Dataset[[model@pos]])
  
  n <- nrow(dataset)

  output <- .C("RCLRMVNORM",
    n = n,
    X = as.double(dataset),
    d = as.integer(d),
    c = as.integer(unlist(c)),
    w = as.double(unlist(w)),
    pdf = as.character(unlist(pdf)),
    theta1 = as.double(unlist(theta1)),
    theta2 = as.double(unlist(theta2)),
    Z = integer(n),
    error = integer(1),
    PACKAGE = "rebmix")

  if (output$error == 1) {
    stop("in RCLRMIX!", call. = FALSE); return(NA)
  }

  model@Zp <- as.factor(output$Z)

  rm(list = ls()[!(ls() %in% c("model"))])

  return(model)
}) ## RCLRMIX

setMethod("RCLRMIX",
          signature(model = "ANY"),
function(model,
  x,
  pos, 
  Zt, ...)
{
  digits <- getOption("digits"); options(digits = 15)
  
  message("RCLRMIX Version 2.8.1")
 
  flush.console()
  
  model <- new(model,
    x = x,
    pos = pos,
    Zt = Zt)
     
  model <- RCLRMIX(model = model, ...)
  
  options(digits = digits)  

  rm(list = ls()[!(ls() %in% c("model"))])

  return(model)
}) ## RCLRMIX

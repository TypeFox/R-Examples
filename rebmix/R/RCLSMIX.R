setMethod("RCLSMIX",
          signature(model = "RCLSMIX"),
function(model, ...)
{
  o <- length(model@x)

  d <- array(data = NA, dim = o, dimnames = NULL)

  c <- array(data = list(NULL), dim = c(o, model@s), dimnames = NULL)

  w <- array(data = list(NULL), dim = c(o, model@s), dimnames = NULL)

  pdf <- array(data = list(NULL), dim = c(o, model@s), dimnames = NULL)

  theta1 <- array(data = list(NULL), dim = c(o, model@s), dimnames = NULL)

  theta2 <- array(data = list(NULL), dim = c(o, model@s), dimnames = NULL)

  for (io in 1:o) {
    for (is in 1:model@s) {
      Names <- names(model@x[[io]]@Theta[[is]])
    
      pdf[[io, is]] <- unlist(model@x[[io]]@Theta[[is]][grep("pdf", Names)])
    
      theta1[[io, is]] <- unlist(model@x[[io]]@Theta[[is]][grep("theta1", Names)])
      
      theta1[[io, is]][is.na(theta1[[io, is]])] <- 0

      theta2[[io, is]] <- unlist(model@x[[io]]@Theta[[is]][grep("theta2", Names)])
      
      theta2[[io, is]][is.na(theta2[[io, is]])] <- 0

      c[[io, is]] <- length(model@x[[io]]@w[[is]])

      w[[io, is]] <- model@x[[io]]@w[[is]]
      
      d[io] <- length(pdf[[io, is]]) / c[[io, is]]
    }
  }

  output <- .C("RCLSMIX",
    n = model@ntest,
    X = as.double(unlist(model@Dataset)),
    s = as.integer(model@s),
    o = as.integer(o),
    d = as.integer(d),
    c = as.integer(unlist(c)),
    w = as.double(unlist(w)),
    pdf = as.character(unlist(pdf)),
    theta1 = as.double(unlist(theta1)),
    theta2 = as.double(unlist(theta2)),
    P = as.double(unlist(model@P)),
    Z = integer(model@ntest),
    error = integer(1),
    PACKAGE = "rebmix")

  if (output$error == 1) {
    stop("in RCLSMIX!", call. = FALSE); return(NA)
  }

  model@Zp <- as.factor(output$Z)

  rm(list = ls()[!(ls() %in% c("model"))])

  return(model)
}) ## RCLSMIX

setMethod("RCLSMIX",
          signature(model = "RCLSMVNORM"),
function(model, ...)
{
  o <- length(model@x)

  d <- array(data = NA, dim = o, dimnames = NULL)

  c <- array(data = list(NULL), dim = c(o, model@s), dimnames = NULL)

  w <- array(data = list(NULL), dim = c(o, model@s), dimnames = NULL)

  pdf <- array(data = list(NULL), dim = c(o, model@s), dimnames = NULL)

  theta1 <- array(data = list(NULL), dim = c(o, model@s), dimnames = NULL)

  theta2 <- array(data = list(NULL), dim = c(o, model@s), dimnames = NULL)

  for (io in 1:o) {
    for (is in 1:model@s) {
      Names <- names(model@x[[io]]@Theta[[is]])
    
      pdf[[io, is]] <- unlist(model@x[[io]]@Theta[[is]][grep("pdf", Names)])
    
      theta1[[io, is]] <- unlist(model@x[[io]]@Theta[[is]][grep("theta1", Names)])
      
      theta1[[io, is]][is.na(theta1[[io, is]])] <- 0

      theta2[[io, is]] <- unlist(model@x[[io]]@Theta[[is]][grep("theta2", Names)])
      
      theta2[[io, is]][is.na(theta2[[io, is]])] <- 0

      c[[io, is]] <- length(model@x[[io]]@w[[is]])

      w[[io, is]] <- model@x[[io]]@w[[is]]
      
      d[io] <- length(pdf[[io, is]]) / c[[io, is]]
    }
  }
  
  output <- .C("RCLSMVNORM",
    n = model@ntest,
    X = as.double(unlist(model@Dataset)),
    s = as.integer(model@s),
    o = as.integer(o),
    d = as.integer(d),
    c = as.integer(unlist(c)),
    w = as.double(unlist(w)),
    pdf = as.character(unlist(pdf)),
    theta1 = as.double(unlist(theta1)),
    theta2 = as.double(unlist(theta2)),
    P = as.double(unlist(model@P)),
    Z = integer(model@ntest),
    error = integer(1),
    PACKAGE = "rebmix")

  if (output$error == 1) {
    stop("in RCLSMIX!", call. = FALSE); return(NA)
  }

  model@Zp <- as.factor(output$Z)

  rm(list = ls()[!(ls() %in% c("model"))])

  return(model)
}) ## RCLSMIX

setMethod("RCLSMIX",
          signature(model = "ANY"),
function(model,
  x,
  Dataset, 
  Zt, ...)
{
  digits <- getOption("digits"); options(digits = 15)
  
  message("RCLSMIX Version 2.8.1")
 
  flush.console()
  
  model <- new(model,
    x = x,
    Dataset = Dataset,
    Zt = Zt)
     
  model <- RCLSMIX(model = model, ...)
  
  model@CM <- table(model@Zt, model@Zp)
  
  model@Accuracy <- sum(diag(model@CM)) / model@ntest
  
  model@Error <- 1.0 - model@Accuracy
  
  model@Precission <- diag(model@CM) / apply(model@CM, 1, sum)
  
  model@Sensitivity <- diag(model@CM) / apply(model@CM, 2, sum)
  
  model@Specificity <- (model@ntest - apply(model@CM, 1, sum)) / (model@ntest - apply(model@CM, 2, sum))
  
  options(digits = digits)  

  rm(list = ls()[!(ls() %in% c("model"))])

  return(model)
}) ## RCLSMIX

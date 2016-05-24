setMethod("boot", 
          signature(x = "REBMIX"),
function(x,
  pos,
  Bootstrap,
  B, 
  n,
  replace, 
  prob, ...)
{
  digits <- getOption("digits"); options(digits = 15)
  
  model <- new("REBMIX.boot", 
    x = x, 
    pos = pos,
    Bootstrap = Bootstrap,
    B = B,
    n = n,
    replace = replace,
    prob = prob)
    
  if (model@Bootstrap == .rebmix.boot$Bootstrap[1]) {
    bsample <- RNGMIX(Dataset.name = paste("bsample_", 1:model@B, sep = ""),
      n = ceiling(model@n * as.numeric(model@x@w[[model@pos]])),
      Theta = model@x@Theta[[model@pos]], ...)
      
    Dataset <- bsample@Dataset
  }
  else
  if (model@Bootstrap == .rebmix.boot$Bootstrap[2]) {
    dataset <- as.matrix(model@x@Dataset[[model@pos]])
    
    Dataset <- list()
    
    for (i in 1:model@B) {
      R1 <- sample.int(n = nrow(dataset), size = model@n, replace = replace, prob = if (length(prob) == 0) NULL else prob)
      
      Dataset[[i]] <- as.data.frame(dataset[R1,], stringsAsFactors = FALSE)
    }
  }
  
  if (length(model@x@theta1) > 0) {
    theta1 <- model@x@theta1; theta1[is.na(theta1)] <- 0
  }
  else {
    theta1 <- numeric()
  }    
    
  if (length(model@x@theta2) > 0) {
    theta2 <- model@x@theta2; theta2[is.na(theta2)] <- 0
  }
  else {
    theta2 <- numeric()
  }  
    
  d <- length(model@x@pdf)
  
  bsampleest <- REBMIX(Dataset = Dataset,
    Preprocessing = as.character(model@x@summary[model@pos, "Preprocessing"]),
    cmax = as.numeric(model@x@summary[model@pos, "cmax"]),
    Criterion = as.character(model@x@summary[model@pos, "Criterion"]),
    Variables = model@x@Variables,
    pdf = model@x@pdf,
    theta1 = theta1,
    theta2 = theta2,
    K = eval(parse(text = as.character(model@x@summary[model@pos, "K"]))),
    y0 = model@x@y0,    
    ymin = model@x@ymin,
    ymax = model@x@ymax,
    ar = as.numeric(model@x@summary[model@pos, "ar"]),
    Restraints = as.character(model@x@summary[model@pos, "Restraints"]))

  freq <- table(as.numeric(bsampleest@summary$c))
  
  c <- as.integer(names(freq)[which.max(freq)])
  
  model@c <- as.numeric(bsampleest@summary$c)
  model@c.se <- sd(model@c)
  model@c.cv <- model@c.se / mean(model@c)  
  
  w <- bsampleest@w[model@c == c]
  
  Theta <- bsampleest@Theta[model@c == c]
  
  model@c.mode <- c
  model@c.prob <- length(w) / model@B
  
  model@w <- matrix(unlist(w), ncol = c, byrow = TRUE)
  
  colnames(model@w) <- paste("comp", if (c > 1) 1:c else "", sep = "")
  rownames(model@w) <- paste(which(bsampleest@summary$c == c), sep = "")
  
  model@w.se <- apply(model@w, 2, sd)
  model@w.cv <- model@w.se / apply(model@w, 2, mean)

  for (i in 1:model@c.mode) {
    theta1 <- paste("theta1.",  i, sep = "")

    model@Theta[[theta1]] <- NULL

    for (j in 1:length(Theta)) {
      model@Theta[[theta1]] <- c(model@Theta[[theta1]], Theta[[j]][[theta1]])
    }

    model@Theta[[theta1]] <- matrix(model@Theta[[theta1]], ncol = d, byrow = TRUE)

    colnames(model@Theta[[theta1]]) <- paste(1:d, sep = "")
    rownames(model@Theta[[theta1]]) <- paste(which(bsampleest@summary$c == c), sep = "")

    theta1.se <- paste("theta1.",  i, ".se", sep = "")
    theta1.cv <- paste("theta1.",  i, ".cv", sep = "")

    model@Theta.se[[theta1.se]] <- apply(model@Theta[[theta1]], 2, sd)
    model@Theta.cv[[theta1.cv]] <- model@Theta.se[[theta1.se]] / apply(model@Theta[[theta1]], 2, mean)
  }
  
  for (i in 1:model@c.mode) {
    theta2 <- paste("theta2.",  i, sep = "")

    model@Theta[[theta2]] <- NULL

    for (j in 1:length(Theta)) {
      model@Theta[[theta2]] <- c(model@Theta[[theta2]], Theta[[j]][[theta2]])
    }

    model@Theta[[theta2]] <- matrix(model@Theta[[theta2]], ncol = d, byrow = TRUE)

    colnames(model@Theta[[theta2]]) <- paste(1:d, sep = "")
    rownames(model@Theta[[theta2]]) <- paste(which(bsampleest@summary$c == c), sep = "")

    theta2.se <- paste("theta2.",  i, ".se", sep = "")
    theta2.cv <- paste("theta2.",  i, ".cv", sep = "")

    model@Theta.se[[theta2.se]] <- apply(model@Theta[[theta2]], 2, sd)
    model@Theta.cv[[theta2.cv]] <- model@Theta.se[[theta2.se]] / apply(model@Theta[[theta2]], 2, mean)
  }
  
  options(digits = digits)  

  rm(list = ls()[!(ls() %in% c("model"))]) 
   
  return(model)  
}) ## boot

setMethod("boot", 
          signature(x = "REBMVNORM"),
function(x,
  pos,
  Bootstrap,
  B, 
  n,
  replace, 
  prob, ...)
{
  digits <- getOption("digits"); options(digits = 15)
  
  model <- new("REBMVNORM.boot", 
    x = x, 
    pos = pos,
    Bootstrap = Bootstrap,
    B = B,
    n = n,
    replace = replace,
    prob = prob)
    
  if (model@Bootstrap == .rebmix.boot$Bootstrap[1]) {
    bsample <- RNGMIX(model = "RNGMVNORM", Dataset.name = paste("bsample_", 1:model@B, sep = ""),
      n = ceiling(model@n * as.numeric(model@x@w[[model@pos]])),
      Theta = model@x@Theta[[model@pos]], ...)
      
    Dataset <- bsample@Dataset
  }
  else
  if (model@Bootstrap == .rebmix.boot$Bootstrap[2]) {
    dataset <- as.matrix(model@x@Dataset[[model@pos]])
    
    Dataset <- list()
    
    for (i in 1:model@B) {
      R1 <- sample.int(n = nrow(dataset), size = model@n, replace = replace, prob = if (length(prob) == 0) NULL else prob)
      
      Dataset[[i]] <- as.data.frame(dataset[R1,], stringsAsFactors = FALSE)
    }
  }
  
  if (length(model@x@theta1) > 0) {
    theta1 <- model@x@theta1; theta1[is.na(theta1)] <- 0
  }
  else {
    theta1 <- numeric()
  }    
    
  if (length(model@x@theta2) > 0) {
    theta2 <- model@x@theta2; theta2[is.na(theta2)] <- 0
  }
  else {
    theta2 <- numeric()
  }  
    
  d <- length(model@x@pdf)
  
  bsampleest <- REBMIX(model = "REBMVNORM",
    Dataset = Dataset,
    Preprocessing = as.character(model@x@summary[model@pos, "Preprocessing"]),
    cmax = as.numeric(model@x@summary[model@pos, "cmax"]),
    Criterion = as.character(model@x@summary[model@pos, "Criterion"]),
    Variables = model@x@Variables,
    pdf = model@x@pdf,
    theta1 = theta1,
    theta2 = theta2,
    K = eval(parse(text = as.character(model@x@summary[model@pos, "K"]))),
    y0 = model@x@y0,    
    ymin = model@x@ymin,
    ymax = model@x@ymax,
    ar = as.numeric(model@x@summary[model@pos, "ar"]),
    Restraints = as.character(model@x@summary[model@pos, "Restraints"]))

  freq <- table(as.numeric(bsampleest@summary$c))
  
  c <- as.integer(names(freq)[which.max(freq)])
  
  model@c <- as.numeric(bsampleest@summary$c)
  model@c.se <- sd(model@c)
  model@c.cv <- model@c.se / mean(model@c)  
  
  w <- bsampleest@w[model@c == c]
  
  Theta <- bsampleest@Theta[model@c == c]
  
  model@c.mode <- c
  model@c.prob <- length(w) / model@B
  
  model@w <- matrix(unlist(w), ncol = c, byrow = TRUE)
  
  colnames(model@w) <- paste("comp", if (c > 1) 1:c else "", sep = "")
  rownames(model@w) <- paste(which(bsampleest@summary$c == c), sep = "")
  
  model@w.se <- apply(model@w, 2, sd)
  model@w.cv <- model@w.se / apply(model@w, 2, mean)

  for (i in 1:model@c.mode) {
    theta1 <- paste("theta1.",  i, sep = "")

    model@Theta[[theta1]] <- NULL

    for (j in 1:length(Theta)) {
      model@Theta[[theta1]] <- c(model@Theta[[theta1]], Theta[[j]][[theta1]])
    }

    model@Theta[[theta1]] <- matrix(model@Theta[[theta1]], ncol = d, byrow = TRUE)

    colnames(model@Theta[[theta1]]) <- paste(1:d, sep = "")
    rownames(model@Theta[[theta1]]) <- paste(which(bsampleest@summary$c == c), sep = "")

    theta1.se <- paste("theta1.",  i, ".se", sep = "")
    theta1.cv <- paste("theta1.",  i, ".cv", sep = "")

    model@Theta.se[[theta1.se]] <- apply(model@Theta[[theta1]], 2, sd)
    model@Theta.cv[[theta1.cv]] <- model@Theta.se[[theta1.se]] / apply(model@Theta[[theta1]], 2, mean)
  }
  
  for (i in 1:model@c.mode) {
    theta2 <- paste("theta2.",  i, sep = "")

    model@Theta[[theta2]] <- NULL

    for (j in 1:length(Theta)) {
      model@Theta[[theta2]] <- c(model@Theta[[theta2]], Theta[[j]][[theta2]])
    }

    model@Theta[[theta2]] <- matrix(model@Theta[[theta2]], ncol = d * d, byrow = TRUE)
    
    colnames(model@Theta[[theta2]]) <- paste(rep(1:d, each = d), rep(1:d, d), sep = "-")
    rownames(model@Theta[[theta2]]) <- paste(which(bsampleest@summary$c == c), sep = "")

    theta2.se <- paste("theta2.",  i, ".se", sep = "")
    theta2.cv <- paste("theta2.",  i, ".cv", sep = "")

    model@Theta.se[[theta2.se]] <- apply(model@Theta[[theta2]], 2, sd)
    model@Theta.cv[[theta2.cv]] <- model@Theta.se[[theta2.se]] / apply(model@Theta[[theta2]], 2, mean)
  }
  
  options(digits = digits)  

  rm(list = ls()[!(ls() %in% c("model"))]) 
   
  return(model)  
}) ## boot

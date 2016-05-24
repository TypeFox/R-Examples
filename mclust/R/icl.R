##
## Integrated Complete-data Likelihood (ICL) Criterion
##

icl <- function(object, ...) UseMethod("icl")

icl.Mclust <- function(object, ...)
{
  n <- object$n
  # G <- object$G + ifelse(is.na(object$hypvol),0,1)
  z <- object$z
  if(is.null(z)) z <- matrix(1, nrow = n, ncol = 1)
  C <- matrix(0, n, ncol(z))
  for(i in 1:n) 
    C[i, which.max(z[i,])] <- 1
  object$bic + 2*sum(C * ifelse(z > 0, log(z), 0))
}

icl.MclustDA <- function(object, ...)
{
  n <- object$n
  z <- predict(object)$z
  df <- object$df
  if(is.null(z)) z <- matrix(1, nrow = n, ncol = 1)
  C <- matrix(0, n, ncol(z))
  for(i in 1:n) 
    C[i, which.max(z[i,])] <- 1
  object$bic + 2*sum(C * ifelse(z > 0, log(z), 0))
}

mclustICL <- function(data, G = NULL, modelNames = NULL, 
                      initialization = list(hcPairs=NULL, subset=NULL, noise=NULL),  
                      x = NULL, ...)
{
  call <- match.call()
  data <- data.matrix(data)
  n <- nrow(data)
  d <- ncol(data)
  mc <- match.call(expand.dots = TRUE)
  mc[[1]] <- as.name("mclustBIC")
  mc[[2]] <- data
  BIC <- eval(mc, parent.frame())
  # browser()
  class(BIC) <- "mclustBIC"
  G <- attr(BIC, "G")
  modelNames <- attr(BIC, "modelNames")
  ICL <- matrix(NA, nrow = length(G), ncol = length(modelNames))
  mostattributes(ICL) <- attributes(BIC)
  if(!is.null(x))
    { 
      r <- match(as.character(G), rownames(x), nomatch = 0)
      c <- match(modelNames, colnames(x), nomatch = 0)
      ICL[r,c] <- BIC[r,c]
  }
  for(i in 1:nrow(ICL))
     { for(j in 1:ncol(ICL))
          { if(is.na(BIC[i,j])) next() # not fitted
            if(!is.na(ICL[i,j])) next() # already available
            Sumry <- summary(BIC, data, G = G[i], modelNames = modelNames[j])
            ICL[i,j] <- icl.Mclust(Sumry)
       }
  }
  class(ICL) <- "mclustICL" # "mclustBIC"
  attr(ICL, "criterion") <- "ICL"
  return(ICL)
}

print.mclustICL <- function (x, pick = 3, ...) 
{
  subset <- !is.null(attr(x, "subset"))
  oldClass(x) <- attr(x, "args") <- NULL
  attr(x, "criterion") <- NULL
  attr(x, "control") <- attr(x, "initialization") <- NULL
  attr(x, "oneD") <- attr(x, "warn") <- attr(x, "Vinv") <- NULL
  attr(x, "prior") <- attr(x, "G") <- attr(x, "modelNames") <- NULL
  ret <- attr(x, "returnCodes") == -3
  n <- attr(x, "n")
  d <- attr(x, "d")
  attr(x, "returnCodes") <- attr(x, "n") <- attr(x, "d") <- NULL
  
  oldClass(x) <- attr(x, "args") <- attr(x, "criterion") <- NULL 
  cat("Integrated Complete-data Likelihood (ICL) criterion:\n")
  print(x)
  cat("\n")
  cat("Top", pick, "models based on the ICL criterion:\n")
  print(pickBIC(x, pick), ...)
  invisible()
}

summary.mclustICL <- function(object, G, modelNames, ...)
{
  if(!missing(G)) 
    object <- object[rownames(object) %in% G,,drop=FALSE]
  if(!missing(modelNames)) 
    object <- object[,colnames(object) %in% modelNames,drop=FALSE]
  structure(pickBIC(object, ...),
            class = "summary.mclustICL")
}

print.summary.mclustICL <- function(x, digits = getOption("digits"), ...)
{
  cat("Best ICL values:\n")
  x <- drop(as.matrix(x))
  x <- rbind(ICL = x, "ICL diff" = x - max(x))
  print(x, digits = digits)
  invisible()
}


plot.mclustICL <- function(x, ylab = "ICL", ...) 
{
  plot.mclustBIC(x, ylab = ylab, ...)  
}


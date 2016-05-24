# Remove r() from fromula (possibly convert to lmer formula)
rparse <- function (f, REML = FALSE) {
  if (class(f) != "formula") 
    stop("'f' must be a formula")
  right <- attr(terms(f),"term.labels") # Let R split short-hand notations first
  if( length(right) == 0){
    return(f)
  }
  n <- length(right)
  result      <- character(n)
  result.REML <- character(n)
  
  # Main recursive loop extracting effects without r()
  for(i in 1:n){
    parsecall <- function(x) {
      if (class(x) != "call") 
        stop("'x' must be a call")
      if (length(x[[1]]) == 1) {
        return(x)
      }
      if (length(x[[1]]) == 2 && x[[1]][[1]] == "r") {
        return(x[[1]][2])
      }
      for (i in 2:length(x[[1]])) x[[1]][i] <- parsecall(x[[1]][i])
      return(x)
    }
    result[[i]] <- as.character(parsecall(formula(paste("~",paste(right[i],sep="+")))[2]))
  }
  f[3] <- formula(paste("~", paste(result, sep="", collapse="+")))[2]
  
  # Recursive loop adding (1 | x) notation for REML estimation
  if(REML){
    for(i in 1:n){
      parsecall <- function(x) {
        if (class(x) != "call") 
          stop("'x' must be a call")
        if (length(x[[1]]) == 1) {
          return(FALSE)
        }
        if (length(x[[1]]) == 2 && x[[1]][[1]] == "r") {
          return(TRUE)
        } else {
          tmp <- logical(length(x[[1]]))
          for (j in 2:length(x[[1]])) tmp[j-1] <- parsecall(x[[1]][j])
          return(any(tmp))
        }
      }
      ran <- parsecall(formula(paste("~",paste(right[i],sep="+")))[2])
      result.REML[i] <- ifelse(ran,as.character(formula(paste("~(1 | ",result[[i]],")",sep=""))[2]), result[[i]])
    }
    f[3] <- formula(paste("~", paste(result.REML, sep="", collapse="+")))[2]
  }
  f
}

# Extract variance components from lmer model
# VarComps <- function(model){
# varcorr <- lme4::VarCorr(model)
# ran.names <- names(varcorr)
# lr <- length(ran.names)
# vc <- numeric(lr+1)
# for(i in 1:lr){
# vc[i] <- varcorr[[i]][1]
# }
# vc[lr+1] <- attr( varcorr, "sc")^2
# names(vc) <- c(ran.names,"residuals")
# vc
# }


# Utility function for checking if a model has balanced data
is.balanced <- function(object){
  if(missing(object) && "package:Rcmdr"%in%search()){
    eval(parse(text=paste("effects <- as.character(attr(terms(formula(", ActiveModel(),")),'variables')[-c(1,2)])", sep="")))
    eval(parse(text=paste("balanced <- length(unique(xtabs(~", paste(effects,sep="",collapse="+"), ", data=", ActiveDataSet(), ")))==1", sep="")))
  } else {
    effects <- as.character(attr(terms(formula(object)),'variables')[-c(1,2)])
    eval(parse(text=paste("balanced <- length(unique(xtabs(~", paste(effects,sep="",collapse="+"), ", data=object$model)))==1",sep="")))
  }
  balanced
}


# Extract variable names from formula
fparse <- function(f) {
  if (class(f) != "formula") stop("'f' must be a formula")
  right <- f[3]
  parsecall <- function(x) {
    if (class(x) != "call") stop("'x' must be a call")
    if (length(x[[1]]) == 1) {
      if(is.numeric(x[[1]])) {
        return()
      } else {
        return(as.character(x[[1]]))
      }
    }
    res <- list()
    for (i in 2:length(x[[1]]))
      res[[i-1]] <- parsecall(x[[1]][i])
    return(unlist(res))
  }
  unique(parsecall(right))
}

# Studentized range for 1 df
qtukey1df <- matrix(c(8.929,13.437,16.358,18.488,20.15,21.504,22.642,23.621,24.477,25.237,25.918,26.536,27.1,27.618,28.097,28.542,28.958,29.347,29.713,
                      17.969,26.976,32.819,37.082,40.408,43.119,45.397,47.357,49.071,50.592,51.957,53.194,54.323,55.361,56.32,57.212,58.044,58.824,59.558,
                      35.99,54,65.69,74.22,80.87,86.29,90.85,94.77,98.2,101.3,104,106.5,108.8,110.8,112.7,114.45,116.2,117.7,119.2,
                      90.024,135.041,164.258,185.575,202.21,215.769,227.166,236.966,245.542,253.151,259.979,266.165,271.812,277.003,281.803,286.263,290.426,294.328,297.997), 19,4)
dimnames(qtukey1df) <- list(k = 2:20, P = c(0.9, 0.95, 0.975, 0.99))

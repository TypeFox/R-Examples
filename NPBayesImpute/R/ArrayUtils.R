###############################################################################
#Functions for manipulating matrices and dataframes with categorical values
# (c) Daniel Manrique-Vallier 2012
###############################################################################

.fn_df_nlevels <- function(d){
    sapply(names(d), FUN = function(x)nlevels(d[,x]))
}

.fn_df_discretize <- function(d, cols = 1:NCOL(d)){
  for (i in cols){
    d[,i] <- factor(d[,i])
  }
  return(d)
}
.fn_df_levels <- function(d, cols = 1:NCOL(d)){

  res <- list()
  for (i in cols){
    res[[names(d)[i]]] <- levels(d[1,i])
  }
  return(res)
}

.fn_numeric_to_factor <- function(
  x, list_levls = rep(list(1:2), NCOL(x)), cols = 1:NCOL(x), missing = -1, offset = 0){
  x <- as.data.frame(x)
  for (c in cols){
    lv <- list_levls[[c]]
    nlv <- length(lv)
    x[,c] <- factor(x[,c], levels = 1:nlv + offset)
    levels(x[,c]) <- lv
  }
  return(x)
}

.fn_apply_levels_from <- function(dest, from, cols = 1:NCOL(from)){
  dest <- as.data.frame(dest)
  lev <- sapply(cols, function(x)levels(from[1,x]), simplify=FALSE)
  for (i in cols){
    dest[,i] <- factor(dest[,i], levels = lev[[i]])
  }
  names(dest[,cols]) <- names(from[,cols])
  return(dest)
}

.fn_dataframe2num_matrix <- function(d, offset = -1, missing = -1, C_style= FALSE){
	if (is.null(d)) {
		d<- matrix(1,1)
		names(d)[1] <- "Dummy"
	}
  	r <- data.matrix(d) + offset
  r[is.na(r)] <- missing
  if(C_style){
    r <- t(r)
  }
  class(r) <- c('matrix', 'DMVmatrix')
  attr(r, 'levels') <- .fn_df_levels(d)
  attr(r, 'na.code') <- missing
  attr(r, 'offset') <- offset
  attr(r, 'C_style') <- C_style
  return(r)
}

.dataframe2matrix <- function(d){
	return(.fn_dataframe2num_matrix(d, offset = -1, missing=-1, 	C_style=T))
}

.checkx <- function(X) {
  # the same check as in CreateModel. 
  if (!is.data.frame(X)) {
    stop(simpleError("Input data must be a data.frame object"));
  } else {
    col_names <- names(X)
    is_factor <- sapply(X[,col_names], is.factor)
    if (any(!is_factor)) {
      stop(simpleError("All columns in input dataframe must be factors"));
    }
    if (dim(X)[2] > dim(X)[1]) {
      stop(simpleError("Please make sure each be an observation in input data"));
    }
  }
}
UpdateX <- function(model,X) {
  .checkx(X)
  x <- .dataframe2matrix(X)
  model$UpdateX(x)
}

CreateModel <- function(X,MCZ,K, Nmax, aalpha, balpha) {
  .checkx(X)
  if (is.null(MCZ)) {
    Nmax <- 0
  } else {
    if (!is.data.frame(MCZ)) {
      stop(simpleError("Input MCZ must be a data.frame object"));
    } else {
      col_names <- names(MCZ)
      is_factor <- sapply(MCZ[,col_names], is.factor)
      if (any(!is_factor)) {
        stop(simpleError("All columns in input MCZ dataframe must be factors"));
      }
    } 
    if (dim(MCZ)[2] != dim(X)[2]) {
      stop(simpleError("Input Data and MCZ dimensions do not match"));
    }
  }
  
  if (K <=1) {
    stop(simpleError("Number of components must be at least 2"));
  }
  
  x <- .dataframe2matrix(X)
  mcz <- .dataframe2matrix(MCZ)
  model <- new(Lcm,x,mcz,K, Nmax, aalpha, balpha)
  model$SetXAsDataframe(X)
  return(model)
}
GetDataFrame <- function(dest, from, cols = 1:NCOL(from)){
  t <- dest + 1
  y <- as.data.frame(t)
  lev <- sapply(cols, function(x)levels(from[1,x]), simplify=FALSE)
  for (i in cols){
    li <- lev[[i]] 
    y[,i] <- factor(li[t[,i]],levels = lev[[i]])
  }
  names(y) <- names(from)
  rownames(y) <- rownames(from)  
  return(y)
}

GetMCZ <- function(dest, from, mcz, cols = 1:NCOL(from)){
  t <- dest + 1
  y <- as.data.frame(t)
  lev <- sapply(cols, function(x)levels(from[1,x]), simplify=FALSE)
  for (i in cols){
    li <- lev[[i]] 
    temp <- t[,i]
    temp[temp ==0] <- NA
    y[,i] <- li[temp]
  }
  names(y) <- names(mcz)
  return(y)
}



## =============================================================================
## =============================================================================
## Matrix plotting - most of these function were slightly modified from 
## similar functions in the deSolve and rootSolve packages.
## =============================================================================
## =============================================================================


## =============================================================================
## Describes a (list of) matrices
## =============================================================================

Mdescribe <- function (M, ..., 
                       select = NULL, which = select, 
                       subset = NULL) {

  getnames <- function(x) 
    if (is.null (cn <- colnames(x))) return (1:ncol(x)) else return(cn)
  
  #dirty trick to get ALL names of ellipsis 
  NN <- deparse(substitute(x(...)))
  NN <- gsub("x(","",NN,fixed=TRUE)
  NN <- gsub(")","",NN)
  NN <- gsub(" ","",NN)    
  dotnames <- unlist(strsplit(NN, ","))
  xnames <- c(deparse(substitute(M)), dotnames)                

  # The ellipsis
  ldots   <- list(...)
  
  Dots    <- splitobject(ldots, M, xnames)
  x2      <- Dots$x2
  nother  <- Dots$nother
  nx      <- nother + 1 # total number of objects to be plotted
  varnames <- getnames(x2[[1]])

 # variables to be described
  Which <- which
  if (is.null(Which)) {
    for (i in 1: length(x2))
      Which <- c(Which,getnames(x2[[i]]))
    Which <- unique(Which)
  }

  np      <- length(Which)  
 # Position of variables to be described in "M" and other matrices
  xWhich <- list()

  for (i in 1: length(x2))
    xWhich[[i]] <- selectvar(Which, getnames(x2[[i]]))

  if (! is.character(Which)) 
    Which <- varnames[xWhich[[1]]]

  if (!missing(subset)){
    isub <- list()
    for (i in 1:nx) {
      e <- substitute(subset)
      r <- eval(e, as.data.frame(x2[[i]]), parent.frame())
      if (!is.logical(r))
        stop("'subset' must evaluate to logical")
      isub[[i]] <- r & !is.na(r)
    }  
  } else isub <- rep(TRUE, nx)

 desc <- data.frame()
 # LOOP for each selected variable 

  for (ip in 1 : np) {
   for (j in 1:nx) {
      ix <- xWhich[[j]][[ip]]      # position of variable in 'x2'
      if (!is.na(ix)) {
        y <- x2[[j]][isub[[j]], ix]
        if (is.factor(y) | is.character(y)) {
        desc <- rbind(desc, data.frame(variable = getnames(x2[[j]])[ix], 
          object = names(x2)[j], factor_or_char = TRUE, 
          n = length(y), missing = sum(is.na(y)), 
          unique = length(unique(y)), Mean = NA,
          Sd = NA, Min = NA,  p0.05 = NA, p0.1 = NA, 
          p0.5 = NA, p0.9 = NA, p0.95 = NA, Max = NA))
        } else {
          y <- as.numeric(y)
          Quant <- as.vector(quantile(y, prob = c(0.05, 0.1, 0.5, 0.9, 0.95), na.rm = TRUE))
          ndesc <- data.frame(variable = getnames(x2[[j]])[ix], 
            object = names(x2)[j], factor_or_char = FALSE, 
            n = length(y), missing = sum(is.na(y)), 
            unique = length(unique(y)), Mean = mean(y, na.rm = TRUE),
            Sd = sd(y, na.rm = TRUE), Min = min(y, na.rm = TRUE), 
            p0.05 = Quant[1], p0.1 = Quant[2], p0.5 = Quant[3], 
            p0.9 = Quant[4], p0.95 = Quant[5], Max = max(y, na.rm = TRUE))
          desc <- rbind(desc, ndesc)
       }
    }
  }
 }
 if (nx  == 1) 
 desc <- desc [, -2]
 desc    
}

Msummary <- function (M, ..., 
                       select = NULL, which = select, 
                       subset = NULL) {

  getnames <- function(x) 
    if (is.null (cn <- colnames(x))) return (1:ncol(x)) else return(cn)
  
  #dirty trick to get ALL names of ellipsis 
  NN <- deparse(substitute(x(...)))
  NN <- gsub("x(","",NN,fixed=TRUE)
  NN <- gsub(")","",NN)
  NN <- gsub(" ","",NN)    
  dotnames <- unlist(strsplit(NN, ","))
  xnames <- c(deparse(substitute(M)), dotnames)                

  # The ellipsis
  ldots   <- list(...)
  
  Dots    <- splitobject(ldots, M, xnames)
  x2      <- Dots$x2
  nother  <- Dots$nother
  nx      <- nother + 1 # total number of objects to be plotted
  varnames <- getnames(x2[[1]])

 # variables to be described
  Which <- which
  if (is.null(Which)) {
    for (i in 1: length(x2))
      Which <- c(Which,getnames(x2[[i]]))
    Which <- unique(Which)
  }

  np      <- length(Which)  
 # Position of variables to be described in "M" and other matrices
  xWhich <- list()

  for (i in 1: length(x2))
    xWhich[[i]] <- selectvar(Which, getnames(x2[[i]]))

  if (! is.character(Which)) 
    Which <- varnames[xWhich[[1]]]

  if (!missing(subset)){
    isub <- list()
    for (i in 1:nx) {
      e <- substitute(subset)
      r <- eval(e, as.data.frame(x2[[i]]), parent.frame())
      if (!is.logical(r))
        stop("'subset' must evaluate to logical")
      isub[[i]] <- r & !is.na(r)
    }  
  } else isub <- rep(TRUE, nx)

 desc <- data.frame()
 # LOOP for each selected variable 

  for (ip in 1 : np) {
   for (j in 1:nx) {
      ix <- xWhich[[j]][[ip]]      # position of variable in 'x2'
      if (!is.na(ix)) {
        y <- x2[[j]][isub[[j]], ix]
        if (is.factor(y) | is.character(y)) {
        desc <- rbind(desc, data.frame(variable = getnames(x2[[j]])[ix], 
          object = names(x2)[j], factor_or_char = TRUE, 
          Min. = NA,  "1st Qu." = NA, Median = NA, Mean = NA, "3rd Qu." = NA, Max. = NA))
        } else {
          y <- as.numeric(y)
          Quant <- as.vector(quantile(y, na.rm = TRUE))
          ndesc <- data.frame(variable = getnames(x2[[j]])[ix], 
            object = names(x2)[j], factor_or_char = FALSE, 
             Min. = Quant[1],  "1st Qu." = Quant[2], Median = Quant[3], 
             Mean = mean(y, na.rm = TRUE), "3rd Qu." = Quant[4], Max. = Quant[5])
          desc <- rbind(desc, ndesc)
       }
    }
  }
 }
 if (nx  == 1) 
 desc <- desc [, -2]
 desc    
}

splitobject <- function(ldots, x, xnames){
  x2      <- list()
  nother <- 0
  islist <- (! is.data.frame(x) & is.list(x))
  
  if (! islist) {
    x2[[1]] <- x
    names(x2)[1] <- xnames[1]
  } else {
    x2 <- x
    nother <- length(x) - 1
  }

  dots   <- list()
  nd     <- 0
  ndots <- xnames[-1]
    
  if (length(ldots) > 0)
    for ( i in 1:length(ldots))
      if ("matrix" %in% class(ldots[[i]]) | "data.frame" %in% class(ldots[[i]])) { 
        nother <- nother + 1        
        x2[[nother + 1]] <- ldots[[i]]
        if (is.null(ndots[i]))
          names(x2)[nother+1] <- nother +1
        else 
          names(x2)[nother+1] <- ndots[i]
        # a list of matrix objects
      } else if (is.list(ldots[[i]]) & 
        ("matrix" %in% class(ldots[[i]][[1]]) | 
         "data.frame" %in% class(ldots[[i]][[1]]))) {
        for (j in 1:length(ldots[[i]])) {
          nother <- nother + 1        
          x2[[nother+1]] <- ldots[[i]][[j]]
          nn <- names(ldots[[i]])[[j]]
          if (is.null(nn)) 
            nn <- nother +1
          names(x2)[nother+1] <- nn
        }
      } 
  if (is.null(names(x2)))
    names(x2) <- 1:length(x2)    
  ii <- which(is.na(names(x2)))
  if (length(ii) > 0)
    names(x2)[ii] <- ii
  list (nother = nother, x2 = x2)
}

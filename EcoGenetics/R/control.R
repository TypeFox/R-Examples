################################################################################
## CHECKPOINT PROGRAMS----------------------------------------------------------
################################################################################
# These functions are used as control points in other programs.

#------------------------------------------------------------------------------#
#' Check a connection network
#' @param con Connection network.
#' @author Leandro Roser \email{leandroroser@@ege.fcen.uba.ar} 
#' @keywords internal

int.check.con <- function(con) {
  
  ccon <- class(con)[1]
  
  if(ccon == "listw") {
    listwg <- sapply(con$neighbours, c, simplify = FALSE)
    weig <- sapply(con$weights, c, simplify = FALSE)
    Z<- 1:length(con$weights)
    wg <- outer(Z, Z)
    wg[] <- 0
    for(i in 1:nrow(wg)) {
      wg[i, ][listwg[[i]]] <- weig[[i]]
    }
  } else if(ccon == "matrix"){
    wg <- con
  } else if(ccon == "eco.weight"){ 
    wg <- con@W
  } else {
    stop("weight object provided is not of class listw, matrix or eco.weight")
  }
  wg
}

#------------------------------------------------------------------------------#
#' Check numeric format in a data frame
#' @param x Matrix or data frame. 
#' @author Leandro Roser \email{leandroroser@@ege.fcen.uba.ar}
#' @keywords internal

int.check.numeric <- function(mat) {
  
  x <- mat
  clases <- character()
  for(i in 1:ncol(x)) {
    clases[i] <- class(x[, i])
  }
  
  if(any(clases != "numeric" | clases != "integer")) {
    x <- as.matrix(x)
    colhier <- ncol(x)
    rowhier <- nrow(x)
    x <- matrix(as.numeric(x), ncol = colhier, nrow= rowhier)
    if(class(mat) == "data.frame") {
      x <- as.data.frame(x)
    }
    colnames(x) <- colnames(mat)
    rownames(x) <- rownames(mat)
    
    x
  }
}

#------------------------------------------------------------------------------#
#' Check row names
#' @param X Matrix or data frame.
#' @param lab String used as label.
#' @author Thibaut Jombart. Adapted by Leandro Roser. 
#' @keywords internal

int.check.rownames <- function(X, lab = "") {
  rnames <- rownames(X)
  rnames <- aue.rmspaces(rnames)
  if (is.null(rnames) || any(duplicated(rnames))){
    message("Note: null or duplicated row names. using generic labels.")
    rownames(X) <- aue.genlab(lab, nrow(X))
  }
  X
  
}


#------------------------------------------------------------------------------#
#' Check column names
#' @param X Matrix or data frame.
#' @param lab String used as label.
#' @author Thibaut Jombart. Adapted by Leandro Roser. 
#' @keywords internal

int.check.colnames <- function(X, lab = "L") {
  cnames <- colnames(X)
  cnames <- aue.rmspaces(cnames)
  if (is.null(cnames) || any(duplicated(cnames))){
    message("Note: null or duplicated column names. using generic labels.")
    colnames(X) <- aue.genlab(lab, ncol(X))
  }
  X
  
}


#------------------------------------------------------------------------------#
#' Check a vector of names
#' @param X Vector of names.
#' @param len.X Expected length of X.
#' @param lab String used as label.
#' @author Leandro Roser \email{leandroroser@@ege.fcen.uba.ar}
#' @keywords internal

int.check.vnames <- function(X, len.X, lab = "V") {
  X <- aue.rmspaces(X)
  if(any(duplicated(X)) || length(X) != len.X){
    message("Note: null or duplicated column names. using generic labels.")
    X <- aue.genlab(lab, len.X)
  }
  X
  
}


#------------------------------------------------------------------------------#
#' Check ploidy and number of digits per allele 
#' @param X Matrix to check.
#' @param ploidy Ploidy level in X.
#' @param ncod Number of digits coding each allele.
#' @param sep Character string separating alleles.
#' @param numeric.dat Numeric data checks. Default FALSE.
#' @author Leandro Roser \email{leandroroser@@ege.fcen.uba.ar}
#' @keywords internal


int.check.ncod <- function(X, ploidy, ncod = NULL, sep = "", numeric.dat = FALSE) {
  
  X <- as.matrix(X)
  mode(X) <- "character"
  
  #ploidy checks 
  
  if(ploidy < 1) {
    stop("ploidy can not be less than 1")
  }
  
  
  #control characters for numeric data
  if(numeric.dat) {
    sep.control <- gsub("[[:digit:]]", "", X)
    sep.control <- sep.control[!is.na(sep.control) & sep.control != sep]
    if(length(sep.control) != 0) {
      stop("non numeric (non-missing, non \"sep\") characters found
           with <numeric> option  = TRUE.
           Character data can be converted into numeric with the
           function \"eco.format\". See help(eco.format)")
    }
  }
  
  
  #---check ncod and ploidy-----------#
  
  X.sub <- gsub(meta2char(sep), "", X)
  X.sub <- X.sub[!is.na(X.sub)]
  n.control <- as.numeric(unique(nchar(X.sub)))
  
  ## more than one character length
  if(length(n.control) != 1) {
    stop("non unique character length found for alleles")
  }
  
  ## check that ncontrol(mod = ploidy) = 0
  if(n.control %% ploidy != 0) {
    stop("incongruence found between the number of (non-missing)
         characters in some cells and the ploidy level")
  }
  
  if(sep != "") {
    ## check that <sep> appears ploidy-1 times
    sep.rep <- gsub(paste("[^", meta2char(sep), "]", sep = ""), "", X)
    sep.rep <- sep.rep[!is.na(sep.rep)]
    sep.rep <- nchar(sep.rep) + 1
    if(any(sep.rep != ploidy)) {
      stop("incongruence between the number of alleles 
           determined by <sep> and the ploidy level 
           in some cells")
    }
    }
  
  # when ncod is NULL, determine its value using the ploidy and the number
  # of non "sep" characters.
  if(is.null(ncod)) {
    
    ncod <- n.control / ploidy 
    
  } else {
    
    if((n.control / ploidy)  != ncod) {
      stop(paste("all (non <0>) cells must have", "a length of", 
                 paste("(", ncod, ")", sep = ""), "non <sep> characters, but
                 seems to have", paste("(", n.control, ")", sep =  "")))
    }
    }
  
  ncod
  
  }


#------------------------------------------------------------------------------#
#' Check factor name consistency in a data frame and returns the corresponding column
#' @param X Matrix or data frame
#' @param lab String used as label
#' @author Thibaut Jombart. Adapted by Leandro Roser. 
#' @keywords internal

int.check.group <- function(X, grp = NULL, dummy = TRUE, exp.l = NULL) {
  
  X <- as.data.frame(X)
  
  
  #----basic control---#
  #check grp
  if(!is.null(grp)) {
    cond1 <- !is.character(grp) || length(grp) != 1
    if(cond1) {
      stop("invalid argument \"grp\" (non <character> or <null>, or length(grp) != 1)")
    }
  }
  #check X class
  cond2 <-  !is.matrix(X) && !is.data.frame(X)
  if(cond2) {
    stop("X is not of class <matrix> or <data.frame>")
  }
  #--------------------#
  
  # empty matrix or data frame
  if(any(dim(X) == 0)) {
    if(is.null(exp.l)) {
      stop("X is an object of dimension 0")
    }
    return(factor(rep(1, exp.l)))
  }
  
  #control the number of rows, if exp.l is passed
  if(!is.null(exp.l)) {
    if(exp.l != nrow(X)) {
      stop(paste("X has not <exp.l>", paste("(", exp.l, ")", sep = ""), "row(s)"))
    }
  }
  
  #else, exp.l is the output size, and exp.l == nrow(X)
  exp.l <- nrow(X)
  
  #if no group defined, return error or dummy variable
  if(is.null(grp)) {
    if(!dummy) {
      stop("no group defined")
    }
    return(factor(rep(1, exp.l)))
  }
  
  #group defined case--------
  pop <- colnames(X) %in% grp
  #control multiple matches
  if(sum(pop) > 1) {
    stop("grp matches multiple colnames of X")
  }
  #no match
  if(sum(pop) == 0) {
    if(!dummy) {
      return(NULL)
    } 
    #create a dummy variable
    dummy.fact <- factor(rep(1, exp.l))
    return(dummy.fact)
  }
  return(X[grp][[1]])
}


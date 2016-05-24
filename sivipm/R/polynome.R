###################################################################
# sivipm R package
# Copyright INRA 2016
# INRA, UR1404, Research Unit MaIAGE
# F78352 Jouy-en-Josas, France.
#
# URL: http://cran.r-project.org/web/packages/sivipm
#
# This file is part of sivipm R package.
# sivipm is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# See the GNU General Public License at:
# http://www.gnu.org/licenses/
#
###################################################################
#--------------------------------------------
# Class poly
# SLOTS:
# P: a list of length equal to the number of monomials.
#     Each component is a vector of length equal to the
#     polynomial degree: its contains the
#     number of the variables in the monomial
# indic: the indicatrice matrix. It is a matrix with as
#        many rows as monomials and nvar columns.
#       The element \code{(i, j)} is 1 when the variable \code{j}
#        is in the monomial, 0 otherwise.
# METHODS:
# bind (or bind.poly), expand (or expand.poly), print, show, summary
# CREATORS:
# new,  crpoly, crindic
# Internal
#--------------------------------------------
poly <- setClass("poly",
                     slots=c(P="list", indic="matrix"),
                 # default value
                 prototype = list(P=list(), indic=matrix(NA)),
          # Function to check:
           validity = function(object) {
             P=slot(object,"P")
             indic=slot(object,"indic")
             retour=
	        is.list(P)  && is.matrix(indic) && is.numeric(indic)
             if (!retour)
                return("P must be a numeric list and indic a numeric matrix")
             # Length of each component of P= number of monomials
             nvar <- ncol(indic)
             retour <- lapply(P, function(X, nvar) {
               return(  any(X<0) || any(X >nvar))
             }, nvar)
             if (any(unlist(retour)==TRUE))
               return(paste("All values in P must be positive and less than the number of X-input, i.e", nvar))
             degree <- length(P[[1]])
             retour <- lapply(P, function(X, degree) {
               return( (length(X) != degree)) }
                              , degree)
             if (any(unlist(retour)==TRUE))
               return(paste("The length of all components in P must be equal to the polynomial degree, i.e ", degree))
#
retour <- all(sapply(P[1:nvar], function(X) { a= (X>0); length(a[a==T]); })==1)
             if (!retour) {
               return("The first monomials should be the X-inputs")
             }
               

             return(TRUE)
           }
          )


#--------------------------------------------
#  Bind two objects of class "poly"
# @title Bind two structures coding for polynomial into a single one
# @param p1 an object of class "poly"
# @param p2 an object of class "poly"
# @return  an object of class "poly" coding for the polynomial catenation of the polys in argument
#--------------------------------------------
  

bindtwopoly <-            function(p1, p2) {

  nvar <- NULL
  P <- NULL
  indic <- NULL
  nbargs <- nargs()

  for (iarg in 1:nbargs) {
    argu <- eval(match.call()[[iarg+1]], envir = parent.frame()) # la valeur du i-eme argu

    if (class(argu) != "poly")       
stop(paste("bind: Bad argument ", argu, ": must be an object of class 'poly'"))

    
    if (!is.null(nvar) && (ncol(argu@indic)!=nvar))
      stop(paste("bind: argument", iarg, "is invalid: all polynomials must have the same number of variables, i.e:", nvar))
    
      if (iarg==1) {
        P <- argu@P
        varnames <- colnames(argu@indic)
        nvar <- ncol(argu@indic)
      } else {
        P <- c(P, argu@P)
      }

  
  # Set the polynomial degree to the maximal degree
  degree <- max(unlist(lapply(P, length)))
  P <- lapply(P, function(X, degree) {
    z <- c(X, rep(0, (degree-length(X))))
    z
  },degree)
#  cat("Polynomial degree: ", degree, "\n")

 # Remove the duplicated monomials from P
  Pz <- lapply(P, unname) # don't compare the variables names
  Pd <- duplicated(Pz)
  P <- P[!Pd]

  

                
  indic <- crindic(P, nvar)@indic
  colnames(indic) <- varnames
# Remove the duplicated monomials from indic
  if (any(Pd)) {
    indic <- indic[-which(Pd),, drop=FALSE]
  }
  } # end (iarg in 1:nbargs)

  lePoly <- poly(P=P, indic=indic)
             return(list(Pindic=lePoly, Pd=Pd ))
} # fin bindtwopoly



#--------------------------------------------
# Method 'print':
#--------------------------------------------
print.poly <- function (x, all=FALSE, ...) {
    # all=T : print of all the monomials. Otherwise, only the number
    # of monomials.
    if (all) {
  descr <- "" # description par nom de variables
  descrnum <- "" # description par no de variables
  
  label <- colnames(x@indic)
  if (length(x@P) >0) {
  for (imon in 1:length(x@P)) {
    for (ivar in 1:length(x@P[[imon]])) {
      if (x@P[[imon]][ivar] ==0)
          break
      descr <- paste(descr, label[ x@P[[imon]][ivar]], sep="")
      descrnum <- paste(descrnum, x@P[[imon]][ivar],  sep="")
      if (ivar < length(x@P[[imon]]) && (x@P[[imon]][ivar+1] !=0)) {
        descr <- paste(descr, "*", sep="")
        descrnum <- paste(descrnum, "*", sep="")
      }
    } # fin ivar
    if (imon < length(x@P)) {
      descr <- paste(descr, "+ ")
       descrnum <- paste(descrnum, "+ ")
    }
      # aller a la ligne pour pas que ce soit coupé n'importe ou
    quot <- nchar(descr) %/% options()$width
      if ( quot  != 0) {
        cat(descr,"\n" )
        descr <- ""
      }
    # pour descrnum, c'est un plus compliqué car on ne l'écrit
    # pas au fur et a mesure
    # trouver la position du dernier passage a la ligne
    icr <- gregexpr("\n",descrnum)[[1]]
    ipos <- icr[length(icr)]
    # trouver la longueur de la ligne courante
    lglig <- nchar(descrnum) - ipos
    quot <- lglig %/% options()$width
      if ( quot  != 0) {
        descrnum <- paste(descrnum, "\n")
      }
  } # fin imon
} else {
  cat("Polynomial unknown")
}
  
  cat(descr , "\n")
  cat("Polynomial description using variable numbers:\n")
  cat(descrnum , "\n")

# ##     
# ##              cat("*** Polynomial list: ***\n")
# ##            print(x@P, ...)
# ##            cat("*** Indicatrice matrix: ***\n")
# ##            print(x@indic, ...)
}
    if (length(x@P) >0) {
            cat("Polynomial degree: ", length(x@P[[1]]), "\n")
            cat("Number of monomials: ", length(x@P), "\n")
          }
    else {
      cat("There is no polynomial description\n")
    }
    if (!is.null(x@indic) && any(!is.na(x@indic))) {
            cat("Number of variables: ", ncol(x@indic), "\n")
          }
    
    return(invisible())
} # end print.poly
               
         
 setMethod("print",  signature(x="poly"),
           definition=print.poly)
#--------------------------------------------
# Method 'show':
#--------------------------------------------
 setMethod("show",  signature(object="poly"),
           function(object) {
           print(object)
         })

#--------------------------------------------
# Method 'summary'
#--------------------------------------------
summary.poly <- function(object, ...) {
  if (length(object@P) >0 ) {
      print.poly(object, all=TRUE)
} else {
  cat("Polynomial unknown\n")
}
  return(invisible())             
         }

 setMethod("summary",  signature(object="poly"),
           definition = summary.poly)
             
#--------------------------------------------
# Method 'expand':
#--------------------------------------------
 expand.poly <-          function(Pindic, dataX) {
# construction de la data.frame avec toutes les interactions prevues
# par le poly: creation d'un object of class 'polyX'
# @title   creation de la data.frame comportant toutes les interactions contenues dans un polynome.
# @expand
# @param Pindic object of class 'poly'
# @param dataX une data.frame (the X data not expanded)
# @return an object of class polyX

    # Verif des arguments
    if (class(Pindic)!="poly") 
        stop("expand: First argument must be an object 'poly'")
    if (!is.data.frame(dataX) && !is.matrix(dataX)) 
        stop("expand: Second argument must be a data.frame or a matrix")
    if (length(Pindic@P) ==0)
        stop("expand: First argument is an unknown polynomial\n")
    if (max(unlist(Pindic@P))>ncol(dataX)) 
      stop("expand: the polynomial cannot have values greater than the maximum number of variables in the data")
    P <- Pindic@P
    dataX.exp <- matrix(nrow = nrow(dataX), ncol = length(P))
    for (i in 1:length(P)) {
        dataX.trans <- dataX[, P[[i]]]
        if (sum(as.integer(P[[i]] > 0)) > 1) {
            dataX.exp[, i] <- apply(dataX.trans, 1, prod)
        } else {
            dataX.exp[, i] <- dataX.trans
        }
        
    }
    onelig <- function(x, label, ret) {
  for (ix in 1:length(x)) {
        unx <- x[ix]
        if (ix>1)
          lesep="*"
        else
          lesep=""
        
        if (unx>0) ret <- paste(ret, label[unx], sep=lesep)
      }
      ret
} # end onelig 
    ret<-""
    colnames(dataX.exp) <- lapply(P, onelig, colnames(dataX), ret)
    
    rownames(dataX.exp) <- seq(1, nrow(dataX))
    colnames(Pindic@indic) <- colnames(dataX)
    lePolyX <- new("polyX",Pindic=Pindic, dataX.exp=as.data.frame(dataX.exp))
    return(lePolyX)

} # end expand
  
setGeneric("expand",
           function(Pindic, dataX){ value <- standardGeneric("expand") })

setMethod("expand",  signature(Pindic="poly"),
          definition=expand.poly)



# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#--------------------------------------------
# Class polyX: a poly object + the extended data
# SLOTS:
# Pindic: an object of class 'poly'
# dataX.exp: the extended data
# METHODS:
# summary, bind, print
# CREATORS:
# new, crpolyX,crpolyXT,  expand.poly

polyX <- setClass("polyX",
           slots=c(dataX.exp="data.frame", Pindic="poly"),
                  # default value
          # Function to check:
           validity = function(object) {
             if (length(object@Pindic@P) >0) {
             Pindic=slot(object,"Pindic")
             dataX.exp=slot(object,"dataX.exp")
             retour=
               (ncol(dataX.exp) >= length(Pindic@P))
             if (!retour) {
                stop("polyX creator: Polynomial description have more monomials as expanded data columns")
                return(FALSE)
              }
             # Cannot have variables with one level
             nvar <- ncol(Pindic@indic)
             retour <-any(apply(as.data.frame(dataX.exp[, 1:nvar]), 2,
                                function(X) {length(unique(X))})==1)
             if (retour) {
                stop("polyX creator: Some X-inputs have only one value")
               return(FALSE)
              } 

           } # fin is.null
                        return(TRUE)
           }
          )

#--------------------------------------------
# Method 'bind'
#--------------------------------------------
bind.polyX <-            function(x, ...) {
  nbargs <- nargs()
  # verif
  for (i in 1:nbargs) {
    argu <- eval(match.call()[[i+1]], envir = parent.frame()) # la valeur du i-eme argu

    if (class(argu) != "polyX")       
stop(paste("bind: Bad argument ", argu, ": must be an object of class 'polyX'"))
    
     if (i==1) {
       lePolyX <- argu
     } else {
        dataX.exp <- lePolyX@dataX.exp
        Pindic <- lePolyX@Pindic
       if (nrow(argu@dataX.exp) != nrow(dataX.exp)) {
         stop("bind: all arguments must have the same number of observations\n")
       }
       dataX.exp <- cbind(dataX.exp, argu@dataX.exp)
       ret <- bindtwopoly(Pindic, argu@Pindic)
       # Oter les colonnes de dataX.exp qui correspondent aux monomials en double
       if (any(ret$Pd)) {
         dataX.exp <- dataX.exp[, -which(ret$Pd), drop=FALSE]
       }
        lePolyX <- polyX(dataX.exp=dataX.exp, Pindic=ret$Pindic)
      } # fin i !=1

    
  } # fin (i in 1:nbargs)
return(lePolyX)
} # fin bind.polyX

setGeneric("bind",
           function(x,...){ value <- standardGeneric("bind") })

setMethod("bind", signature(x="polyX"),
           definition=bind.polyX)



  
#--------------------------------------------
# Method 'print'
#--------------------------------------------
print.polyX <- function(x, all=FALSE, ...) {

 
  if (!is.null(x@Pindic))
  print(x@Pindic, all)
  else
    cat("unknown")
#  cat("Dimension of the expanded data-frame:\n")
#  print(dim(x@dataX.exp))
# Test if the polynomial description is given.
  # If it is not given, the 'poly' object, Pindic, does not
  # know the number of variables.
  if (is.null(x@Pindic@indic) || all(is.na(x@Pindic@indic))) {
            cat("Number of variables: ", ncol(x@dataX.exp), "\n")
          }

  cat("Number of observations:", nrow(x@dataX.exp), "\n")
  return(invisible())
}
setMethod("print", signature(x="polyX"),
           definition=print.polyX)

#--------------------------------------------
# Method 'show':
#--------------------------------------------
 setMethod("show",  signature(object="polyX"),
           function(object) {
           print(object)
         })


#--------------------------------------------
# Method 'summary'
#--------------------------------------------
summary.polyX <- function(object, ...) {
             summary(object@Pindic)
            cat("Number of observations: ", nrow(object@dataX.exp), "\n")
 return(invisible())              
         }

 setMethod("summary",  signature(object="polyX"),
           definition=summary.polyX)

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Creators of objects "poly"
#--------------------------------------------

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# @crpoly
#  Build an object of class 'poly'
# @title Build an object of class 'poly' (list coding for a multivariate polynomial and indicatrice matrix)
# @param nvar number of variables in the polynomial
# @param d  maximal degree of the variables  in the polynomial
# @param type type of the required polynomial. Character string among:
# "full": complete polynomial with all terms
# "power": only quadratic terms (d=2), cubic (d=3), etc..
# "interact":  only interactions terms of degree d (not included quadratic, cubic, ... terms)
# @return   an object of class 'poly'
crpoly <- function(nvar, d, type="full") {

  # Build the polynomial list
  switch(type,
         full = {
    P <- polynomecomplet(nvar, d)
 # P est une liste; chaque elt est un vecteur de longueur <=d
  },
         
         power = { P= polysquares(nvar, d)},
         interact = { P= polyinteract(nvar, d)},
         stop("crpoly: invalid argument 'type'. Valid values are 'full', 'power, 'interact'")
       ) # end switch
  # labeller les composants de P
  P <- lapply(P, function(X) { names(X) <- paste("V", 1:length(X), sep="")
                return(X) })
  
    # Build the   indicatrice matrix
  retour <- crindic(P,nvar)
  return(retour)
  
} # end crpoly
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# Build the indicatrice matrix useful to calculate the TSIVIP
# @title build the indicatrice matrix
# @param P the list coding a multivariate polynomial
# @param nvar number of variables in the model
# @return an object of class 'poly'
# Dimension:  number of monomials x nvar 

crindic <- function (P, nvar) {
  indic1 <- function(i, j, P) {
    a <- (j == P[[i]])
    a <- as.integer(a)
    a <- sum(a)
    a <- (a > 0)
    a <- as.integer(a)
    return(a)
} # end indic1 

  if (!is.list(P)) {
       stop("crindic: the first argument must be a list ")
        } else {
          # Verifier que les nos de variables dans P <= nvar
          if (any(sapply(P, function(x) any(x>nvar))))
              stop("crindic: values in the polynomial list greater than the number of variables") 
          
indic <- matrix(nrow = length(P), ncol =  nvar)
    for (i in 1:nrow(indic)) {
        for (j in 1:ncol(indic)) {
            indic[i, j] <- indic1(i, j, P)
        }
    }

lePoly <- poly(P=P, indic=indic)
  return(lePoly)
}
} # fin crindic
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Creators of objects "polyX"
#--------------------------------------------
# @crpolyX
#  Build an object of class 'polyX'
# @title Build an object of class 'polyX' (list coding for a multivariate polynomial and indicatrice matrix+ expanded data)
# @param dataX data.frame (not expanded)
# @param d  maximal degree of the variables  in the polynomial
# @param type type of the required polynomial. Character string among:
# "full": complete polynomial with all terms
# "power": only quadratic terms (d=2), cubic (d=3), etc..
# "interact":  only interactions terms of degree d (not included quadratic, cubic, ... terms)
# @return   an object of class 'polyX'
crpolyX <- function( dataX,  d, type="full") {
  if (!is.data.frame(dataX))
    stop("crpolyX: the first argument must be a data.frame")
  nvar <- ncol(dataX)
 Pindic <- crpoly( nvar, d, type)
  lePolyX <- expand(Pindic, dataX)
  return(lePolyX)
} # end crpolyX
#--------------------------------------------
# @crpolyXT
#  Build an object of class 'polyX'
# @title Build an object of class 'polyX' (list coding for a multivariate polynomial and indicatrice matrix+ expanded data)
# @param dataXT data.frame ( expanded)
# @param varnames names of the inputs ("raw" variables)
# @param d  maximal degree of the variables  in the polynomial
# @param type type of the required polynomial. Character string among:
# "full": complete polynomial with all terms
# "power": only quadratic terms (d=2), cubic (d=3), etc..
# "interact":  only interactions terms of degree d (not included quadratic, cubic, ... terms)
# @return   an object of class 'polyX'
crpolyXT <- function(varnames, dataXT, d, type="full") {
  if (!is.data.frame(dataXT))
    stop("crpolyXT: second argument must be a data.frame")
  nvar <- length(varnames)
 Pindic <- crpoly( nvar, d, type)
  colnames(Pindic@indic) <- varnames
  lePolyX <- new("polyX", dataX.exp=dataXT, Pindic=Pindic)
  return(lePolyX)
} # end crpolyXT

#--------------------------------------------
# Remove monomials from an object of class 'polyX'
takeoff.polyX <- function(P, monomials) {
  if (is.character(monomials)) {
    indices <- match(monomials, colnames(P@dataX.exp))
    if (any(is.na(indices))) {
      faux <- which(is.na(indices))
      stop(paste("takeoff: The monomial numbers,",  faux, ", are not in the polynomial\n"))
    }
  } else {
    if (any(monomials > ncol(P@dataX.exp)))
      stop("takeoff: Some monomials are not in the polynomial")
    indices <- monomials
  }
  if (length(indices) >= ncol(P@dataX.exp))
     stop("takeoff: Too many monomials are to be removed")
  
  # Remove from dataX.exp
  dataX.exp <- P@dataX.exp[, -indices, drop=FALSE]
    # Remove from P
  PP <- P@Pindic@P[-indices]
  # Remove from indic
  indic <- P@Pindic@indic[ -indices, , drop=FALSE]
  lePoly <- new("poly", P=PP, indic=indic)
  lePolyX <- new("polyX",dataX.exp=dataX.exp, Pindic=lePoly)
   return(lePolyX)
} # end takeoff.polyX



setGeneric("takeoff",
           function(P, monomials){ value <- standardGeneric("takeoff") })
 setMethod("takeoff",  signature(P="polyX"),
           definition=takeoff.polyX)


# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  
# Build a list coding for a multivariate polynomial which only includes
# interactions of degree d
# Internal function
polyinteract<- function(nvar, d) {
  if (d<2)
    stop("polyinteract: degree should be greater than 1")
  
    P <- NULL
  ret <- polynomecomplet(nvar, d)
    # ret est une liste
  # Oter les squares
  otesq <- function(X) {
    if (all(X == X[1]))
      return( NULL)
    else 
      return( X)
  }
  retour <- lapply(ret,  otesq)
    # retour est une liste
  
    # Garder les interactions 
    gardeinteract<- function(X) {
      if (length(unique(X[X>0])) >1 || length(X[X>0])==1 )
        return(X)
      else
      return (NULL)
    }

    retour <- lapply(retour, gardeinteract)
    # Oter les elements nuls

    s <- length(retour)
    j <- 1
    for (i in 1:s) {
      if (!is.null(retour[[i]])) {
        P[[j]] <- retour[[i]]
        j <- j+1
      }
    }
    
    return(P)
  
} # end polyinteract
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Build a list coding for a multivariate polynomial which only includes
# squares (cubic, ...) terms
# Internal function
polysquares<- function(nvar, d) {
  retour <- list()
  l <- 1
  for (j in 1:d) {
  for (i in 1:nvar) {
      tab <- rep(i, j)
      retour[[l]] <- c(tab, rep(0, (d-j)))
      l <- l+1
    }
}
  return(retour)
} # end polysquares


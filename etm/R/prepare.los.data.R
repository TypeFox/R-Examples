"prepare.los.data" <-
function(x) {
## --------------------------------------------------------------------------------
## Title: R-function prepare.los.data()
## ---------------------------------------------------------------------------------
## Author: Matthias Wangler
##         mw@imbi.uni-freiburg.de
##         Institute of Med. Biometry and Med. Computer Science
##         Stefan-Meier-Strasse 26, D-79104 Freiburg,
##         http://www.imbi.uni-freiburg.de
## ---------------------------------------------------------------------------------
## Description: Read and prepare a data set which can be passed to the function clos
## ---------------------------------------------------------------------------------
## Required Packages: -
## ---------------------------------------------------------------------------------
## Usage: prepare.los.data( x )
##
## x: data.frame of the form data.frame( id, j.01, j.02, j.03, j.12, j.13, cens):
##
## id:      id (patient id, admision id, ...)
## j.01:    observed time for jump from "0" to "1"
## j.02:    observed time for jump from "0" to "2"
## j.03:    observed time for jump from "0" to "3"
## j.12:    observed time for jump from "1" to "2"
## j.13:    observed time for jump from "1" to "3"
## cens:    observed time for censoring
## ---------------------------------------------------------------------------------
## Value: data.frame of the form data.frame(id, from, to, time ):
##
##        id:   id (patient id, admision id)
##        from: the state from where a transition occurs
##        to:   the state to which a transition occurs
##        time: the time a transition occurs
##        oid:  the observation id
## ---------------------------------------------------------------------------------
## Notes: It's possible that the same patient, person or object was observed several
##        times (e.g. bootstrap).
##        So for each observation the same id recieves different observation id's. 
## ---------------------------------------------------------------------------------
## Example: > data(los.data)
##          > my.observ <- prepare.los.data(x=los.data)
## ---------------------------------------------------------------------------------
## License: GPL 2
##----------------------------------------------------------------------------------
## History:    20.06.2004, Matthias Wangler
##                         first version
## ---------------------------------------------------------------------------------
  
  ## check the passed parameters
  if( missing(x) )
  {
    stop("Argument 'x' is missing, with no defaults.")
  }
  if( !is.data.frame(x) )
  {
    stop("Argument 'x' must be a 'data.frame'.")
  }

  ## check the number of columns of the passed data.frame x
  if( dim(x)[2] != 7 )
  {
    stop("The passed data.frame 'x' doesn't include 7 columns.")      
  }

  ## compute variables cens.0 for admissions censored in the initial state 0
  ## and cens.1 for admissions censored in state 1

  x$cens.0 <- x$cens
  x$cens.0[is.finite(x[,2])] <- Inf
    
  x$cens.1 <- x$cens
  x$cens.1[is.infinite(x[,2])] <- Inf

    
  x <- x[,c(1,2,3,4,5,6,8,9)]

  
  id   <- c(x[,1][x[,2] != Inf], x[,1][x[,3] != Inf],x[,1][x[,4] != Inf],
            x[,1][x[,5] != Inf], x[,1][x[,6] != Inf],x[,1][x[,7] != Inf],
            x[,1][x[,8] != Inf])
  
  from <- c(rep("0",length(x[,2][x[,2] != Inf])), rep("0",length(x[,3][x[,3] != Inf])),
            rep("0",length(x[,4][x[,4] != Inf])), rep("1",length(x[,5][x[,5] != Inf])),
            rep("1",length(x[,6][x[,6] != Inf])), rep("0",length(x[,7][x[,7] != Inf])),
            rep("1",length(x[,8][x[,8] != Inf])))

  to   <- c(rep("1",length(x[,2][x[,2] != Inf])), rep("2",length(x[,3][x[,3] != Inf])),
            rep("3",length(x[,4][x[,4] != Inf])), rep("2",length(x[,5][x[,5] != Inf])),
            rep("3",length(x[,6][x[,6] != Inf])), rep("cens",length(x[,7][x[,7] != Inf])),
            rep("cens",length(x[,8][x[,8] != Inf])))

  time <- c(x[,2][x[,2] != Inf], x[,3][x[,3] != Inf],x[,4][x[,4] != Inf],
            x[,5][x[,5] != Inf], x[,6][x[,6] != Inf],x[,7][x[,7] != Inf],
            x[,8][x[,8] != Inf])  

  ## observation id
  x$oid <- 1:length(x[,1])

  oid   <- c(x[,9][x[,2] != Inf], x[,9][x[,3] != Inf],x[,9][x[,4] != Inf],
             x[,9][x[,5] != Inf], x[,9][x[,6] != Inf],x[,9][x[,7] != Inf],
             x[,9][x[,8] != Inf])

  observ <- data.frame(id, from, to, time, oid)
  
  return(observ)
}


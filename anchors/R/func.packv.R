## File:   func.pack.R
## Author  : Jonathan Wand <jwand@latte.harvard.edu>
## Date:   2002-05-30
##
## Created :  2002-05-30
## Modified:  2003-06-14 : made rm() target only things that exist
## Modified:  $Date: 2005/08/10 19:06:50 $
## Revision:  $Revision: 1.1 $
## RCS-ID:    $Id: func.packv.R,v 1.1 2005/08/10 19:06:50 jwand Exp $
##
## Overview:
##
## The goal of these functions is to be able to pass in
## a single list() of objects which contain small things we find easy to handle
## (individual vectors of names and parameters),
## and have this translated accurately to a single vector
## which optim(), genound() and other functions like to use,
## AND be able to translate back again, without thinking about cutting up
## the vectors ourselves each time.
##
## packv() takes the numerical values in omega$start (input)
## and translates it to a single vector omega$pvec
##
## unpackv() takes omega$pvec as input and returns a revised omega$start
##


## Functions: packv()
## Authors:   Jonathan Wand and Gary King
## Date:      2002-05-31
##
## Purpose:  make a list of parameters a vector
##
## Input:
##   list() with following members:
##     labels: list of vectors, containing parameter names (strings)
##     start: list of vectors, containing parameter values (numerical)
##     estimated: list of vectors, containing flags to make parms static/exclude (logical)
##
## NOTE: labels and start must have member vectors with matching length
##       estimated can omit vectors, or assign NULL, if nothing is excluded
##       
## Output:
##   list() with members as input PLUS:
##     pvec : vector representation of start -- this is result of PACKING
##     nvec : vector representation of labels -- this is result of PACKING
##   and some stuff for later unpacking...
##     len  : vector of lengths of vectors in start
##     parse: expression() which will be used to unpack 'v' into start!
##      
packv <- function(omega) {
  nn <- names(omega$labels)
  dnames <- names(omega$estimated) 

  omega$nitems <- length(nn)
  ## clean out stuff we are going to update
  omega$parse    <- omega$len      <- omega$pvec     <- omega$nvec <- NULL
  omega$all.pvec <- omega$all.nvec <- omega$all.fvec <- NULL
  
  for (i in 1:omega$nitems) {
    name <- nn[i]
    len.p <- eval(parse(text=paste("length(omega$start$",name,")",sep="")))
    len.n <- eval(parse(text=paste("length(omega$labels$",name,")",sep="")))
    len.d <- eval(parse(text=paste("length(omega$estimated$",name,")",sep="")))
    if (len.p!=len.n) {
      cat("Parameters and Names of unequal length for",
          name,":",len.p,"!=",len.n,"\n")
      return(NULL)
    } else {
      tmp <- NA
      eval(parse(text=paste("tmp <- length(omega$estimated$",name,")",sep="")))
      ## do we have exclusions? 
      if (!any(name == dnames) || tmp < 1) {
        ## NO: 
        dtmp <- rep(TRUE,len.p)
      } else {
        ## possibly:
        dtmp <- eval(parse(text=paste("omega$estimated$",name,sep="")))
      }

      omega$all.pvec <-
        c(omega$all.pvec,eval(parse(text=paste("omega$start$",name,sep=""))))
      omega$all.nvec <-
        c(omega$all.nvec,eval(parse(text=paste("omega$labels$",name,sep=""))))
      omega$all.fvec <- c(omega$all.fvec , dtmp)

      ## take out stuff from beta if necessary...
      if (any(dtmp==FALSE)) {
        omega$len   <- c(omega$len, sum(dtmp))
        dtmps <- paste(dtmp,collapse=",")

        txt   <-
          paste("omega$start$",name,"[c(",dtmps,")] <- tmp.extract",sep="" )
        omega$pvec <-
          c(omega$pvec,eval(parse(text=paste("omega$start$",name,"[c(",dtmps,")]",sep=""))))
        omega$nvec <-
          c(omega$nvec,eval(parse(text=paste("omega$labels$",name,"[c(",dtmps,")]",sep=""))))

        
      } else {
        omega$len   <- c(omega$len, len.p)
        txt <- paste("omega$start$",name,"<- tmp.extract",sep="" )
        omega$pvec <-
          c(omega$pvec,eval(parse(text=paste("omega$start$",name,sep=""))))

        omega$nvec <-
          c(omega$nvec,eval(parse(text=paste("omega$labels$",name,sep=""))))

      }
      omega$parse <- c(omega$parse,
                       parse(text=txt))
    }
  }
  return(omega)
}

## Functions: unpackv()
## Authors:   Jonathan Wand and Gary King
## Date:      2002-05-31
##
## Purpose:  make a vector a list of parameters
##           eg., can attach list ( $start ) to get parameters
##           for different parts of a likelihood function
##
## Input:
##   omega: list() with following members
##     pvec : vector representation of start
##     len  : vector of lengths of vectors in start
##     parse: expression() which will be used to unpack 'v' into start!
##     AND (optionally) you could also have this stuff...
##     labels: list of parameter names 
##     start: list of parameter values (which will overwriten)
##
## Output:
##   list() with members as above, plus revised value of:
##     start: list of parameter values -- this is result of UNPACKING
##

unpackv <- function(omega) {

  cur.max <- 0
  for (i in 1:omega$nitem) {
    if (omega$len[i] > 0) {
      idx <- c(1:omega$len[i])+cur.max
      cur.max <- cur.max+omega$len[i]
      tmp.extract <- omega$pvec[idx]
      eval(omega$parse[i])
    }
  }
  return(omega)
}

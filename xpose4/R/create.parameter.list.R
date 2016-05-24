# Xpose 4
# An R-based population pharmacokinetic/
# pharmacodynamic model building aid for NONMEM.
# Copyright (C) 1998-2004 E. Niclas Jonsson and Mats Karlsson.
# Copyright (C) 2005-2008 Andrew C. Hooker, Justin J. Wilkins, 
# Mats O. Karlsson and E. Niclas Jonsson.
# Copyright (C) 2009-2010 Andrew C. Hooker, Mats O. Karlsson and 
# E. Niclas Jonsson.

# This file is a part of Xpose 4.
# Xpose 4 is free software; you can redistribute it and/or
# modify it under the terms of the GNU Lesser General Public License
# as published by the Free Software Foundation, either version 3
# of the License, or (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.

# You should have received a copy of the GNU Lesser General Public License
# along with this program.  A copy can be cound in the R installation
# directory under \share\licenses. If not, see http://www.gnu.org/licenses/.

"create.parameter.list"  <- function(listfile)
{
  ## Read list file
  if(is.readable.file(listfile)) {
    x <- read.lst(listfile)
    if(length(x)==0 && x == 0) {
      cat("The output file does not contain any information")
      return()
    }
  } else {
    cat("The output file couldn't be found in the current directory.\n")
    return()
  }

  ## Count the number of parameters
  npar.list   <- calc.npar(x)

  #attach(npar.list, warn.conflicts=F)

  npar <- npar.list$npar
  nth <- npar.list$nth
  nseth <- npar.list$nseth
  nom <- npar.list$nom
  nseom <- npar.list$nseom
  nsi <- npar.list$nsi
  nsesi <- npar.list$nsesi



  seenterm <- seenobj <- seenth <- seenom <- seensi <- seenseth <- seenseom <-
    seensesi <- seennth <- seennom <- 0
  #attach(x, warn.conflicts=F)

  term <- x$term
  ofv <- x$ofv
  thetas <- x$thetas
  omega <- x$omega
  sigma <- x$sigma
  sethetas <- x$sethetas
  seomegas <- x$seomegas
  sesigmas <- x$sesigmas

  
  if(!any(is.null(term)))
    seenterm <- 1
  if(!any(is.null(ofv)))
    seenobj <- 1
  if(!any(is.null(thetas)))
    seenth <- 1
  if(!any(is.null(omega)) && nom !=0)
    seenom <- 1
  if(!any(is.null(sigma)) && nsi !=0)
    seensi <- 1
  if(!any(is.null(sethetas)))
    seenseth <- 1
  if(!any(is.null(seomegas)) && nseom !=0)
    seenseom <- 1
  if(!any(is.null(sesigmas)) && nsesi !=0)
    seensesi <-1

  ## Add parameters

  ## Construct the parameter information that is to be added to the screen
  ## This will be in the form of 1 or two arrays, the first being the array
  ## with parameter estimates and the second with SEs
  ## The nonparametrics are added to the parameter array afterwards
  ## The SE array will be constructed wether or not there are SE but filled
  ## with 0.

  ## If we haven't got any parameters in the lst-file
  if(seenth != 1) {
    cat("No parameters found, check your NONMEM output file.\n")
    return()
  }

  ## If we have thetas (we should have them!)
  if(seenth == 1) {
    nth           <- length(thetas)
    parnam        <- 0
    parval        <- 0
    parnam[1:nth] <- paste("TH",1:nth,sep="")
    parval[1:nth] <- format.default(thetas,digits=2)
  }

  ## If we have SE for the thetas
  if(seenseth == 1) {
    nseth           <- length(sethetas)
    separnam        <- 0
    separval        <- 0
    separnam[1:nseth] <- paste("RSE TH",1:nseth,sep="")

    ## To avoid division by zero or NA
    selzero           <- thetas == 0
    selna             <- is.na(sethetas)
    sel <- !selzero & !selna
    sel1              <- !sel
    cvthetas          <- 0
    cvthetas[sel] <- format.default(sethetas[sel]/abs(thetas[sel]),digits=2)
    cvthetas[sel1]    <- ""
    separval[1:nseth] <- cvthetas
  } else {
    nth             <- length(thetas)
    separnam        <- 0
    separval        <- 0
    separnam[1:nth] <- paste("RSE TH",1:nth,sep="")
    separval[1:nth] <- 0
  }

  ## If we have omegas
  if(seenom == 1 ) {
    nomega <- length(omega)
    for(i in 1:nomega) {
      sel <- omega[[i]] != 0
      if(nomega == 1) sel <- T
      if(i == 1) sel <- T
      if(length(omega[[i]][sel]) == 1 || i == 1) {        # Must be a diagonal omega

	parnam[length(parnam)+1]  <- paste("OM",i,":",i,sep="")

        ## If the first omega is fixed to 0
	if(omega[[i]][sel] == 0) {
          parval[length(parval)+1]  <-
            format.default(0,digits=2)
          ##separnam[length(separnam)+1]  <- paste("CT OM",i,":",i,sep="")
          separval[length(separval)+1]  <- ""
        } else {

          parval[length(parval)+1]  <-
            format.default(sqrt(omega[[i]][sel]),digits=2)
          separnam[length(separnam)+1]  <- paste("CT OM",i,":",i,sep="")
          separval[length(separval)+1]  <- ""
        }


      } else { # There are off-diagonals or the whole row is zero

	## The diagonal element
	parnam[length(parnam)+1]  <- paste("OM",i,":",i,sep="")
	parval[length(parval)+1]  <-
	  format.default(sqrt(omega[[i]][i]),digits=2)

	separnam[length(separnam)+1]  <- paste("CT OM",i,":",i,sep="")
	separval[length(separval)+1]  <- ""

        ## Loop over the off-diagonals
	for(j in 1:(length(omega[[i]])-1)) {

	  if(omega[[i]][j] == 0) next
	  if(omega[[i]][j] != 0) {
	    parnam[length(parnam)+1]  <- paste("OM",i,":",j,sep="")
	    parval[length(parval)+1]  <-
	      format.default(omega[[i]][j]/sqrt(omega[[i]][i]*omega[[j]][j]),digits=2)

	    separnam[length(separnam)+1]  <- paste("CT OM",i,":",j,sep="")
	    separval[length(separval)+1]  <- ""
	  }
	}
      }
    }
  }

  ## If we have SE for the omegas -- fill in their values
  if(seenseom == 1) {
    n <- length(thetas)           # n - is a flag that makes the values appear
                                  #     in the right place
    nseomega <- length(seomegas)
    for(i in 1:nseomega) {
      ## Select the non-zero omegas
      sel <- omega[[i]] != 0
      
      if(nseomega == 1) sel <- T
      if(i == 1) sel <- T

      if(length(seomegas[[i]][sel]) == 1) {
        ## Must be a diagonal omega

        n <- n+1
        separnam[n] <-
          paste("RSE OM",i,":",i,sep="")

        if(omega[[i]][sel] == 0) { ## If first omega is fixed to 0
          separval[n]  <- ""
        } else {
          if(is.na(seomegas[[i]][sel])){
            separval[n]  <- ""
          } else {
            if(seomegas[[i]][sel] == 0){
              separval[n]  <- ""
            } else {
              separval[n]  <-
                format.default(seomegas[[i]][sel]/abs(omega[[i]][sel]),digits=2)
            }
          }
        }

      } else {
        ## There are off-diagonals
        n <- n+1
	separnam[n]  <- paste("RSE OM",i,":",i,sep="")
        if(is.na(seomegas[[i]][i])){
          separval[n]  <- ""
        } else {
          if(seomegas[[i]][i] == 0){
            separval[n]  <- ""
          } else {
            separval[n]  <-
              format.default(seomegas[[i]][i]/abs(omega[[i]][i]),digits=2)
          }
        }

        ## Loop over the off-diagonals
	for(j in 1:(length(seomegas[[i]])-1)) {
	  if(omega[[i]][j] == 0) next
	  if(omega[[i]][j] != 0) {
	    n <- n +1
	    separnam[n]  <-
	      paste("RSE OM",i,":",j,sep="")
            if(is.na(seomegas[[i]][j])){
              separval[n]  <- ""
            } else {
              if(seomegas[[i]][j] == 0){
                separval[n]  <- ""
              } else {
                separval[n]  <-
                  format.default(seomegas[[i]][j]/abs(omega[[i]][j]),digits=2)
              }
            }
	  }
	}
      }
    }
  }
  
  ## capture the length of the parameter vector before sigma added
  th.om.par.length <- length(parnam)

  ## If we have sigmas
  if(seensi == 1) {
      nsigma <- length(sigma)
      for(i in 1:nsigma) {
      sel <- sigma[[i]] != 0

      if(length(sigma[[i]][sel]) == 1) {        # Must be a diagonal omega

	parnam[length(parnam)+1]  <- paste("SI",i,":",i,sep="")
	parval[length(parval)+1]  <-
	  format.default(sqrt(sigma[[i]][sel]),digits=2)

	separnam[length(separnam)+1]  <- paste("RSE SI",i,":",i,sep="")
	separval[length(separval)+1]  <- ""

      } else {                                  # There are off-diagonals

	## The diagonal element
	parnam[length(parnam)+1]  <- paste("SI",i,":",i,sep="")
	parval[length(parval)+1]  <-
	  format.default(sqrt(sigma[[i]][i]),digits=2)

	separnam[length(separnam)+1]  <- paste("RSE SI",i,":",i,sep="")
	separval[length(separval)+1]  <- ""

        ## Loop over the off-diagonals
	for(j in 1:(length(sigma[[i]])-1)) {
          if(sigma[[i]][j] == 0) next
	  if(sigma[[i]][j] != 0) {
	    parnam[length(parnam)+1]  <- paste("SI",i,":",j,sep="")
	    parval[length(parval)+1]  <-
	      format.default(sqrt(sigma[[i]][j]),digits=2)

	    separnam[length(separnam)+1]  <- paste("RSE SI",i,":",j,sep="")
	    separval[length(separval)+1]  <- ""
	  }
	}
      }
    }
  }

  ## If we have SE for the sigmas -- fill in their values
  if(seensesi == 1) {
    n <- th.om.par.length

    nsesigma <- length(sesigmas)
    for(i in 1:nsesigma) {
      sel <- sigma[[i]] != 0

      if(nsesigma == 1) sel <- T

      if(length(sesigmas[[i]][sel]) == 1) {
        ## Must be a diagonal omega
        n <- n+1
        separnam[n]<- paste("RSE SI",i,":",i,sep="")
        if(sigma[[i]][sel] == 0) {
          ## If first sigma is fixed to 0
          separval[n]<- ""
        } else {
          if(is.na(sesigmas[[i]][sel])){
            separval[n]  <- ""
          } else {
            if(sesigmas[[i]][sel] == 0){
              separval[n]  <- ""
            } else {
              separval[n]<-
                format.default(sesigmas[[i]][sel]/abs(sigma[[i]][sel]),digits=2)
            }
          }
        }
      } else {
        ## There are off-diagonals
        n <- n+1
        separnam[n]<- paste("RSE SI",i,":",i,sep="")
        if(is.na(sesigmas[[i]][sel])){
          separval[n]  <- ""
        } else {
          if(seomegas[[i]][sel] == 0) {
            separval[n]  <- ""
          } else {
            separval[n]<-
              format.default(sesigmas[[i]][i]/abs(sigma[[i]][i]),digits=2)
          }
        }
        ## Loop over the off-diagonals
	for(j in 1:(length(sesigmas[[i]])-1)) {

	  if(sigma[[i]][j] == 0) next
	  if(sigma[[i]][j] != 0) {
	    n <- n+1
	    separnam[n] <-
	      paste("RSE SI",i,":",j,sep="")
            if(is.na(sesigmas[[i]][j])){
              separval[n]  <- ""
            } else {
              if(sigma[[i]][j] == 0){
                separval[n]  <- ""
              } else {
                separval[n] <-
                  format.default(sesigmas[[i]][j]/abs(sigma[[i]][j]),digits=2)
              }
            }
	  }
	}
      }
    }
  }


  ret.list <- list(term = term, ofv = ofv,
                   seenterm = seenterm,
                   seenobj = seenobj,
                   seenth = seenth,
                   seenom = seenom,
                   seensi = seensi,
                   seenseth = seenseth,
                   seenseom = seenseom,
                   seensesi = seensesi,
                   npar = npar,
                   parnam = parnam,
                   parval = parval,
                   separnam = separnam,
                   separval = separval
                   )



  #detach(x)
  #detach(npar.list)

  return(ret.list)

}

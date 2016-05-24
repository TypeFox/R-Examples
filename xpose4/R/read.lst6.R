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

read.lst6 <- function(filename) {

  ## Function to split character strings
  string2num <- function(x)
    {
      oldopts <- options(warn = -1)
      on.exit(options(oldopts))
      nc <- nchar(x)
      tmp <- substring(x, 1:nc, 1:nc)
      spc <- tmp == " "
      st <- !spc & c(T, spc[ - nc])
      end <- !spc & c(spc[-1], T)
      as.numeric(substring(x, (1:nc)[st], (1:nc)[end]))
    }

  listfile <- scan(filename, sep = "\n", what = character(),quiet=TRUE)	
  
  ## Find termination messages
  minim <- pmatch("0MINIMIZATION", listfile)
  if(!is.na(minim)) {
    fin.minim <- pmatch("1", listfile[minim:length(listfile)], 
			duplicates.ok = T)
    termes <- listfile[minim:(minim + fin.minim - 2)]
    termes <- substring(termes, 2)
  } else {
    termes <- NULL
    #ret.list <- list(term = termes, 
    #                 ofv = NULL, 
    #                 thetas = NULL, 
    #                 omega = NULL, 
    #                 sigma = NULL, 
    #                 sethetas = NULL, 
    #                 seomegas = NULL, 
    #                 sesigmas = NULL)
    #return(ret.list)
  }


  ## Find ofv
  if(is.null(version$language)){
    cat("need to use R for this version of Xpose")
    ##&&
    ## platform() == "WIN386" &&
    ## version$major <6) {
    ##minvalpat <- "*MINIMUM*VALUE*"
  } else {
    minvalpat <- "MINIMUM VALUE"
  }
  line <- grep(minvalpat, listfile)
  ofvline <- listfile[line + 4] ## Need +3 for S-PLUS?
  ofv <- as.numeric(substring(ofvline, 52, 71))	
  

  ## Find parameter estimates
  if(is.null(version$language)) {
    cat("need to use R for this version of Xpose")
    ## &&
    ##  platform() == "WIN386" &&
    ##  version$major < 6) {
    ## finalparpat <- "*FINAL*PARAMETER*"
    ## sepat <- "*STANDARD*ERROR*OF"
    ## tmatpat <- "*T MATRIX*"
    ## thvecpat <- "*THETA*"
    ## omegapat <- "*OMEGA*"
    ## sigmapat <- "*SIGMA*"
    ## pluspat <- "*+*"
    ## etpat <- "*ET*"
    ## eppat <- "*EP*"
    ## covmatpat <- "*COVARIANCE*MATRIX*OF*ESTIMATE*"
    ## tablepat <- "*TABLES*OF*DATA*"
    ## notepat   <- "*1 Note:*"
  } else {
    finalparpat <- "FINAL PARAMETER"
    sepat <- "STANDARD ERROR OF"
    tmatpat <- "\\*+ +T MATRIX +\\*+"
    rmatpat <- "\\*+ +R MATRIX +\\*+"
    smatpat <- "\\*+ +S MATRIX +\\*+"
    thvecpat <- "THETA"
    omegapat <- "OMEGA"
    sigmapat <- "SIGMA"
    pluspat <- "\\+"
    etpat <- "ET"
    eppat <- "EP"
    covmatpat <- "COVARIANCE MATRIX OF ESTIMATE"
    tablepat <- "TABLES OF DATA"
    notepat   <- "1 Note" # Fix for c255
  }
  finline  <- grep(finalparpat, listfile)
  seline   <- grep(sepat, listfile)
  tmatline <- grep(tmatpat, listfile)
  rmatline <- grep(rmatpat, listfile)
  smatline <- grep(smatpat, listfile)
  noteline <- grep(notepat, listfile)
  tableline <- grep(tablepat, listfile)

  if(length(seline) == 0 && length(tmatline) == 0
     && length(noteline) == 0 && length(tableline) == 0
     && length(rmatline) == 0
     && length(smatline) == 0) {
    
     if(length(grep(pluspat, listfile[length(listfile)])) == 0) {
      final.block <- listfile[finline:(length(listfile) - 1)]
    } else  {
      final.block <- listfile[finline:length(listfile)]
    }
  } else if(length(seline) !=0) {
    final.block <- listfile[finline:seline[1]]
  } else if (length(noteline)!=0) {
    ## If the last line of the lst file does not include a line
    ## beginning with a plus, i.e. an omega or sigma estimate
    ## This should always be true if length(noteline) >0
    if(length(grep(pluspat, listfile[length(listfile)])) == 0) {
      g <- 1
      final.block <- listfile[finline:(length(listfile) - (g+1))]
      ## This is tricky. The while loop is dangerous.
      while(length(grep(pluspat, listfile[length(listfile)-g])) == 0) {
        final.block <- listfile[finline:(length(listfile) - (g+1))]
        g <- g+1
      }
    }
  } else if (length(tmatline)!=0){
    final.block <- listfile[finline:(tmatline-3)]
  } else if (length(rmatline)!=0){
    final.block <- listfile[finline:(rmatline-3)]
  } else if (length(smatline)!=0){
    final.block <- listfile[finline:(rmatline-3)]
  } else if (length(tableline)!=0){
    final.block <- listfile[finline:(tableline-3)]
  } else {
    stop("the NONMEM output file has a strange format and cannot be read")
  }

  ## Check if we have sigmas. If not set sigmaline to length(final.block)
  sigmaline <- grep(sigmapat, final.block)
  nosigma <- 0
  if(length(sigmaline) == 0) {
    nosigma <- 1
    sigmaline <- length(final.block)
  }
  
  ## Find thetas
  nthlines <- grep(omegapat, final.block) - 4 - 1
  nthlines <- nthlines/2
  thetas <- NULL
  for(i in (4 + 1 + nthlines):(grep(omegapat, final.block) - 1))
    thetas <- paste(thetas, final.block[i], sep = " ")
  thetas <- string2num(thetas)

  ## Find omegas
  omega.block <- final.block[(grep(omegapat, final.block) + 1):
                             (sigmaline - 1)]

  omega.block <- omega.block[ - grep(etpat, omega.block)]
  omegas <- substring(omega.block, 2)
  starlines <- grep("\\*\\*\\*\\*",omegas)
  if(length(starlines)!=0){
    omegas <- omegas[-starlines]
  }
  omegas <- omegas[sapply(omegas, nchar) != 0]
  omega <- list()
  for(i in 1:length(omegas))
    omega[[i]] <- string2num(omegas[i])
  omega <- fix.wrapped.lines(omega)
 
  ## Find sigmas
  if(!nosigma) {
    if(length(seline) == 0) {
      sigma.block <- final.block[(grep(sigmapat, final.block) + 1):
                                 length(final.block)]
    } else {
      sigma.block <- final.block[(grep(sigmapat, final.block) + 1):
                                 (length(final.block) - 4)]
    }

    ## check to make sure that there is no extra text at end of block
    pluslines <- grep(pluspat, sigma.block) # find the lines with '+' at the start
    lastplusline <- pluslines[length(pluslines)]  # last line with '+' at the start
    nextline <- lastplusline+1
    while (((nextline+1) < length(sigma.block)) &&
           length(grep("[[:alnum:]]", sigma.block[nextline]))!=0 ) {
      nextline <- nextline+1
    }
    lastSigmaLine <- nextline-1
    sigma.block <- sigma.block[1:lastSigmaLine]

        
    ## now extract sigmas
    sigma.block <- sigma.block[ - grep(eppat, sigma.block)]
    sigmas <- substring(sigma.block, 2)
    sigmas <- sigmas[sapply(sigmas, nchar) != 0]
    sigma <- list()
    for(i in 1:length(sigmas))
      sigma[[i]] <- string2num(sigmas[i])
    sigma <- fix.wrapped.lines(sigma)
  } else {
    sigma <- NULL
  }
  ##
  ## Find standard errors
  ##
  if(length(seline) == 0) {
    sethetas <- NULL
    seomega <- NULL
    sesigma <- NULL
  } else {
    covmatline <- grep(covmatpat, listfile)[1]
    se.block <- listfile[seline:(covmatline - 4)]
    
    sigmaline <- grep(sigmapat, se.block)
    nosigma <- 0
    if(length(sigmaline) == 0) {
      nosigma <- 1
      sigmaline <- length(se.block)
    }
  
    ## Find sethetas
    nthlines <- grep(omegapat, se.block) - 4 - 1
    nthlines <- nthlines/2
    sethetas <- NULL
    for(i in (4 + 1 + nthlines):(grep(omegapat, se.block) - 1))
      sethetas <- paste(sethetas, se.block[i], sep = " ")
    sethetas <- string2num(sethetas)
    na2zero <- function(x)
      {
        if(is.na(x))
          return(0)
        else return(x)
      }
    ## Find omegas
    omega.block <- se.block[(grep(omegapat, se.block) + 1):
                            (sigmaline - 1)]
    omega.block <- omega.block[ - grep(etpat, omega.block)]
    seomegas <- substring(omega.block, 2)
    seomegas <- seomegas[sapply(seomegas, nchar) != 0]
    seomega <- list()
    for(i in 1:length(seomegas)) {
      ##seomega[[i]] <- sapply(string2num(seomegas[i]), na2zero)
      seomega[[i]] <- string2num(seomegas[i])
    }
    seomega <- fix.wrapped.lines(seomega)
  

    ## Find sigmas
    if(!nosigma) {
      sigma.block <- se.block[(sigmaline + 1):
                              length(se.block)]

      sigma.block <- sigma.block[ - grep(eppat, sigma.block)]
      sesigmas <- substring(sigma.block, 2)
      sesigmas <- sesigmas[sapply(sesigmas, nchar) != 0]
      sesigma <- list()
      for(i in 1:length(sesigmas))
        sesigma[[i]] <- string2num(sesigmas[i])
      sesigma <- fix.wrapped.lines(sesigma)
    } else {
      sesigma <- NULL
    }

  }
  
  ret.list <- list(term = termes, ofv = ofv, thetas = thetas, omega = 
                   omega, sigma = sigma, sethetas = sethetas,
                   seomegas = seomega, 
                   sesigmas = sesigma)
  return(ret.list)
}

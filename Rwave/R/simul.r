#########################################################################
#
#               (c) Copyright  1997                             
#                          by                                   
#      Author: Rene Carmona, Bruno Torresani, Wen-Liang Hwang   
#                  Princeton University 
#                  All right reserved                           
#########################################################################

#*********************************************************************#
# mnpval
#*********************************************************************#

mnpval <- function(inputdata, maxresoln, wl=128, scale=FALSE)
{
  x <- adjust.length(inputdata)
  s <- x$signal
  np <- x$length
  num.of.windows <- (np/wl - 1) * 4 + 1

  pval <- matrix(0, nrow=maxresoln, ncol=np)
  pval <- t(pval)
  dim(pval) <- c(length(pval), 1)

  z <- .C("normal_pval_compute",
          a=as.double(pval),
          as.double(s),
          as.integer(maxresoln),
          as.integer(np),
          as.integer(num.of.windows),
          as.integer(wl),
           PACKAGE="Rwave")

  pval <- t(z$a)
  dim(pval) <- c(np, maxresoln)

  plotResult(pval, s, maxresoln, scale)
  list( original=s, pval=pval, maxresoln=maxresoln, np=np )
}

#*********************************************************************#
# mbpval
#*********************************************************************#

mbpval <- function(inputdata, maxresoln, wl=128, scale=FALSE)
{
  x <- adjust.length(inputdata)
  s <- x$signal
  np <- x$length
  num.of.windows <- (np/wl - 1) * 4 + 1

  pval <- matrix(0, nrow=maxresoln, ncol=np)
  pval <- t(pval)
  dim(pval) <- c(length(pval), 1)

  z <- .C("compute_mallat_bootstrap_pval",
	a=as.double(pval),
	as.double(s),
	as.integer(maxresoln),
	as.integer(np),
	as.integer(num.of.windows),
	as.integer(wl),
           PACKAGE="Rwave")

  pval <- t(z$a)
  dim(pval) <- c(np, maxresoln)

  plotResult(pval, s, maxresoln, scale)
  list( original=s, pval=pval, maxresoln=maxresoln, np=np )
}

#*********************************************************************#
# mntrim
#*********************************************************************#

mntrim <- function(extrema, scale=FALSE, prct=.95)
{
  s <- extrema$original
  maxresoln <- extrema$maxresoln
  np <- extrema$np
  sample.size <- 128
  nthresh <- 1:maxresoln

  # Find the threshold for each resoln

  z <- .C("nthresh_compute",
	a=as.double(nthresh),
	as.double(s),
	as.integer(maxresoln),
	as.integer(sample.size),
	as.double(prct),
           PACKAGE="Rwave")

  nthresh <- z$a

  trim <- matrix(0, nrow=np, ncol=maxresoln)

  for (j in 1:maxresoln)
  {
    # Keep the extrema if the absolute value of the extrema >= threshold
    temp <- (abs(extrema$extrema[,j]) >= nthresh[j])
    trim[,j] <- temp * extrema$extrema[,j]
  }

  cat("number of extrema left =", sum(trim!=0), "\n")

  plotResult(trim, s, maxresoln, scale)
  list( original=s, extrema=trim, Sf=extrema$Sf, maxresoln=maxresoln, np=np )
}

#*********************************************************************#
# mbtrim
#*********************************************************************#

mbtrim <- function(extrema, scale=FALSE, prct=.95)
{
  s <- extrema$original
  maxresoln <- extrema$maxresoln
  np <- extrema$np
  sample.size <- 128
  bthresh <- 1:maxresoln

  # Find the threshold for each resoln

  z <- .C("bthresh_compute",
	a=as.double(bthresh),
	as.double(s),
	as.integer(maxresoln),
	as.integer(sample.size),
	as.double(prct),
           PACKAGE="Rwave")

  bthresh <- z$a

  trim <- matrix(0, nrow=np, ncol=maxresoln)
  for (j in 1:maxresoln)
  {
    # Keep the extrema if the absolute value of the extrema >= threshold
    temp <- (abs(extrema$extrema[,j]) >= bthresh[j])
    trim[,j] <- temp * extrema$extrema[,j]
  }

  cat("number of extrema left =", sum(trim!=0), "\n")

  plotResult(trim, s, maxresoln, scale)
  list( original=s, extrema=trim, Sf=extrema$Sf, maxresoln=maxresoln, np=np )
}











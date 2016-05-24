######################################################################
# neqtl.R
#
# Karl W Broman
# Ported from http://github.com/kbroman/neqtl on 27 apr 2012
#
#     This program is free software; you can redistribute it and/or
#     modify it under the terms of the GNU General Public License,
#     version 3, as published by the Free Software Foundation.
# 
#     This program is distributed in the hope that it will be useful,
#     but without any warranty; without even the implied warranty of
#     merchantability or fitness for a particular purpose.  See the GNU
#     General Public License, version 3, for more details.
# 
#     A copy of the GNU General Public License, version 3, is available
#     at http://www.r-project.org/Licenses/GPL-3
#
# Contains: neqtl, smoothall, smoothchr
######################################################################

neqtl <- function(sigpos.out,chr,pos,win=5)
     smoothall(sigpos.out,chr,pos,window=win)

smoothall <- function(themax, thechr, thepos, window=5)
{
  thesmooth <- vector("list", length(themax))
  names(thesmooth) <- names(themax)
  for(i in names(themax))
      thesmooth[[i]] <- smoothchr(themax[[i]], thepos[thechr==i], window=window)
   out <- NULL
  for(i in 1:length(thesmooth))
    out <- rbind(out, data.frame(chr=rep(names(themax)[i], nrow(thesmooth[[i]])),
                    pos=thesmooth[[i]][,1], nqtl=thesmooth[[i]][,2]))
  class(out) <- c("scanone", "data.frame")

  ## This chokes right here!

  rownames(out) <- paste("c", out[,1], ".loc", 1:nrow(out), sep="")

  out
}

## Uses positions from thepos for smoothing: ATB 9/10/09 ##
## In theloc by=0.2 was outside the seq() function--moved it inside  ATB 12/15/09 ##
smoothchr <- function(themax, thepos, window=5)
{
  ## theloc <- sort(unique(c(thepos, seq(0, max(thepos), by=0.2))))
  theloc <- thepos

  temploc <- c(themax, theloc)
  tempval <- c(rep(1, length(themax)), rep(0, length(theloc)))
  o <- order(temploc)
  temploc <- temploc[o]
  tempval <- tempval[o]
  smoothed <- runningmean(temploc, tempval, at=thepos, window=window, what="sum")
  ## NB: This differs from R/neqtl, where at=theloc.
  cbind(thepos, smoothed)
}

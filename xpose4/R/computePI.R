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

"computePI" <-
  function(x,y,object,limits=object@Prefs@Graph.prefs$PIlimits,
           logy=FALSE,logx=FALSE,onlyfirst=FALSE,
           inclZeroWRES=FALSE,PI.subset=NULL) {

    ## this prediction interval should be from data passed by the calling function
    ## should not be computed here!
    data <- SData(object,inclZeroWRES=inclZeroWRES,onlyfirst=onlyfirst,subset=PI.subset)
    if(is.null(data)) return(NULL)

    ## x is not in data
    if(!any(names(data)==x)) {
      return(NULL)
    }

    ## y is not in data
    if(!any(names(data)==y)) {
      return(NULL)
    }

    if(logy) data[,y] <- log10(data[,y])

    data <- data[order(data[,x]),]

    ## Find the number of intervals
    tids <- data[,x]
    tims <- unique(tids)
    nint <- 12
    if(length(tims) <= 12) nint <- length(tims)-1

    bins <- eq.xpose(tids,number=nint,overlap=0)

    ## Add bin indicator to data set
    data[,"bin"] <- rep(0,nrow(data))
    for(b in 1:nint) {
      #cat(b,bins$lower[b],"-",bins$upper[b],"\n")
      if(b == 1) {
        data[,"bin"] <- ifelse(tids <= bins$upper[b],b,data[,"bin"])
      } else {
        data[,"bin"] <- ifelse(tids >bins$lower[b] & tids <=bins$upper[b],b,data[,"bin"])
      }
    }


    qua  <- sapply(1:nint,function(xx,dat)
                   {quantile(dat[dat[,"bin"]==xx,y],probs=limits)},data)
    qua    <- data.frame(t(qua))
    names(qua) <- c("lower","upper")
    qua$median <- sapply(1:nint,function(xx,dat)
                         {median(dat[dat[,"bin"]==xx,y])},data)
    qua$mean   <- sapply(1:nint,function(xx,dat)
                         {mean(dat[dat[,"bin"]==xx,y])},data)

    qua$Xmiddle <- bins$middle
    qua$Xlower  <- bins$lower
    qua$Xupper  <- bins$upper
    return(qua[!is.na(qua[,1]),])

  }


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

create.strat.rand <- function(data,object,x,y,frac,dilci,seed=NULL) {

  if(!is.null(seed)) {
    set.seed(seed)
  } else {
    seed <- "noseed"
  }

  if(missing(frac))  frac  <- object@Prefs@Graph.prefs$dilfrac
  if(missing(dilci)) dilci <- object@Prefs@Graph.prefs$dilci
  
  facnam <- paste("RS",seed,sep="")

  yres <- residuals(loess(formula(paste(y,"~",x)),data))
  lims <- quantile(yres,
                   probs=c(1-dilci,dilci))

  ## Can sample from the 1s
  mp3fun <- function(x,lims) {
    ret <- 1
    if(any(x < lims[1])) ret <- 0
    if(any(x > lims[2])) ret <- 0
    ret
  }

  RR     <- tapply(yres,data[,xvardef("id",object)],FUN="mp3fun",lim=lims)
  ids    <- unique(data[,xvardef("id",object)])
  RRR    <- sample(c(0,1),size=length(ids),replace=TRUE,
                         prob=c(frac,1-frac))

  tmp    <- data.frame(ID=as.factor(names(RR)),
                       R=RR*RRR)

  names(tmp) <- c(xvardef("id",object),facnam)
  newdata       <- merge(tmp,data,by=xvardef("id",object),sort=FALSE)
  
  newdata[,facnam] <- as.factor(newdata[,facnam])
  xlabel(object)    <- c(facnam,facnam)
  data[,facnam]    <- newdata[,facnam]

  names(data)[length(names(data))] <- facnam
  invisible()
  return(data)
}

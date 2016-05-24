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

create.rand <- function(data,object,frac,seed=NULL) {

  if(!is.null(seed)) {
    set.seed(seed)
  } else {
    seed <- "noseed"
  }

  if(missing(frac)) frac <- object@Prefs@Graph.prefs$dilfrac

  facnam <- paste("R",seed,sep="")
  ids    <- unique(data[,xvardef("id",object)])
  
  tmp    <- data.frame(ID=ids,R=sample(c(0,1),size=length(ids),replace=TRUE,
                         prob=c(frac,1-frac)))  

  names(tmp)        <- c(xvardef("id",object),facnam)
  newdata           <- merge(data,tmp,by=xvardef("id",object))
  newdata[,facnam]  <- as.factor(newdata[,facnam])

  xlabel(object) <- c(facnam,facnam)
  data[,facnam] <- newdata[,facnam]
  
  names(data)[length(names(data))] <- facnam
  invisible()
  return(data)
}

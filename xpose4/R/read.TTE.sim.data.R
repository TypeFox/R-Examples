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

read.TTE.sim.data <-
  function(sim.file,
           subset=NULL,
           headers=c("REP", "ID", "DV", "TIME", "FLAG2","DOSE"),
           xpose.table.file=FALSE,
           ...)
{
  if(xpose.table.file){
    sim.data <- read.nm.tables(sim.file,quiet=T,...)
  } else {
    sim.data <- read.table(sim.file,skip=0,header=F)
    names(sim.data) <- headers
  }
  d.sim.sub <- sim.data[eval(parse(text=paste("sim.data$", subset))),]

  ## add a unique ID identifier
  rle.result <- rle(d.sim.sub$ID)
  rle.result$values <- 1:length(rle.result$values)
  new.id <- inverse.rle(rle.result)
  d.sim.sub$new.id <- new.id

  ## old way to add a unique ID identifier
  ## d.sim.sub$new.id <- max(unique(d.sim.sub$ID))*(d.sim.sub$REP-1)+d.sim.sub$ID

  tmp <- d.sim.sub
  tmp2 <- tmp[!duplicated(tmp$new.id, fromLast = TRUE),]
  d.sim <- tmp2
  return(d.sim)
}

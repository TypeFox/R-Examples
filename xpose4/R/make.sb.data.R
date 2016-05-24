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

"make.sb.data"  <- function(data,idv,dv,nbins=6,by=NULL,
                            by.nbins  = 6,
                            #ordby     = NULL,
                            #byordfun  = "median",
                            #shingnum  = 6,
                            #shingol   = 0.5,
                            ...)
{

  if(is.null(idv)){
    idv <- "all.values"
    data[idv] <- 1
  }
  if(nbins < length(unique(data[,idv]))){
    data$bins.tmp <- xpose.bin(data,idv,bins=nbins)
    idv <- "bins.tmp"
  }
  doses <- unique(data[,idv])
  #doses <- as.vector(doses[order(doses)],"numeric")
  doses <- doses[order(doses)]
  dvs   <- unique(data[,dv])
  #dvs   <- as.vector(dvs[order(dvs)],"numeric")
  dvs   <- dvs[order(dvs)]

 
  ## get conditioning
  if(is.null(by)){## No conditioning 
    nlevels <- 1
    levs <- 1
    by <- "all.values"
    data[by] <- 1
  } else {
    if(by.nbins < length(unique(data[,by]))){
      data$by.bins.tmp <- xpose.bin(data,by,bins=by.nbins)
      by <- "by.bins.tmp"
    }
    levs <- unique(data[,by])
    levs <- levs[order(levs)]
    nlevels <- length(levs)
    
    ##     ##could be done with shingles like this...kinda
    ##     if(!is.factor(data[,by])) {
    ##       data[,by] <- equal.count(data[,by],number=shingnum,overl=shingol)
    ##     } else {
    ##       if(!is.null(ordby)) {
    ##         data[,by] <- reorder.factor(data[,by],data[,ordby],byordfun)
    ##       }
    ##       if(names(data[,by,drop=F])!="ind") {
    ##         levels(data[,by]) <-
    ##           paste(xlabel(names(data[,by,drop=F]),object),":",   ## Needs to be fixed
    ##                 levels(data[,by]),sep="")
    ##       }
    ##     }
    ##     ## end shingle stuff  
  }

  ## Set up the data frame for the data to be plotted
  num.row <- length(dvs)
  num.col <- length(doses)
  ##if(!is.null(by)) num.col <- num.col+1
  ret   <- data.frame(matrix(nrow = num.row,
                             ncol = num.col))
  wdths <- rep(1,length(doses))
  #row.names(ret) <- paste("P", as.vector(dvs,"numeric"), sep = "")
  row.names(ret) <- paste(dv,"=",dvs, sep = "")
  names(ret) <- doses

  ## Set up the data frame for the data to be plotted
  num.col.new <- 5
  num.row.new <- length(doses)*length(dvs)*nlevels
  ret.new   <- data.frame(matrix(nrow = num.row.new,
                             ncol = num.col.new))
  names(ret.new) <- c("idv","dv","proportion","by.var","wdth")
  if(!is.null(levels(doses))){
    ret.new["idv"] <- factor(ret.new["idv"],levels=levels(doses))
  }
  ret.new["dv"] <- factor(ret.new["dv"],levels=levels(dvs))
  if(!is.null(levels(levs))){
    ret.new["by.var"] <- factor(ret.new["by.var"],levels=levels(levs))
  }

  ## add loop here to go through all the levels
  i <- 1
  for(LEVS in 1:nlevels){
    tmp.by=levs[LEVS]
    dat1 <- data[data[,by] == levs[LEVS], ,drop=F ]
    for(DOS in 1:length(doses)){
      tmp.idv <- doses[DOS] 
      dat2 <- dat1[dat1[,idv] == doses[DOS], ,drop=F ]
      tmp.wdth <- nrow(dat2)
      for(DV in 1:length(dvs)){
        tmp.dv <- dvs[DV]
        if(nrow(dat2) == 0) {
          tmp.proportion <- 0
        } else {
          if(is.null(nrow(dat2[dat2[, dv] == dvs[DV],,drop=F]))) {
            tmp.proportion <- 0
          } else {
            tmp.proportion <- nrow(dat2[dat2[, dv] == dvs[DV],,drop=F])/
            nrow(dat2)
          }
        }
        ret.new[i,"idv"] <- tmp.idv
        ret.new[i,"dv"] <- tmp.dv
        ret.new[i,"proportion"] <- tmp.proportion
        ret.new[i,"by.var"] <- tmp.by
        ret.new[i,"wdth"] <- tmp.wdth
        i <- i+1
      }      
    }
  }
    #ret$idv[(ii-1)*] <- rep(paste(dv,"=",dvs, sep = ""),nlevels)

   ## Fill in the data frame
  for(dos in 1:length(doses)) {
    dat1 <- data[data[,idv] == doses[dos], ,drop=F ]
    wdths[dos] <- nrow(dat1)

    for(d in 1:num.row) {

      if(nrow(dat1) == 0) {
        ret[d, dos] <- 0
        next
      }

      if(is.null(nrow(dat1[dat1[, dv] == dvs[d],,drop=F]))) {
        ret[d, dos] <- 0
      } else {
        ret[d, dos] <- nrow(dat1[dat1[, dv] == dvs[d],,drop=F])/
          nrow(dat1)
      }
    }
  }
    
  retlist <- list(ret=ret,wdths=wdths)
  retlist.new <- list(ret=ret.new)
  
  ##return(retlist)
  return(retlist.new)
}

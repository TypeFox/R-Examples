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

"continuous.table"  <- function(object, 
                                vars,
                                onlyfirst=TRUE,
                                subset=xsubset(object),
                                inclZeroWRES=FALSE,
                                miss=object@Prefs@Miss) # can be a number
{
  
  data <- Data(object,onlyfirst=onlyfirst,subset=subset,inclZeroWRES=inclZeroWRES)
  usemiss <- FALSE
  
  for (nam in vars) {
    if(any(is.na(data[[nam]]) | data[[nam]]==miss)) {
      usemiss <- TRUE
    }
  }
  
  if (usemiss == TRUE) {
    ret.mat <- matrix(0,ncol=9,nrow=1+length(vars))
    ret.mat[1,] <- c("","Mean","SD","Q1","Median","Q3","Range","N","Missing")
    
  } else {
    ret.mat <- matrix(0,ncol=8,nrow=1+length(vars))
    ret.mat[1,] <- c("","Mean","SD","Q1","Median","Q3","Range","N")
    
  }
  
  i <- 1
  
  for(nam in vars) {
    i <- i+1
    
    micov <- subset(data[[nam]], is.na(data[[nam]]) | data[[nam]]==miss)
    nomicov<- subset(data[[nam]], !is.na(data[[nam]]) & data[[nam]]!=miss)
    suma <- summary(nomicov)[c(4,2,3,5,6,1)]

    if (usemiss==TRUE) {
      ret.mat[i,] <- c(nam,
                       suma[1],
                       signif(sd(nomicov), digits=4),
                       suma[2:4],
                       paste(suma[6],"-",suma[5],sep=""),
                       length(nomicov),
                       paste(length(micov)," (",
                             sprintf("%.1f",
                                     100*length(micov)/(length(micov)+length(nomicov))),
                             "%)",sep="")
                       )

    } else {
      ret.mat[i,] <- c(nam,
                       suma[1],
                       signif(sd(nomicov), digits=4),
                       suma[2:4],
                       paste(suma[6],"-",suma[5],sep=""),
                       length(nomicov)
                       )

    }
  }

  class(ret.mat) <- "char.matrix"

  return(ret.mat)
}




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

SData <- function(object,
                  inclZeroWRES=FALSE,
                  onlyfirst=FALSE,
                  subset=NULL,
                  samp=NULL) {

  data <- object@SData
  if(is.null(data)) return(NULL)

  if(!is.null(samp)) {
    data <- data[data$iter==samp,]       # something is wrong with this line
    #cat(samp)
    if(is.null(data)) return(NULL)
  }
  
  if(!inclZeroWRES) {
    data <- data[data[,xvardef("wres",object)]!=0,]
  }

  if(onlyfirst) {
    id  <- xvardef("id",object)
    ind <- paste(data[,id],data$iter,sep="")
    data<- data[!duplicated(ind),]
  }
  
  if(!is.null(subset)) {
    #on.exit(detach(data))
    #attach(data)
    # data <- data[eval(parse(text=subset)),]
    #data <- data[eval(parse(text=paste("data$", subset))),] # fix subsets 22/3/06
    data<-with(data,data[eval(parse(text=subset)),])
    
    if(dim(data)[1]==0) return(NULL)
  }
  
  return(data)

}

"SData<-" <- function(object,value) {

  Snro <- dim(value)[1]
  Dnro <- dim(Data(object,inclZeroWRES=TRUE))[1]

  if(Dnro == 0) return("Data should be set before SData!")
  
  ## Check to see if the length of the SData is an even multiplier
  ## of the xData.
  if (!is.null(Snro)) {
    if(regexpr("\\.",as.character(Snro/Dnro)) !=-1) {
      cat("The length of the Data and the SData do not match!\n")
      return(object)
    }

    nams <- names(Data(object))
    for(n in nams) {
      #class(value[,n]) <- class(Data(object)[,n])
      if (is.factor(Data(object)[,n])) {
        value[,n] <- as.factor(value[,n])
      }
    }
                            
    nsim(object)      <- Snro/Dnro

    ## Check to see if WRES is all zero. This would indicate that
    ## ONLYSIM was used during the simulation in NONMEM. In this case
    ## many plots will fail. A fix to this is to replace the WRES
    ## column in SData with a column that is zero where the WRES in
    ## Data is zero and 1 otherwise.
    if(!any(value$WRES !=0)) {
      Data.wres  <- Data(object,inclZeroWRES=TRUE)[,"WRES"]
      SData.wres <- rep(Data.wres,Snro/Dnro)
      SData.wres <- ifelse(SData.wres==0,0,1)
      value$WRES <- SData.wres
    }

    ## Add a column with a number indicating each simulated data set
    value[,"iter"]   <- sort(rep(1:nsim(object),Dnro))
    object@SData      <- value
  
    return(object)
  } else {
    return(NULL)
  }
}

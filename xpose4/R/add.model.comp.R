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

"add.model.comp" <-
  function(object, 
           object.ref = NULL,
           onlyfirst = FALSE,
           inclZeroWRES = FALSE,
           subset = xsubset(object),
           main="Default",
           force.wres=FALSE,
           #ref.default = ".ref.db",
           ...) {

    if (is.null(object.ref)) {
      ref.list <- get.refrunno()
      if(exists(".ref.db")){
        object.ref <- eval(parse(text=".ref.db"))
      } else {
        return()
      }
      if(any(is.null(ref.list)))
        return()
    } 

    
    if(dim(object@Data)[1] != dim(object.ref@Data)[1]) {
      cat("The current database and the reference database do not have\n")
      cat("the same number of lines!\n")
      invisible()
      return()
    }

    if(is.null(check.vars(c("idlab","pred","ipred","iwres","idv"),
                          object,silent=FALSE))) {      
      return()
    }

    
    use.cwres=TRUE
    if(force.wres){
      use.cwres=FALSE
      if(is.null(check.vars(c("wres"),object,silent=FALSE))) return()
    } else {
      if(is.null(check.vars(c("cwres"),object,silent=TRUE))) {
        use.cwres=FALSE
        if(is.null(check.vars(c("wres"),object,silent=FALSE))) return()
      }
    }

    
    object@Data$PRED.REF <- abs(object.ref@Data[,xvardef("pred", object.ref)])
    object@Data$IPRED.REF <- abs(object.ref@Data[,xvardef("ipred", object.ref)])
    #object@Data$WRES.REF <- abs(object.ref@Data[,xvardef("wres", object.ref)])
    object@Data$IWRES.REF <- abs(object.ref@Data[,xvardef("iwres", object.ref)])
    #object@Data$dWRES <- abs(object@Data[,xvardef("wres", object)] - object.ref@Data[,xvardef("wres", object.ref)])
    object@Data$dIWRES <- abs(object@Data[,xvardef("iwres", object)] - object.ref@Data[,xvardef("iwres", object.ref)])

    if(use.cwres){
      object@Data$CWRES.REF <- abs(object.ref@Data[,xvardef("cwres", object.ref)])
      object@Data$dCWRES <- abs(object@Data[,xvardef("cwres", object)] - object.ref@Data[,xvardef("cwres", object.ref)])
      object@Data[,xvardef("cwres", object)] <- abs(object@Data[,xvardef("cwres", object)])
    } else {
      object@Data$WRES.REF <- abs(object.ref@Data[,xvardef("wres", object.ref)])
      object@Data$dWRES <- abs(object@Data[,xvardef("wres", object)] - object.ref@Data[,xvardef("wres", object.ref)])
      object@Data[,xvardef("wres", object)] <- abs(object@Data[,xvardef("wres", object)])
    }
    
    object@Data[,xvardef("pred", object)] <- abs(object@Data[,xvardef("pred", object)])
    object@Data[,xvardef("ipred", object)] <- abs(object@Data[,xvardef("ipred", object)])
    #object@Data[,xvardef("wres", object)] <- abs(object@Data[,xvardef("wres", object)])
    object@Data[,xvardef("iwres", object)] <- abs(object@Data[,xvardef("iwres", object)])
    

    ## |PRED| vs |PRED|
    if(!any(is.null(xvardef("pred", object))) &&
       !any(is.null(xvardef("pred", object.ref)))) {
      xlb <- paste("|",xlabel(xvardef("pred",object),object),
                   "| (Run ", object@Runno, ")",sep="")
      ylb <- paste("|",xlabel(xvardef("pred",object.ref),object.ref),
                   "| (Run ", object.ref@Runno, ")",sep="")
#      main <- paste(ylb, " vs ", xlb, sep="")

      xplot1 <- xpose.plot.default(xvardef("pred", object),
                                   "PRED.REF",
                                   object,
                                   xlb = xlb,
                                   ylb = ylb,
                                   main = NULL,
                                   abline=c(0,1),
                                   onlyfirst = onlyfirst,
                                   inclZeroWRES = inclZeroWRES,
                                   subset = subset,
                                   pass.plot.list = TRUE,
                                   ...)
    }

    ## |IPRED| vs |IPRED|
    if(!any(is.null(xvardef("ipred", object))) && !any(is.null(xvardef("ipred", object.ref)))) {
      xlb <- paste("|",xlabel(xvardef("ipred",object),object), "| (Run ", object@Runno, ")",sep="")
      ylb <- paste("|",xlabel(xvardef("ipred",object.ref),object.ref), "| (Run ", object.ref@Runno, ")",sep="")
#      main <- paste(ylb, " vs ", xlb, sep="")
      
      xplot2 <- xpose.plot.default(xvardef("ipred", object),
                                   "IPRED.REF",
                                   object,
                                   xlb = xlb,
                                   ylb = ylb,
                                   main = NULL,
                                   abline=c(0,1),
                                   onlyfirst = onlyfirst,
                                   inclZeroWRES = inclZeroWRES,
                                   subset = subset,
                                   pass.plot.list = TRUE,
                                   ...)
    }
    if(use.cwres){
                                        # |dCWRES| vs IDV
      if(!any(is.null(xvardef("cwres", object))) && !any(is.null(xvardef("cwres", object.ref)))) {
        xlb <- paste(xlabel(xvardef("idv",object),object),sep="")
        ylb <- paste("|dCWRES| (Run ", object@Runno, " - Run ",object.ref@Runno,")",sep="")
                                        #      main <- paste(ylb, " vs ", xlb, sep="")
        
        xplot3 <- xpose.plot.default(xvardef("idv", object),
                                     "dCWRES",
                                     object,
                                     xlb = xlb,
                                     ylb = ylb,
                                     main = NULL,
                                     onlyfirst = onlyfirst,
                                     inclZeroWRES = inclZeroWRES,
                                     subset = subset,
                                     pass.plot.list = TRUE,
                                     ...) 
        
      }

    } else {
                                        # |dWRES| vs IDV
      if(!any(is.null(xvardef("wres", object))) && !any(is.null(xvardef("wres", object.ref)))) {
        xlb <- paste(xlabel(xvardef("idv",object),object),sep="")
        ylb <- paste("|dWRES| (Run ", object@Runno, " - Run ",object.ref@Runno,")",sep="")
                                        #      main <- paste(ylb, " vs ", xlb, sep="")
        
        xplot3 <- xpose.plot.default(xvardef("idv", object),
                                     "dWRES",
                                     object,
                                     xlb = xlb,
                                     ylb = ylb,
                                     main = NULL,
                                     onlyfirst = onlyfirst,
                                     inclZeroWRES = inclZeroWRES,
                                     subset = subset,
                                     pass.plot.list = TRUE,
                                     ...) 
        
      }
    }
    
                                        # |dIWRES| vs IDV
    if(!any(is.null(xvardef("iwres", object))) && !any(is.null(xvardef("iwres", object.ref)))) {
      xlb <- paste(xlabel(xvardef("idv",object),object),sep="")
      ylb <- paste("|diWRES| (Run ", object@Runno, " - Run ",object.ref@Runno,")",sep="")
#      main <- paste(ylb, " vs ", xlb, sep="")
      
      xplot4 <- xpose.plot.default(xvardef("idv", object),
                                   "dIWRES",
                                   object,
                                   xlb = xlb,
                                   ylb = ylb,
                                   main = NULL,
                                   onlyfirst = onlyfirst,
                                   inclZeroWRES = inclZeroWRES,
                                   subset = subset,
                                   pass.plot.list = TRUE,
                                   ...)
    }

    ## create enpty list for plots
    num.of.plots <- 4
    plotList <- vector("list",num.of.plots)

    plotList[[1]] <- xplot1
    plotList[[2]] <- xplot2
    plotList[[3]] <- xplot3
    plotList[[4]] <- xplot4      

        
    
    default.plot.title <- "Additional model comparison plots"
    plotTitle <- xpose.multiple.plot.title(object=object,
                                           plot.text = default.plot.title,
                                           subset=subset,
                                           main=main,
                                           ...)
  
    
    obj <- xpose.multiple.plot(plotList,plotTitle,...)
    return(obj)

  }


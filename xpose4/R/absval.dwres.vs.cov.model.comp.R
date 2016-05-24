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

"absval.dwres.vs.cov.model.comp" <-
  function(object, 
           object.ref = NULL,
           type = NULL,
           ylb=expression(paste("|", Delta, "WRES|")),
           main="Default",
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
    
                                        #ref.db <- ref.list$ref.db
                                        #ref.runno <- ref.list$ref.runno
    
    if(dim(object@Data)[1] != dim(object.ref@Data)[1]) {
      cat("The current database and the reference database do not have\n")
      cat("the same number of lines!\n")
      return()
    }

    if ((is.null(xvardef("idlab",object))) || (is.null(xvardef("wres",object)))) {
      cat("The required variables (ID label, WRES) are not set in the database!\n")
      return()
    } 
    
    if(any(is.null(xvardef("covariates",object)))) {
      return(cat("No covariates found in the current database!\n"))
    }

    object@Data$dWRES <- abs(object@Data[,xvardef("wres", object)] -
                             object.ref@Data[,xvardef("wres", object.ref)])
    ##object@Data[,xvardef("pred", object)] <- abs(object@Data[,xvardef("pred", object)])

    ## create list for plots
    number.of.plots <- 0
    for (i in xvardef("covariates", object)) {
      number.of.plots <- number.of.plots + 1
    }
    plotList <- vector("list",number.of.plots)
    plot.num <- 0 # initialize plot number
    
    for (i in xvardef("covariates", object)) {
      
            
      if (is.null(type)) {
        if (!is.factor(object@Data[,i])) {
          type <- "p"
        } else {
          type = object@Prefs@Graph.prefs$type
        }
      }
      
      xplot <- xpose.plot.default(i,
                                  "dWRES",
                                  object,
                                  #xlb = xlb,
                                  ylb = ylb,
                                  #main = main,
                                  type = type,
                                  pass.plot.list=TRUE,,
                                  main=NULL,
                                  ...)

      plot.num <- plot.num+1
      plotList[[plot.num]] <- xplot
    }
    
    
    ## |dPRED| vs covariates
    default.plot.title <- paste("|WRES_(Run", object@Runno,
                                ") - WRES_(Run",object.ref@Runno,
                                ")| vs. Covariates",sep="")

    plotTitle <- xpose.multiple.plot.title(object=object,
                                           plot.text = default.plot.title,
                                           no.runno=TRUE,
                                           main=main,
                                           ...)
    
    xpose.multiple.plot.default(plotList,plotTitle=default.plot.title,...)

    invisible()
}


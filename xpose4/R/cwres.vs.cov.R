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

"cwres.vs.cov" <-
  function(object,
                                        #xlb  = NULL,
           ylb  = "CWRES",
                                        #onlyfirst=FALSE,
                                        #inclZeroWRES=FALSE,
                                        #subset=xsubset(object),
                                        # abline=c(0,1),
           smooth=TRUE,
                                        #abllwd=2,
           type="p",
                                        #mirror=FALSE,
                                        #seed  = NULL,
                                        #prompt = TRUE,
           main="Default",
           ...) {

      ## check for arguments in function
      if(is.null(check.vars(c("covariates","cwres"),
                            object,silent=FALSE))) {
          return()
      }

      ## create list for plots
      number.of.plots <- 0
      for (i in xvardef("covariates", object)) {
          number.of.plots <- number.of.plots + 1
      }
      plotList <- vector("list",number.of.plots)
      plot.num <- 0 # initialize plot number

      ## loop (covs)
      for (j in xvardef("covariates", object)) {

          xplot <- xpose.plot.default(j,
                                      xvardef("cwres",object),
                                      object,
                                      main=NULL,
                                        #xlb = xlb,
                                      ylb = ylb,
                                        #abline=abline,
                                        #abllwd=abllwd,
                                      smooth=smooth,
                                      type=type,
                                        #subset=subset,
                                      pass.plot.list=TRUE,
                                      ...)
          plot.num <- plot.num+1
          plotList[[plot.num]] <- xplot
      }


      default.plot.title <- paste(xlabel(xvardef("cwres",object),object),
                                  " vs ",
                                  "Covariates",
                                  sep="")
      plotTitle <- xpose.multiple.plot.title(object=object,
                                             plot.text = default.plot.title,
                                             main=main,
                                             ...)
      obj <- xpose.multiple.plot(plotList,plotTitle,...)
      return(obj)

  }


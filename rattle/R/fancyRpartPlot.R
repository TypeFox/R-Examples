# Rattle: A GUI for Data Mining in R
#
# Time-stamp: <2015-07-26 11:46:14 gjw>
#
# Copyright (c) 2009-2014 Togaware Pty Ltd
#
#' Plot rpart decision trees nicely.
#'
#' @param model an rpart object
#' @param main title for the plot
#' @param sub sub title for the plot (default is a Rattle string with
#' date, time and username)
#' @param palettes a list of sequential palettes names as supported by
#' RColorBrewer::brewer.pal including Blues BuGn BuPu
#' GnBu Greens Greys Oranges OrRd PuBu PuBuGn PuRd Purples RdPu Reds
#' YlGn YlGnBu YlOrBr YlOrRd.
#' @param ... additional arguments passed on to rpart.plot::prp
#
# This files is part of Rattle.
#
# Rattle is free software: you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#
# Rattle is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with Rattle. If not, see <http://www.gnu.org/licenses/>.

fancyRpartPlot <- function(model,
                           main="",
                           sub,
                           palettes,
                           ...)
{
  if (missing(sub))
    sub <- paste("Rattle",
                 format(Sys.time(), "%Y-%b-%d %H:%M:%S"), 
                 Sys.info()["user"])
  
  num.classes <- length(attr(model, "ylevels"))

  # Generate a colour palette, with a range of 5 (palsize) colours for
  # each of the 6 (numpals) palettes. The palette is collapsed into
  # one list. We index it according to the class. Keep to the lighter
  # end of the palette to ensure printing is okay otherwise the black
  # text is hard to read.

  default.palettes <- c("Greens", "Blues", "Oranges", "Purples", "Reds", "Greys")
  if (missing(palettes))
    palettes <- default.palettes
  missed <- setdiff(1:6, seq(length(palettes)))
  palettes <- c(palettes, default.palettes[missed])

  numpals <- 6
  palsize <- 5
  pals <- c(RColorBrewer::brewer.pal(9, palettes[1])[1:5],
            RColorBrewer::brewer.pal(9, palettes[2])[1:5],
            RColorBrewer::brewer.pal(9, palettes[3])[1:5],
            RColorBrewer::brewer.pal(9, palettes[4])[1:5],
            RColorBrewer::brewer.pal(9, palettes[5])[1:5],
            RColorBrewer::brewer.pal(9, palettes[6])[1:5])
  
  # Extract the scores/percentages for each of the nodes for the
  # majority decision.  The decisions are in column 1 of yval2 and the
  # percentages are in the final num.classes columns.

  # 121106 Need to handle regression as pointed out by Yana
  # Kane-Esrig, 26 October 2012.

  if (model$method == "class")
  {
    yval2per <- -(1:num.classes)-1
    per <- apply(model$frame$yval2[,yval2per], 1, function(x) x[1+x[1]])
  }
  else
  {
    # 130329 This is the deviance relative the the total deviance measured at
    # the root node. We use this to colour the strength of the node -
    # so more intense colour means less relative deviance.
    
    #per <- 1 - (model$frame$dev/model$frame$dev[1])

    # 130329 Perhaps instead we want to use the yval as the intensity
    # of the predicted value. Currently not handling negative values.

    per <- model$frame$yval/max(model$frame$yval)
    
  }
  
  # The conversion of a tree in CORElearn to an rpart tree results in these
  # being character, so ensure we have numerics.
  
  per <- as.numeric(per)
  
  # Calculate an index into the combined colour sequence. Once we go
  # above numpals * palsize (30) start over.

  if (model$method == "class")
    col.index <- ((palsize*(model$frame$yval-1) +
                   trunc(pmin(1 + (per * palsize), palsize))) %%
                  (numpals * palsize))
  else
    col.index <- round(per * (palsize-1)) + 1

  # Ensure the index is positive. Thanks to John Vorwald, 8 Dec
  # 2014. The bug can arise when model$frame$yval are all
  # negative. The error is:
  #
  #  fancyRpartPlot(rtreeFit,main=paste('RPART:',cName))
  #  Error in pals[col.index] : only 0's may be mixed with negative subscripts
  
  col.index <- abs(col.index)

  # Determine the amount of extra information added to the nodes.

  if (model$method == "class")
    extra <- 104
  else
    extra <- 101
  
  # Generate the plot and title.
 
  rpart.plot::prp(model, type=2, extra=extra,
                  box.col=pals[col.index],
                  nn=TRUE,
                  varlen=0, faclen=0,
                  shadow.col="grey",
                  fallen.leaves=TRUE,
                  branch.lty=3, ...)
  
  title(main=main, sub=sub)
}


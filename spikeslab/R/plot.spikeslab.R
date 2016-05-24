####**********************************************************************
####**********************************************************************
####
####  SPIKE AND SLAB 1.1.5
####
####  Copyright 2013, Cleveland Clinic Foundation
####
####  This program is free software; you can redistribute it and/or
####  modify it under the terms of the GNU General Public License
####  as published by the Free Software Foundation; either version 2
####  of the License, or (at your option) any later version.
####
####  This program is distributed in the hope that it will be useful,
####  but WITHOUT ANY WARRANTY; without even the implied warranty of
####  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
####  GNU General Public License for more details.
####
####  You should have received a copy of the GNU General Public
####  License along with this program; if not, write to the Free
####  Software Foundation, Inc., 51 Franklin Street, Fifth Floor,
####  Boston, MA  02110-1301, USA.
####
####  ----------------------------------------------------------------
####  Written and Developed by:
####  ----------------------------------------------------------------
####    Hemant Ishwaran, Ph.D.
####    Director of Statistical Methodology
####    Professor, Division of Biostatistics
####    Clinical Research Building, Room 1058
####    1120 NW 14th Street
####    University of Miami, Miami FL 33136
####
####    email:  hemant.ishwaran@gmail.com
####    URL:    http://web.ccs.miami.edu/~hishwaran
####    --------------------------------------------------------------
####    J. Sunil Rao, Ph. D.
####    Professor and Director of the Division of Biostatistics, 
####    Department of Epidemiology & Public Health
####    Clinical Research Bldg, R-669
####    1120 NW 14th Street, Room 1056
####    Miami, FL 33136
####    email:  rao.jsunil@gmail.com
####    URL:    http://biostat.med.miami.edu/people/primary-faculty/sunil-rao
####  ----------------------------------------------------------------
####  Maintained by:
####    Udaya B. Kogalur, Ph.D.
####    Adjunct Staff
####    Dept of Quantitative Health Sciences
####    Cleveland Clinic Foundation
####    
####    Kogalur & Company, Inc.
####    5425 Nestleway Drive, Suite L1
####    Clemmons, NC 27012
####
####    email:  ubk@kogalur.com
####    URL:    http://www.kogalur.com
####    --------------------------------------------------------------
####
####**********************************************************************
####**********************************************************************

plot.spikeslab <- function(x, plot.type = c("path", "cv"), breaks = FALSE, ...)
{

  ### Check that object is compatible
  if (!inherits(x, "spikeslab"))
     stop("This function only works for objects of class `spikeslab'")
  
  ###check whether object inherits mixing type
  ###make suitable alterations to merge mixing results
  if (sum(inherits(x, c("spikeslab", "mixing"), TRUE) == c(1, 2)) == 2) {
     x <- x$spikeslab.obj
  }

  ### determine the plot type
  plot.type <- match.arg(plot.type)

  ### plots
  if (plot.type == "path") {
    # plot the gnet solution path (calls plot.lars)
    x$gnet.obj$type <- "gnet"
    plot.lars(x$gnet.obj, breaks = breaks, ...)
  }
  else if (plot.type == "cv" & sum(inherits(x, c("spikeslab", "cv"), TRUE) == c(1, 2)) == 2) {
    # plot the cv curve
    cv <- x$cv
    cv.plot.path <- x$cv.path
    model.size <- unlist(x$gnet.path$model.size)
    K <- length(cv)
    p <- nrow(cv.plot.path) - 1
    cv.plot.mean <- apply(cv.plot.path, 1, mean, na.rm = TRUE)
    cv.plot.se <- apply(cv.plot.path, 1, SD)/sqrt(K)
    matplot(0:p, cv.plot.path, type = c("l", "n")[1 + 1 * (K > 20)], lty = 3, col = "gray", 
         xlim = range(c(0, model.size), na.rm = TRUE),
         ylim = range(c(cv.plot.path, cv.plot.mean + cv.plot.se, cv.plot.mean - cv.plot.se), na.rm = TRUE),
         xlab="Model Size", ylab="Cross-Validated MSE")
    lines(0:p, cv.plot.mean, lty = 1, lwd = 2, col = 4)
    error.bars(0:p, cv.plot.mean + cv.plot.se, cv.plot.mean - cv.plot.se, width = 0.0025, col = 2)
  }

}

grtspts <- function(src.frame="shapefile", shapefilename=NULL, ptsframe,
   samplesize=100, SiteBegin=1, shift.grid=TRUE, do.sample=TRUE, startlev=NULL,
   maxlev=11) {

################################################################################
# Function: grtspts.r
# Purpose: Select a generalized random-tesselation stratified (GRTS) sample of a
#    finite resource
# Programmers: Tony Olsen, Tom Kincaid, Don Stevens, Christian Platt,
#   			Denis White, Richard Remington
# Date: October 8, 2002
# Last Revised: January 27, 2015
# Description:
#   This function select a GRTS sample of a finite resource.  This function uses
#   hierarchical randomization to ensure that the sample will include no more
#   than one point per cell and then picks a point in selected cells.  
# Arguments:
#   src.frame = source of the frame, which equals "shapefile" if the frame is to
#     be read from a shapefile, or "att.frame" if the frame is included in
#     ptsframe.  The default is "shapefile".
#   shapefilename = name of the input shapefile. If src.frame equals "shapefile"
#     and shapefilename equals NULL, then the shapefile or shapefiles in the
#     working directory are used.  The default is NULL.
#   ptsframe = a data frame containing id, x, y, mdcaty, and mdm.
#   samplesize = number of points to select in the sample.  The default is 100.
#   SiteBegin = first number to start siteID numbering.  The default is 1.
#   shift.grid = the option to randomly shift the hierarchical grid.  The
#     default is TRUE.
#   do.sample = option to select a sample, where TRUE means select a sample and
#     FALSE means return the entire sample frame in reverse hierarchical order.
#     The default is TRUE.  
#   startlev = initial number of hierarchical levels to use for the GRTS grid,
#     which must be less than or equal to maxlev (if maxlev is specified) and
#     cannot be greater than 11.  The default is NULL.
#   maxlev = maximum number of hierarchical levels to use for the GRTS grid,
#     which cannot be greater than 11.  The default is 11.
# Results: 
#   A data frame of sample points containing: siteID, id, x, y, mdcaty,
#   and weight.
# Other Functions Required:
#   numLevels - C function to determine the number of levels for hierarchical
#     randomization
#   cell.wt - calculates total inclusion probability for a cell
#   constructAddr - C function to construct the hierarchical address for all
#     points
#   ranho - C function to construct the randomized hierarchical address for all
#     points
#   pickGridCells - C function to select grid cells that get a sample point
#   selectpts - pick sample point(s) from selected cells
################################################################################

# If src.frame is not "shapefile", determine the minimum and maximum values for
# the grid and extent of the grid

   if(src.frame != "shapefile") {
      rx <- range (ptsframe$x)
      ry <- range (ptsframe$y)
      grid.extent <- max(rx[2] - rx[1], ry[2] - ry[1])
      temp <- 0.04*grid.extent
      grid.xmin <- rx[1] - temp
      grid.ymin <- ry[1] - temp
      grid.extent <- 1.08*grid.extent
      grid.xmax <- grid.xmin + grid.extent
      grid.ymax <- grid.ymin + grid.extent
   }

# Determine the number of levels for hierarchical randomization

   if(src.frame == "shapefile") {
      temp <- .Call("numLevels", shapefilename, samplesize, shift.grid,
         startlev, maxlev, ptsframe$id, ptsframe$mdm)
      if(is.null(temp[[1]]))
         stop("\nAn error occured while determining the number of levels for hierarchical \nrandomization.") 
      nlev <- temp$nlev
      dx <- temp$dx
      dy <- temp$dy
      xc <- temp$xc
      yc <- temp$yc
      cel.wt <- temp$cel.wt
      sint <- temp$sint
   } else {
      if(is.null(startlev)) {
         nlev <- ceiling(logb(samplesize, 4))
         if(nlev == 0)
            nlev <- 1
      } else {
         nlev <- startlev
      }
      cel.wt <- 99999
      celmax.ind <- 0
      sint <- 1
      if(shift.grid) {
         roff.x <- runif(1, 0, 1)
         roff.y <- runif(1, 0, 1)
      }
      while (any(cel.wt/sint > 1) && celmax.ind < 2 && nlev <= maxlev) {
         cat( "Current number of levels:", nlev, "\n");
         celmax <- max(cel.wt)
         nlv2 <- 2^nlev
         dx <- dy <- grid.extent/nlv2
         xc <- seq(grid.xmin, grid.xmax, length=nlv2+1)
         yc <- seq(grid.ymin, grid.ymax, length=nlv2+1)
         if(shift.grid) {
            xc <- rep(xc, nlv2+1) + (roff.x * dx)
            yc <- rep(yc, rep(nlv2+1, nlv2+1)) + (roff.y * dy)
         } else {
            xc <- rep(xc, nlv2+1)
            yc <- rep(yc, rep(nlv2+1, nlv2+1))
         }
         
# Determine total inclusion probability for each grid cell and, as necessary, 
# adjust the indicator for whether maximum of the total inclusion probabilities
# is changing   

         cel.wt <- sapply(1:length(xc), cell.wt, xc, yc, dx, dy, ptsframe)
         if(max(cel.wt) == celmax) {
            celmax.ind <- celmax.ind + 1
    	       if(celmax.ind == 2)
    	          warning("\nSince the maximum value of total inclusion probability for the grid cells was \nnot changing, the algorithm for determining the number of levels for \nhierarchical randomization was terminated.\n")
         }

# Adjust sampling interval and number of hierarchical levels

         sint <- sum(cel.wt)/samplesize
         ifelse(nlev == maxlev,
            nlev <- nlev + 1,
            nlev <- nlev + max(1, ceiling(logb(cel.wt[cel.wt > 0]/sint, 4))))
      }

#  Print the final number of levels

      cat( "Final number of levels:", nlev-1, "\n");
   }

# Assign the final number of levels

   endlev <- nlev - 1

# Remove cells with zero weight

   indx <- cel.wt > 0
   xc <- xc[indx]
   yc <- yc[indx]
   cel.wt <- cel.wt[indx]

# Construct the hierarchical address for all cells

   hadr <- .Call("constructAddr", xc, yc, dx, dy, as.integer(nlev))

# Construct randomized hierarchical addresses

   ranhadr <- .C("ranho", hadr, as.integer(length(hadr)))[[1]]

# Determine order of the randomized hierarchical addresses

   rord <- order(ranhadr)

   if(do.sample) {

# Select grid cells that get a sample point
        
      rstrt <- runif(1, 0, sint)
      ttl.wt <- c(0, cumsum(cel.wt[rord]))
      idx <- ceiling((ttl.wt - rstrt)/sint)
      smpdx <- .Call("pickGridCells", samplesize, as.integer(idx))
      rdx <- rord[smpdx]
      n.cells <- length(unique(rdx))
      if(length(rdx) > n.cells) {
         temp <- sum(sapply(split(rdx, rdx), length) > 1)
         warning(paste("\nOf the ", n.cells, " grid cells from which sample points were selected,\n", temp, " (", round(100*temp/n.cells, 1), "%) of the cells contained more than one sample point.\n", sep=""))
      }

# Pick sample point(s) in selected cells

      id <- selectpts(rdx, xc, yc, dx, dy, ptsframe)
      rho <- ptsframe[match(id, ptsframe$id), ]
   
   } else {

# Pick all points in the frame

      id <- selectframe(rord, xc, yc, dx, dy, ptsframe)
      rho <- ptsframe[match(id, ptsframe$id), ]
   }

# Construct sample hierarchical address

   np <- nrow(rho)
   nlev <- max(1, trunc(logb(np,4)))
   ifelse(np == 4^nlev, nlev, nlev <- nlev + 1)
   ad <- matrix(0, 4^nlev, nlev)
   rv4 <- 0:3
   pwr4 <- 4.^(0.:(nlev - 1.))
   for(i in 1:nlev)
      ad[, i] <- rep(rep(rv4, rep(pwr4[i], 4.)),pwr4[nlev]/pwr4[i])
   rho4 <- as.vector(ad%*%matrix(rev(pwr4), nlev, 1))

# Place sample in reverse hierarchical order

   rho <- rho[unique(floor(rho4 * np/4^nlev)) + 1.,]

# Assign Site ID

   siteID <- SiteBegin - 1 + 1:nrow(rho)

# Place Site ID as first column and add weights

   rho <- data.frame(siteID=siteID, id=rho$id, xcoord=rho$x, ycoord=rho$y,
      mdcaty=rho$mdcaty, wgt=1/rho$mdm)
   row.names(rho) <- 1:nrow(rho)

# Assign the final number of levels as an attribute of the output data frame

   attr(rho, "nlev") <- endlev

# Return the sample

   rho
}

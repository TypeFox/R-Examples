spbalance <- function(spsample, spframe = NULL, tess_ind = TRUE,
   sbc_ind = FALSE, nrows = 5, dxdy = TRUE) {

################################################################################
# Function: spbalance
# Purpose: Calculate spatial balance metrics for a survey design
# Programmer: Tom Kincaid
# Date: February 17, 2012
# Last Revised: April 17, 2015
# Description:      
#   This function calculates spatial balance metrics for a survey design.    Two
#  options for calculation of spatial balance metrics are available: (1) use
#  proportions obtained from the intersection of Dirichlet tesselation polygons
#  for the sample points with the frame object and (2) use proportions obtained
#  from a rectangular grid superimposed on the sample points and the frame
#  object.  In both cases the proportions are used to calculate the spatial
#  balance metrics.  Two metrics are calculated: (1) the Pielou evenness measure
#  and (2) the chi-square statistic.
# Arguments:
#   spsample = an object of class SpatialDesign produced by either the grts or
#     irs functions that contains survey design information and additional
#     attribute (auxiliary) variables.
#   spframe = an sp package object of class SpatialPointsDataFrame,
#     SpatialLinesDataFrame, or SpatialPolygonsDataFrame that contains the
#     survey design frame.  The default is NULL.
#   tess_ind = a logical variable indicating whether spatial balance metrics are
#     calculated using proportions obtained from the intersection of Dirichlet
#     tesselation polygons for the sample points with the frame object.  TRUE
#     means calculate the metrics.  FALSE means do not calculate the metrics.
#     The default is TRUE. 
#   sbc_ind = a logical variable indicating whether spatial balance metrics are
#     calculated using proportions obtained from a rectangular grid superimposed
#     on the sample points and the frame.  TRUE means calculate the metrics.
#     FALSE means do not calculate the metrics. The default is FALSE. 
#   nrows = number of rows (and columns) for the grid of cells.  The default is
#     5.
#   dxdy = indicator for equal x-coordinate and y-coordinate grid cell
#     increments, where TRUE means the increments are equal and FALSE means the
#     increments are not equal.  The default is TRUE.
# Results: 
#   A list containing the following components:
#     (1) tess - results for spatial balance metrics using tesselation polygons
#     (2) sbc - results for spatial balance metrics using a rectangular grid
#   If either the tess_ind or sbc_ind arguments are set to FALSE, the
#   corresponding component in the list is set to NULL.  Otherwise, each
#   components of the list is a lists that contains the following components:
#     (1) J_subp - Pielou evenness measure
#     (2) chi_sq - Chi-square statistic
#     (3) extent - frame extent for each  Dirichlet tesselation polygon or
#                  rectangular grid cell
#     (4) prop - frame proportion for each Dirichlet tesselation polygon or
#                rectangular grid cell
# Other Functions Required:
#   deldir - deldir package function that computes the Delaunay triangulation
#     and Dirichlet tesselation of a set of points.
#   tile.list - deldir package function that extracts coordinates of the
#     Dirichlet tesselation polygons from the object produced by the deldir
#     function.
#   gIntersection - rgeos package function that determines the intersection
#     between two sp package objects
#   LinesLength - sp package function that determines length of the line
#     segemnts in a class Lines object
#   sbcframe - function to calculate spatial balance grid cell extent and
#     proportions for a sample frame
#   sbcsamp - function to calculate spatial balance grid cell extent and
#     proportions for a survey design
# Example:
#   design <- list(Stratum1=list(panel=c(PanelOne=50), seltype="Equal",
#      over=10), Stratum2=list(panel=c(PanelOne=50, PanelTwo=50),
#      seltype="Unequal", caty.n=c(CatyOne=25, CatyTwo=25, CatyThree=25,
#      CatyFour=25), over=75))
#   frame <- read.shp("shapefile")
#   samp <- grts(design=design, DesignID="Test.Site", type.frame="area",
#      src.frame="shapefile", in.shape="shapefile", att.frame=frame@data,
#      stratum="stratum", mdcaty="mdcaty", shapefile=TRUE,
#      shapefilename="sample")
#   spbalance(samp, frame, sbc_ind = TRUE)
################################################################################

# Obtain the sample x-coordinates, y-coordinates, survey design weights,
# multidensity category values, and  stratum names from the spsample object
xcoord <- spsample@data$xcoord
ycoord <- spsample@data$ycoord
wgt <- spsample@data$wgt
mdcaty <- spsample@data$mdcaty
stratum <- spsample@data$stratum
n <- nrow(spsample@data)

# Determine the strata names
strata.names <- unique(stratum)

#
# Section for metrics calculted using Dirichlet tesselation polygons
#

if(tess_ind) {

# Determine whether an appropriate frame object was supplied

   if(is.null(spframe))
      stop("\nAn object containing the survey design frame must be supplied as the spframe \nargument.")
   if(!(class(spframe) %in% c("SpatialPointsDataFrame", "SpatialLinesDataFrame", "SpatialPolygonsDataFrame")))
      stop("\nThe spframe argument must be a member of class SpatialPointsDataFrame, \nSpatialLinesDataFrame, or SpatialPolygonsDataFrame.")

# Obtain the bounding box from the spframe object
   bbox <- c(spframe@bbox[1,], spframe@bbox[2,])

# Create an sp object containing the Dirichlet tesselation polygons for the
# sample points
   tiles <- tile.list(deldir(xcoord, ycoord, rw=bbox))
   sptess <- rep(list(NA), n)
   for(i in 1:n) {
      nv <- length(tiles[[i]]$x)
      sptess[[i]] <- SpatialPolygons(list(Polygons(list(Polygon(cbind(
         c(tiles[[i]]$x[1], tiles[[i]]$x[nv:1]),
         c(tiles[[i]]$y[1], tiles[[i]]$y[nv:1])))), 1)))
   }

# Determine the type of frame object
   temp <- class(spframe)
   ftype <- substr(temp, 8, nchar(temp) - 9)

# Intersect each Dirichlet tesselation polygon with the frame object and
# calculate extent and proportion
   extent <- numeric(n)
   for(i in 1:n) {
      if(ftype == "Points") {
         temp <- gIntersection(sptess[[i]], spframe)
         extent[i] <- nrow(temp@coords)
      } else if(ftype == "Lines") {
         temp <- sapply(spframe@lines, function(x) length(x@Lines))
         if(any(temp > 1)) {
            for(j in 1:length(spframe@lines)) {
               temp <- gIntersection(sptess[[i]],
                                     SpatialLines(list(spframe@lines[[j]])))
               extent[i] <- extent[i] + ifelse(is.null(temp), 0,
                  LinesLength(temp@lines[[1]]))
            }
         } else {
            temp <- gIntersection(sptess[[i]], spframe)
            extent[i] <- LinesLength(temp@lines[[1]])
         }
      } else if(ftype == "Polygons") {
         temp <- sapply(spframe@polygons, function(x) length(x@Polygons))
         if(any(temp > 1)) {
            for(j in 1:length(spframe@polygons)) {
               temp <- gIntersection(sptess[[i]],
                                  SpatialPolygons(list(spframe@polygons[[j]])))
               extent[i] <- extent[i] + ifelse(is.null(temp), 0,
                  temp@polygons[[1]]@area)
            }
         } else {
            temp <- gIntersection(sptess[[i]], spframe)
            extent[i] <- temp@polygons[[1]]@area
         }
      } else {
         stop(paste("'Spatial", ftype, "DataFrame' is not a known class of sp object.\n", sep=""))
      }
   }
   prop <- extent/sum(extent)

# Calculate the spatial balance metrics
   prob <- wgt/sum(wgt)
   J_subp <- sum(prop * log(prop))/sum(prob * log(prob))
   chi_sq <- sum(((prop - prob)^2)/prob)

# Print the spatial balance metrics
   cat("\nSpatial Balance Metric using Dirichlet Tesselation Polygons:\n")
   cat(paste("   Pielou evenness measure:", J_subp, "\n"))
   cat(paste("   Chi-square statistic:", chi_sq, "\n"))

# Create the output list
   tess <- list(J_subp=J_subp, chi_sq=chi_sq, extent=extent, prop=prop)

# Metrics calculated using Dirichlet tesselation polygons were not requested
} else {
   tess <- NULL
}

#
# Section for metrics calculted using a rectangular grid
#

if(sbc_ind) {

# Calculate grid cell extent and proportion for the frame
   sbc.frame <- sbcframe(spframe = spframe, nrows = nrows, dxdy = dxdy)

# Calculate grid cell extent and proportion for the sample
   sbc.sample <- sbcsamp(spsample, sbc.frame)      

# Calculate the spatial balance metrics
   prop_f <- sbc.frame$prop[sbc.frame$prop != 0]
   prop_s <- sbc.sample$prop[sbc.sample$prop != 0]
   J_subp <- sum(prop_s * log(prop_s))/sum(prop_f * log(prop_f))
   ind <- sbc.frame$prop != 0
   prop_f <- sbc.frame$prop[ind]
   prop_s <- sbc.sample$prop[ind]
   chi_sq <- sum(((prop_s - prop_f)^2)/prop_f)

# Print the spatial balance metrics
   cat("\n\nSpatial Balance Metric using a Rectangular Grid:\n")
   cat(paste("   Pielou evenness measure:", J_subp, "\n"))
   cat(paste("   Chi-square statistic:", chi_sq, "\n"))

# Create the output list
   sbc <- list(J_subp=J_subp, chi_sq=chi_sq, extent_f=sbc.frame$extent,
               prop_f=sbc.frame$prop, extent_s=sbc.sample$extent,
               prop_s=sbc.sample$prop)

# Metrics calculated using a rectangular grid were not requested
} else {
   sbc <- NULL
}

# Return results
list(tess=tess, sbc=sbc)
}

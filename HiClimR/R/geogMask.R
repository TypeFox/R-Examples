# $Id: geogMask.R, v1.2.3 2015/08/05 12:00:00 hsbadr EPS JHU              #
#-------------------------------------------------------------------------#
# This function is a part of HiClimR R package.                           #
#-------------------------------------------------------------------------#
#  HISTORY:                                                               #
#-------------------------------------------------------------------------#
#  Version  |  Date      |  Comment   |  Author          |  Email         #
#-------------------------------------------------------------------------#
#           |  May 1992  |  Original  |  F. Murtagh      |                #
#           |  Dec 1996  |  Modified  |  Ross Ihaka      |                #
#           |  Apr 1998  |  Modified  |  F. Leisch       |                #
#           |  Jun 2000  |  Modified  |  F. Leisch       |                #
#-------------------------------------------------------------------------#
#   1.0.0   |  03/07/14  |  HiClimR   |  Hamada S. Badr  |  badr@jhu.edu  #
#   1.0.1   |  03/08/14  |  Updated   |  Hamada S. Badr  |  badr@jhu.edu  #
#   1.0.2   |  03/09/14  |  Updated   |  Hamada S. Badr  |  badr@jhu.edu  #
#   1.0.3   |  03/12/14  |  Updated   |  Hamada S. Badr  |  badr@jhu.edu  #
#   1.0.4   |  03/14/14  |  Updated   |  Hamada S. Badr  |  badr@jhu.edu  #
#   1.0.5   |  03/18/14  |  Updated   |  Hamada S. Badr  |  badr@jhu.edu  #
#   1.0.6   |  03/25/14  |  Updated   |  Hamada S. Badr  |  badr@jhu.edu  #
#-------------------------------------------------------------------------#
#   1.0.7   |  03/30/14  |  Hybrid    |  Hamada S. Badr  |  badr@jhu.edu  #
#   1.0.8   |  05/06/14  |  Updated   |  Hamada S. Badr  |  badr@jhu.edu  #
#-------------------------------------------------------------------------#
#   1.0.9   |  05/07/14  |  CRAN      |  Hamada S. Badr  |  badr@jhu.edu  #
#   1.1.0   |  05/15/14  |  Updated   |  Hamada S. Badr  |  badr@jhu.edu  #
#   1.1.1   |  07/14/14  |  Updated   |  Hamada S. Badr  |  badr@jhu.edu  #
#   1.1.2   |  07/26/14  |  Updated   |  Hamada S. Badr  |  badr@jhu.edu  #
#   1.1.3   |  08/28/14  |  Updated   |  Hamada S. Badr  |  badr@jhu.edu  #
#   1.1.4   |  09/01/14  |  Updated   |  Hamada S. Badr  |  badr@jhu.edu  #
#   1.1.5   |  11/12/14  |  Updated   |  Hamada S. Badr  |  badr@jhu.edu  #
#-------------------------------------------------------------------------#
#   1.1.6   |  03/01/15  |  GitHub    |  Hamada S. Badr  |  badr@jhu.edu  #
#-------------------------------------------------------------------------#
#   1.2.0   |  03/27/15  |  MVC       |  Hamada S. Badr  |  badr@jhu.edu  #
#   1.2.1   |  05/24/15  |  Updated   |  Hamada S. Badr  |  badr@jhu.edu  #
#   1.2.2   |  07/21/15  |  Updated   |  Hamada S. Badr  |  badr@jhu.edu  #
#   1.2.3   |  08/05/15  |  Updated   |  Hamada S. Badr  |  badr@jhu.edu  #
#-------------------------------------------------------------------------#
# COPYRIGHT(C) 2013-2015 Earth and Planetary Sciences (EPS), JHU.         #
#-------------------------------------------------------------------------#
# Function: Geographic mask from longitude and latitute                   #
#-------------------------------------------------------------------------#

# Function: Geographic mask for an area from longitude and latitute
geogMask <- function(continent = NULL, region = NULL, country = NULL, 
    lon = NULL, lat = NULL, InDispute = TRUE, verbose = TRUE, 
    plot = FALSE, colPalette = NULL, pch = 15, cex = 1) {
    
    # Get World mask from LazyData
    wMask <- get("WorldMask", envir = .GlobalEnv)

    if (verbose) write("---> Checking geographic masking options...", "")
    if (is.null(continent) && is.null(region) && is.null(country)) {
        gMask <- list()
        gMask$continent <- unique(sort(wMask$info[, 7]))
        gMask$region <- unique(sort(wMask$info[, 6]))
        gMask$country <- unique(as.character(sort(wMask$info[!is.na(wMask$info[, 
            3]), 3])))
    } else {
        if (verbose) write("---> Checking geographic masking data...", "")

        if (is.null(lon) || is.null(lat) || length(lon) != length(lat)) 
            stop("invalid coordinates")
        
        if (!is.null(continent)) {
            if (verbose) write("---> Geographic masking by continent...", "")

            area <- which(wMask$info[, 7] %in% continent)
            
            if (length(area) < 1) 
                stop("invalid continent")
            
        } else if (!is.null(region)) {
            if (verbose) write("---> Geographic masking by region...", "")

            area <- which(wMask$info[, 6] %in% region)
            
            if (length(area) < 1) 
                stop("invalid region")
            
        } else if (!is.null(country)) {
            if (verbose) write("---> Geographic masking by country...", "")

            area <- which(wMask$info[, 3] %in% country)
            
            if (length(area) < 1) 
                stop("invalid country")
            
            # Fix confusing country codes/names
            #for (i in 1:length(area)) {
            #    area <- union(area, which(gregexpr(pattern = wMask$info[area[i], 
            #      3], wMask$info[, 1]) != -1))
            #}
            
            # Areas in dispute
            InDisputeArea <- which(gregexpr(pattern = "In dispute", wMask$info[, 
                1]) != -1)
            if (InDispute) {
                for (i in 1:length(area)) {
                    area <- union(area, InDisputeArea[which(grepl(wMask$info[area[i], 
                    1], wMask$info[InDisputeArea, 1]))])
                }
            }
        }
        
        dx <- as.numeric(rownames(wMask$mask))[2] - as.numeric(rownames(wMask$mask))[1]
        dy <- as.numeric(colnames(wMask$mask))[2] - as.numeric(colnames(wMask$mask))[1]
        
        rx <- abs(as.integer(round(log(dx, 10))))
        ry <- abs(as.integer(round(log(dy, 10))))
        
        i <- (round(lon, rx) - as.numeric(rownames(wMask$mask))[1])/dx + 
            1
        j <- (round(lat, ry) - as.numeric(colnames(wMask$mask))[1])/dy + 
            1
        
        # gMask <- seq(1,length(lon))[-which(diag(wMask$mask[i,j]) %in% area)]
        gMask <- NULL
        for (nn in 1:length(lon)) {
            nnMask <- ifelse(is.na(wMask$mask[i[nn], j[nn]]), -999, wMask$mask[i[nn], 
                j[nn]])
            if (!any(area == nnMask)) {
                gMask <- c(gMask, nn)
            }
        }
    }
    
    if (plot) {
        if (verbose) write("---> Generating geographic mask map...", "")        
        Regions <- rep(1, length(lon))
        Regions[gMask] <- NA
            
        if (is.null(colPalette)) {
            colPalette <- colorRampPalette(c("#8DD3C7", "#FFFFB3", 
              "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", 
              "#FCCDE5", "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F"))
        }
        Longitude <- lon
        Latitude <- lat
        plot(Longitude, Latitude, col = Regions, pch = pch, cex = cex)
    }
    
    #gc()
    return(gMask)
}

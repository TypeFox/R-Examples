rnd_init_maps<-function(){
# Adapted from:
## Simulated categorical habitat/landscape pattern (simulated habitat distribution in the landscape)
## Implementation of the Modified Random Cluster (MRC) method of Saura and Martinez-Millan (2000) Landscape Ecology 15: 661-678. 

## The paper by Saura and Martinez-Millan (2000), the SIMMAP software where the MRC method was initially implemented
## and the SIMMAP user manual are available at http://www2.montes.upm.es/personales/saura/software.html#simmap. 

## Read that paper (particularly the description of steps A-D there available) and the SIMMAP user manual for  
## better understanding this implementation and the impact of the simulation parameters on the resultant patterns.

## The author of this R implementation is Murray Efford (University of Auckland, New Zealand), 2012-04-09,10,11,12
## Santiago Saura (Universidad Politecnica de Madrid, Spain) has made slight modifications (minor aspects only) 2012-06-06.

## The original R implementation by Murray Efford is named "randomHabitat" and is part of 
## the "secr" package available at http://cran.r-project.org/web/packages/secr/index.html

## This implementation Requires 'raster' and 'igraph' packages:

## Robert J. Hijmans & Jacob van Etten (2011). raster: Geographic
##  analysis and modeling with raster data. R package version 1.9-33.
##  http://CRAN.R-project.org/package=raster

## Csardi G, Nepusz T: The igraph software package for complex network
##  research, InterJournal, Complex Systems 1695. 2006.
##  http://igraph.sf.net

## Restrictions in this implementation.
##   1) adjacency for step 'B' is based on a 1-step rook move (direction = 4) 
##   (i.e. 4 neighbours are considered for defining habitat patches)
##   2) Only two classes in the landscape pattern are considered (habitat vs non-habitat)
## It would be easy to adapt the R code to make it more general as related to 1) and 2).
## See the SIMMAP sofware package for generating patters without these restrictions.

## Input parameters (values to be specified below)
## 'Lx' and 'Ly' size of the rectangular pattern as the number of pixels in the vertical and horizontal directions
## 'p' initial probability (controls the degree of habitat fragmentation) 
## 'A' relative area (proportion) of the landscape occupied by habitat 
## 'plt' logical for whether to plot resultant patterns (both intermediate steps and final pattern)

Lx <- 113
Ly <- 114
p <- 0.5
A <- 0.7
plt <- TRUE

## This is fixed in the current code (adjacency for step 'B' is based on a 1-step rook move (direction = 4))
dirs <- 4

n<- Lx*Ly

        ## create rasterLayer
        #require(raster)
	      layer <- raster::raster(nrows = Lx, ncols = Ly, xmn = 688075, xmx = 690925, ymn = 4173950, ymx = 4176775)

        ## step A. Percolation map generation
        raster::values(layer) <- rep(0, Lx*Ly)
        raster::values(layer)[sample.int(Lx*Ly, round(Lx*Ly*p))] <- 1   ## as close as possible
        #if (plt) plot(layer, useRaster = FALSE, col=c("white","black")) ## plots selected cells in black (cells with a value of in the percolation map)
        ## the cells with a value of 1 in the percolation map will be distributed among non-habitat (class 1) and habitat (class 2) in step C
        
        ## step B. Cluster identification (clustering of adjoining pixels)
        clumped <- raster::clump(layer, directions=dirs, gaps = FALSE)
        
        ## step C. Cluster type assignment
        ## ncluster is the number of clusters with a value of 1 in the percolation map
        ncluster <- max(raster::values(clumped), na.rm = TRUE)
        types <- factor(c(0,1))          # non-habitat, habitat
        numTypes <- as.numeric(types)    # 1(non-habitat),2(habitat)
        ## this now assigns the clusters either to 1 (non-habitat) or 2 (habitat)
        clustertype <- sample(numTypes, ncluster, replace = TRUE, prob = c(1-A,A))
        raster::values(clumped) <- clustertype[raster::values(clumped)]
        #if (plt) plot(clumped, useRaster = FALSE, col=c("red","green")) ## in black the provisional non-habitat (class 1) and in green the provisional habitat (class2)
        ## this provisional assignment will be extended to the full image (the pixels with 0 in the inital percolation map) in step D.

        ## step D. Filling in image
        cellsUnassigned <- (1:n)[is.na(raster::values(clumped))]
        cellsAssigned <- (1:n)[!is.na(raster::values(clumped))]
        tempadj <- raster::adjacent(clumped, cellsUnassigned, cellsAssigned, directions = 8)
        tempadj <- split(tempadj[,2], tempadj[,1])
        fillinType <- function (adjcells) {
            type <- raster::values(clumped)[adjcells]
            type <- factor(type)
            freq <- tabulate(type)
            result <- as.numeric(levels(type)[freq == max(freq)])
            if (length(result)>1)
                result <- sample (result, 1)
            result
        }
        # cells with typed neighbours
        filled <- sapply(tempadj, fillinType)
        filledCells <- as.numeric(names(filled))
        raster::values(clumped)[filledCells] <- filled
        # cells with no typed neighbours
        notfilledCells <- cellsUnassigned[!(cellsUnassigned %in% filledCells)]
        randomType <- sample(numTypes, length(notfilledCells), replace = TRUE, prob = c(1-A,A))
        raster::values(clumped)[notfilledCells] <- randomType
        #if (plt) plot(clumped, useRaster = FALSE, col=c("red","green"))

## second map

	      layer <- raster::raster(nrows = Lx, ncols = Ly, xmn = 688075, xmx = 690925, ymn = 4173950, ymx = 4176775)

        ## step A. Percolation map generation
        raster::values(layer) <- rep(0, Lx*Ly)
        raster::values(layer)[sample.int(Lx*Ly, round(Lx*Ly*p))] <- 1   ## as close as possible
        #if (plt) plot(layer, useRaster = FALSE, col=c("white","black")) ## plots selected cells in black (cells with a value of in the percolation map)
        ## the cells with a value of 1 in the percolation map will be distributed among non-habitat (class 1) and habitat (class 2) in step C
        
        ## step B. Cluster identification (clustering of adjoining pixels)
        clumped2 <- raster::clump(layer, directions=dirs, gaps = FALSE)
        
        ## step C. Cluster type assignment
        ## ncluster is the number of clusters with a value of 1 in the percolation map
        ncluster <- max(raster::values(clumped2), na.rm = TRUE)
        types <- factor(c(0,1))          # non-habitat, habitat
        numTypes <- as.numeric(types)    # 1(non-habitat),2(habitat)
        ## this now assigns the clusters either to 1 (non-habitat) or 2 (habitat)
        clustertype <- sample(numTypes, ncluster, replace = TRUE, prob = c(1-A,A))
        raster::values(clumped2) <- clustertype[raster::values(clumped2)]
       # if (plt) plot(clumped2, useRaster = FALSE, col=c("red","green")) ## in black the provisional non-habitat (class 1) and in green the provisional habitat (class2)
        ## this provisional assignment will be extended to the full image (the pixels with 0 in the inital percolation map) in step D.

        ## step D. Filling in image
        cellsUnassigned <- (1:n)[is.na(raster::values(clumped2))]
        cellsAssigned <- (1:n)[!is.na(raster::values(clumped2))]
        tempadj <- raster::adjacent(clumped2, cellsUnassigned, cellsAssigned, directions = 8)
        tempadj <- split(tempadj[,2], tempadj[,1])
        fillinType <- function (adjcells) {
            type <- raster::values(clumped2)[adjcells]
            type <- factor(type)
            freq <- tabulate(type)
            result <- as.numeric(levels(type)[freq == max(freq)])
            if (length(result)>1)
                result <- sample (result, 1)
            result
        }
        # cells with typed neighbours
        filled <- sapply(tempadj, fillinType)
        filledCells <- as.numeric(names(filled))
        raster::values(clumped2)[filledCells] <- filled
        # cells with no typed neighbours
        notfilledCells <- cellsUnassigned[!(cellsUnassigned %in% filledCells)]
        randomType <- sample(numTypes, length(notfilledCells), replace = TRUE, prob = c(1-A,A))
        raster::values(clumped2)[notfilledCells] <- randomType
        #if (plt) plot(clumped2, useRaster = FALSE, col=c("red","green"))

c3<<-NULL
c3<<-clumped + clumped2
c3<-get('c3')
}        

        
   


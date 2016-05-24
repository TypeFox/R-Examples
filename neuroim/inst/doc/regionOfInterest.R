## ---- echo = FALSE, message = FALSE--------------------------------------
knitr::opts_chunk$set(collapse = T, comment = "#>")
library(neuroim)

## ------------------------------------------------------------------------
      # attach MNI BrainSpace instance
      
      data("MNI_SPACE_1MM")
      
      # we create a spherical ROI centered around voxel coordinates [20,20,20] with a 5mm radius, 
      # filling all values in the ROI with 100.

      sphere <- RegionSphere(MNI_SPACE_1MM, c(20,20,20), radius=5, fill=100)
      
      # to extract the voxel coordinates of the sphere:
      
      vox <- coords(sphere)
      
      # to get the values at the coordinate locations
      
      vals <- values(sphere)
      all.equal(vals, rep(100, length(vals)))   

## ------------------------------------------------------------------------

    
    
    rpoint <- c(-50,-28,10)
    
    # Because RegionSphere takes a coordinate in voxel units, 
    # we need to convert to the real-world MNI coordinate to grid coordinates.
    
    vox <- coordToGrid(MNI_SPACE_1MM, rpoint)
    sphere <- RegionSphere(MNI_SPACE_1MM, vox, radius=10, fill=1)
    dim(coords(sphere))
    
    # convert back to MNI coordinates
    
    mnicoords <- indexToCoord(MNI_SPACE_1MM, indices(sphere))
    
    ## compute center of mass of MNI coords in ROI (should be close to original coordinate)
    centerOfMass <- colMeans(mnicoords)
    centerOfMass
    

## ------------------------------------------------------------------------
    
    sphere <- RegionSphere(MNI_SPACE_1MM, c(50,50,50), radius=10, fill=1)
    sparsevol <- SparseBrainVolume(values(sphere),MNI_SPACE_1MM,indices=indices(sphere))
    
    sum(sparsevol) == sum(values(sphere))
    
    all(dim(sparsevol) == dim(MNI_SPACE_1MM))
    
    
    


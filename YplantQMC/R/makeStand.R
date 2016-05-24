#'Make a stand of virtual plants
#'
#'@description Make a stand of plants, for use in \code{runYplant}, and for visualization.
#'See the example below to get started. Support for \code{runYplant} is
#'somewhat experimental, and \code{YplantDay} is not supported yet. Proceed at
#'your own risk.
#'
#'@details
#'The \code{xyz} argument must be a dataframe or matrix with three columns, and
#'it is assumed to be in the order X,Y,Z.
#'
#'The \code{plotbox} argument is optional, if it is not provided the plot
#'boundary will be as a rectangle that just fits around the projected crown
#'area. In some cases, the base of the stem can thus fall outside the plot
#'boundary. For now, the plot boundary is only used to calculate the leaf area
#'index, which has no bearing on any simulation results.
#'
#'@param plants List of plants to be placed in the stand.
#'@param xyz Data frame (or matrix) with x,y,z locations of the stem positions
#'of the plants.
#'@param plotbox Optional. Plot boundary, used for scaling-up purposes.
#'@return An object of class \code{stand3d}, methods exist for \code{print},
#'\code{plot}, \code{runYplant}. And soon, \code{YplantDay}.
#'@author Remko Duursma
#'@keywords misc
#'@examples
#'
#'
#'\dontrun{
#'
#'# Make a stand consisting of three 'toona' plants.
#'toonastand <- makeStand(list(toona,toona,toona),
#'                        xyz=data.frame(x=c(0,200,100),
#'                                       y=c(50,50,300),
#'                                       z=c(0,0,0)))
#'
#'# The print method shows a very short summary:
#'toonastand
#'
#'# Plot the stand
#'plot(toonastand)
#'
#'}
#'
#'@export
makeStand <- function(plants=list(),
                      xyz=data.frame(x=0,y=0,z=0),
                      plotbox=NULL
){
  
  
  if(nrow(xyz) != length(plants))
    stop("Must provide X,Y,Z of stem position for each plant.")
  
  
  # Shift plants.
  for(i in 1:length(plants)){
    plants[[i]] <- shiftplant(plants[[i]], xyz[i,1], xyz[i,2], xyz[i,3])
  }
  
  # Find plotbox, if not provided already.
  if(is.null(plotbox)){
    
    maxx <- maxy <- minx <- miny <- 0
    for(i in 1:length(plants)){
      minx <- min(minx, min(plants[[i]]$leafbasecoor[,1]))
      minx <- min(minx, min(plants[[i]]$leaftipcoor[,1]))    
      
      maxx <- max(maxx, max(plants[[i]]$leafbasecoor[,1]))
      maxx <- max(maxx, max(plants[[i]]$leaftipcoor[,1])) 
      
      miny <- min(miny, min(plants[[i]]$leafbasecoor[,2]))
      miny <- min(miny, min(plants[[i]]$leaftipcoor[,2]))    
      
      maxy <- max(maxy, max(plants[[i]]$leafbasecoor[,2]))
      maxy <- max(maxy, max(plants[[i]]$leaftipcoor[,2]))   
    }
    plotbox <- c(minx,miny,maxx,maxy)
    
  }
  
  # LAI:
  la <- c()
  for(i in 1:length(plants)){
    la[i] <- plants[[i]]$leafarea*10^6
  }
  LA <- sum(la)
  b <- plotbox
  area <- (b[3] - b[1]) * (b[4] - b[2])
  LAI <- LA / area
  
  l <- list()
  l$plants <- plants
  l$nplants <- length(plants)
  l$nleaves <- sapply(plants, "[[","nleaves")
  l$xyz <- xyz
  l$plotbox <- plotbox
  l$LAI <- LAI
  
  ld <- lapply(plants, "[[", "leafdata")
  l$leafdata <- do.call(rbind,ld)
  
  class(l) <- "stand3d"
  return(l)
}

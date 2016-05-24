Df2SpLinesDf <-
function( spLines, df, add.distance=F, add.azimuth=F )
{
  # This function converts an object of class SpatialLines, calculated by the 
  # function Df2SpLines, into an Object of class SpatialLinesDataFrame.
  #
  # Args:
  #  spLines: Object of class SpatialLines calculated by the function Df2SpLines.
  #  df: Data Frame Object created by the function ProcTraj.
  #  add.distance: Logical: If True, it will calculate and include the distance in meters between the first and last point for every line.
  #  add.azimuth: Logical: If True it will calculate and include the azimuth for every line.
  
  # Results:
  #   Returns an object of class SpatialLinesDataFrame.
  
  # get the trajectory lenghth
  # all trajectories have the same length
  
  max.traj.length <- max(abs(df$hour.inc)) + 1
  
  # create a traj ID column to identify each trajctory uniquely 
  df['ID'] <- rep(1:(nrow(df) / max.traj.length), each=max.traj.length )
  
  # apply the function to each subgroup of the data frame
  # the dataframe is divided in subgroups of equal IDs 
  #(each trajectory has an unique ID)
  # the function just get the fist line of each trajectory and returns 
  # a data.frame with that information
  data.list <- ddply(df, 'ID', function(df){ df[1,] }, .inform=TRUE)
  
  if(add.distance == T){
    CalcDistance <- function( line ){ 
      # get the coordinates of the point
      cc <- as.data.frame(coordinates(line))
      
      # get the first and last pair of coordinates
      cc <- cc[-c(2,3),]
      
      # calculate the distance between those two points
      dist <- spDists(as.matrix(cc), longlat=TRUE)[1,2]
    }
    
    data.list$distance <- sapply(slot(spLines, "lines"), FUN=CalcDistance)
    data.list$distance <- data.list$distance * 1000
  }
  
  if(add.azimuth==T){
    CalcAzimuth <- function( line ){ 
      # get the coordinates of the point
      cc <- as.data.frame(coordinates(line))
      
      # get the first and last pair of coordinates
      first.p <- as.matrix(cc[1,])
      second.p <-  as.matrix(cc[nrow(cc),])
      
      # calculate the distance between those two points
      gzAzimuth(first.p, second.p)
    }
    
    data.list$azimuth <- sapply(slot(spLines, "lines"), FUN=CalcAzimuth)
  }
  
  spLinesDataFrame <- SpatialLinesDataFrame(spLines, data = data.list)
  
  spLinesDataFrame
}

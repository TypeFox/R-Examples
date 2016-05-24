Df2SpLines <-
function( df, crs=NA )
  {   
    # This function converts an object of type data frame, calculated by the function
    # ProcTraj, into an object of type Spatial Lines.
    #
    # Args:
    #   df: Data Frame Object created by the function ProcTraj.
    #   crs: String: Valid projection string. An example would be crs= "+proj=longlat +datum=NAD27"
    #
    # Results:
    #  Returns an object of class SpatialLines.
    
    if( !is.na(crs) ){
      crs <- CRS(crs)
    }
    
    max.traj.length <- max(abs(df$hour.inc)) + 1
    
    if(nrow(df) %% max.traj.length != 0) {
      stop("The number of rows in the 'df' argument is not a multiple of the length of an individual trajectory" )
    }
    
    # create a traj ID column to identify each trajctory uniquely 
    df['ID'] <- rep(1:(nrow(df) / max.traj.length), each=max.traj.length )
    
    CreateLines <- function(df) {
      # get the coordinates out of the data.frame
      cc <- df[7:8]
      
      # reverse the order of the columns from [Lat Long] to [Long Lat]
      cc <- cc[, c(2, 1)]
      
      # create a individual line
      line <- Line(cc)
      
      # transfor the line [line] into a Lines object and assign a unique ID
      Lines(line, ID=as.character(df$ID[1]))  
    }
    
    lines.list <- dlply( df, 'ID', CreateLines)
    
    sp.lines <- SpatialLines(lines.list, proj4string = crs)
    
    sp.lines
  }

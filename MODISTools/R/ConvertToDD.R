ConvertToDD <-
function(XY, FileSep = NULL, LatColName, LongColName)
{
    if(!is.object(XY) & !is.character(XY)) stop("XY must be an object in R or a file path character string.")
    if(is.object(XY)) XY<- data.frame(XY)
    if(is.character(XY)){
      if(!file.exists(XY)) stop("Character string input for XY argument does not resemble an existing file path.")
      if(is.null(FileSep)) stop("To load a file as input, you must also specify its delimiter (FileSep).")
      XY<- read.delim(XY, sep = FileSep) 
    }
    
    DMS.lat <- as.character(XY[ ,which(names(XY) == LatColName)])
    DMS.long <- as.character(XY[ ,which(names(XY) == LongColName)])  
    
    which.format.lat <- gregexpr("([^0-9.][0-9])", DMS.lat)
    which.format.long <- gregexpr("([^0-9.][0-9])", DMS.long)
    DM.or.DMS.lat <- rep(NA, nrow(XY))
    DM.or.DMS.long <- rep(NA, nrow(XY))
  
    for(i in 1:nrow(XY)){ 
      DM.or.DMS.lat[i] <- length(which.format.lat[[i]]) 
      DM.or.DMS.long[i] <- length(which.format.long[[i]]) 
    }
    if(any(DM.or.DMS.lat != DM.or.DMS.long)){
      stop("A coordinate has been recognised with inconsistent formatting between lat and long. 
           Check for erroneous non-numeric characters. See the help page for advice on correct formats.")
    }
    if(any(DM.or.DMS.lat != 1 & DM.or.DMS.lat != 2)){
      stop("A coordinate has been found that does not match the required format for degrees minutes seconds or degrees minutes. 
           Check for erroneous non-numeric characters. See the help page for advice on correct formats.")
    }
    
    DD.lat <- rep(NA, nrow(XY))
    DD.long <- rep(NA, nrow(XY))
    D.lat <- rep(NA, nrow(XY))
    D.long <- rep(NA, nrow(XY))
    M.lat <- rep(NA, nrow(XY))
    M.long <- rep(NA, nrow(XY))
    S.lat <- rep(NA, nrow(XY))
    S.long <- rep(NA, nrow(XY))
    D.point.lat <- regexpr("[^0-9][0-9]{1,2}[^0-9]", DMS.lat)
    D.point.long <- regexpr("[^0-9][0-9]{1,2}[^0-9]", DMS.long)
    
    for(i in 1:nrow(XY)){     
      # For degrees minutes seconds coordinates.
      if(DM.or.DMS.lat[i] == 2){
        # Latitude
        D.lat[i] <- as.numeric(substr(DMS.lat[i], 1, D.point.lat[i]-1))
        M.lat[i] <- as.numeric(substr(DMS.lat[i], D.point.lat[i]+1, D.point.lat[i]+attr(D.point.lat, "match.length")[i]-2))
        if(substr(DMS.lat[i], nchar(DMS.lat[i]), nchar(DMS.lat[i])) == 'N'){
          S.lat[i] <- as.numeric(substr(DMS.lat[i], D.point.lat[i]+attr(D.point.lat, "match.length")[i], nchar(DMS.lat[i])-2)) 
        } else {
          if(substr(DMS.lat[i], nchar(DMS.lat[i]), nchar(DMS.lat[i])) == 'S'){ 
            D.lat[i] <- -D.lat[i]
            S.lat[i] <- as.numeric(substr(DMS.lat[i], D.point.lat[i]+attr(D.point.lat, "match.length")[i], nchar(DMS.lat[i])-2)) 
          } else {
            S.lat[i] <- as.numeric(substr(DMS.lat[i], D.point.lat[i]+attr(D.point.lat, "match.length")[i], nchar(DMS.lat[i])-1))
          }
        }
        # Calculate latitude decimal degrees.
        if(D.lat[i] >= 0){
          DD.lat[i] <- D.lat[i] + (M.lat[i] / 60) + (S.lat[i] / 3600)
        } else {
          DD.lat[i] <- -(S.lat[i] / 3600) - (M.lat[i] / 60) + D.lat[i]
        }
        if(substr(DMS.lat[i], nchar(DMS.lat[i]), nchar(DMS.lat[i])) == 'S' & D.lat[i] == 0){ 
          DD.lat[i] <- -DD.lat[i] 
        }
        
        # Longitude
        D.long[i] <- as.numeric(substr(DMS.long[i], 1, D.point.long[i]-1))
        M.long[i] <- as.numeric(substr(DMS.long[i], D.point.long[i]+1, D.point.long[i]+attr(D.point.long, "match.length")[i]-2))
        if(substr(DMS.long[i], nchar(DMS.long[i]), nchar(DMS.long[i])) == 'E'){ 
          S.long[i] <- as.numeric(substr(DMS.long[i], D.point.long[i]+attr(D.point.long, "match.length")[i], nchar(DMS.long[i])-2)) 
        } else {
          if(substr(DMS.long[i], nchar(DMS.long[i]), nchar(DMS.long[i])) == 'W'){ 
            S.long[i] <- as.numeric(substr(DMS.long[i], D.point.long[i]+attr(D.point.long, "match.length")[i], nchar(DMS.long[i])-2))
            D.long[i] <- -D.long[i] 
          } else {
            S.long[i] <- as.numeric(substr(DMS.long[i], D.point.long[i]+attr(D.point.long, "match.length")[i], nchar(DMS.long[i])-1))
          }
        }          
        # Calculate longitude decimal degrees.
        if(D.long[i] >= 0){
          DD.long[i] <- D.long[i] + (M.long[i] / 60) + (S.long[i] / 3600)
        } else {
          DD.long[i] <- -(S.long[i] / 3600) - (M.long[i] / 60) + D.long[i]
        }
        if(substr(DMS.long[i], nchar(DMS.long[i]), nchar(DMS.long[i])) == 'W' & D.long[i] == 0){ 
          DD.long[i] <- -DD.long[i] 
        }
      }
      
      # For degrees minutes coordinates.
      if(DM.or.DMS.lat[i] == 1){
        # Latitude
        D.lat[i] <- as.numeric(substr(DMS.lat[i], 1, D.point.lat[i]-1))
        if(substr(DMS.lat[i], nchar(DMS.lat[i]), nchar(DMS.lat[i])) == 'N'){ 
          M.lat[i] <- as.numeric(substr(DMS.lat[i], D.point.lat[i]+1, nchar(DMS.lat[i])-2)) 
        } else {
          if(substr(DMS.lat[i], nchar(DMS.lat[i]), nchar(DMS.lat[i])) == 'S'){ 
            M.lat[i] <- as.numeric(substr(DMS.lat[i], D.point.lat[i]+1, nchar(DMS.lat[i])-2))
            D.lat[i] <- -D.lat[i] 
          } else {
            M.lat[i] <- as.numeric(substr(DMS.lat[i], D.point.lat[i]+1, nchar(DMS.lat[i])-1))
          }
        }
        # Calculate latitude decimal degrees.
        if(D.lat[i] >= 0){
          DD.lat[i] <- D.lat[i] + (M.lat[i] /60)
        } else {
          DD.lat[i] <- -(M.lat[i] / 60) + D.lat[i]
        }
        if(substr(DMS.lat[i], nchar(DMS.lat[i]), nchar(DMS.lat[i])) == 'S' & D.lat[i] == 0){ 
          DD.lat[i]<- -DD.lat[i] 
        }
        
        # Longitude        
        D.long[i] <- as.numeric(substr(DMS.long[i], 1, D.point.long[i]-1))
        
        if(substr(DMS.long[i], nchar(DMS.long[i]), nchar(DMS.long[i])) == 'E'){   
          M.long[i] <- as.numeric(substr(DMS.long[i], D.point.long[i]+1, nchar(DMS.long[i])-2))
        } else {
          if(substr(DMS.long[i], nchar(DMS.long[i]), nchar(DMS.long[i])) == 'W'){   
            M.long[i] <- as.numeric(substr(DMS.long[i], D.point.long[i]+1, nchar(DMS.long[i])-2))
            D.long[i] <- -D.long[i] 
          } else {
            M.long[i] <- as.numeric(substr(DMS.long[i], D.point.long[i]+1, nchar(DMS.long[i])-1))
          }
        }
        # Calculate longitude decimal degrees.
        if(D.long[i] >= 0){
          DD.long[i] <- (D.long[i]) + ((M.long[i])/60)
        } else {
          DD.long[i] <- -((M.long[i])/60) + (D.long[i])
        }
        if(substr(DMS.long[i], nchar(DMS.long[i]), nchar(DMS.long[i])) == 'W' & D.long[i] == 0){ 
          DD.long[i] <- -DD.long[i] 
        }
      } 
    }
    
    # Checks that lat answers are going to be sensible before returning result.
    if(any(abs(D.lat) > 90)){
      cat("Invalid degrees of latitude entries:", "\n", XY[which(abs(D.lat) > 90), ], "\n")
      stop("Range of valid degrees is from -90 to 90.")
    }
    if(any(M.lat < 0 & M.lat > 60)){
      cat("Invalid minutes entries:", "\n", XY[which(M.lat > 0 & M.lat < 60), ], "\n")
      stop("Range of valid minutes is from 0 to 60.")
    }
    if(any(DM.or.DMS.lat == 2)){
      if(any(S.lat[which(DM.or.DMS.lat == 2)] < 0 & S.lat[which(DM.or.DMS.lat == 2)] > 60)){
        cat("Invalid seconds entries:", "\n",
            XY[which(S.lat[which(DM.or.DMS.lat == 2)] > 0 & S.lat[which(DM.or.DMS.lat == 2)] < 60), ], "\n")
        stop("Range of valid seconds is from 0 to 60.")
      }
    }  
    # Checks that long answers are going to be sensible before returning result.
    if(any(abs(D.long) > 180)){
      cat("Invalid degrees of longitude entries:", "\n", XY[which(abs(D.long) > 180), ], "\n")
      stop("Range of valid degrees longitude is from -180 to 180.")
    }
    if(any(M.long < 0 & M.long > 60)){
      cat("Invalid minutes entries:", "\n", XY[which(M.long > 0 & M.long < 60), ], "\n")
      stop("Range of valid minutes is from 0 to 60.")
    }
    if(any(DM.or.DMS.lat == 2)){
      if(any(S.long[which(DM.or.DMS.lat == 2)] < 0 & S.long[which(DM.or.DMS.lat == 2)] > 60)){
        cat("Invalid seconds entries:", "\n",
            XY[which(S.long[which(DM.or.DMS.lat == 2)] > 0 & S.long[which(DM.or.DMS.lat == 2)] < 60), ], "\n")
        stop("Range of valid seconds is from 0 to 60.")
      }
    }
    
    # Final checks that -90 <= decimal lat <= 90 and -180 <= decimal long <= 180, and then return the result.
    lat.res.check <- all(abs(DD.lat) <= 90)
    long.res.check <- all(abs(DD.long) <= 180)
    if(!lat.res.check & !long.res.check){
      stop("It appears an invalid answer has been calculated. Check for values just beyond the valid ranges of lat and long.")
    } else {
      return(cbind(DD.lat, DD.long))
    }
  }
#
# function for automatic mapping of bounding points for a given figure 
#



automapPts <- function(Splot,
                    fname.root="Splot",
                    boundFileName="SplotDot",
                    dir="./",
                    automap.method="mode"
                    ){
 

   
  # check if tif files where created
  if(dir == ""){
    d = dir("./")
  }else{
    d = dir(dir)
  }
  
  #if(dir=="./"){
    dot.loc = which(d == paste(boundFileName, ".tif", sep=""))
    fin.loc = which(d == paste(fname.root, ".tif", sep=""))
  #}else{
  #  dot.loc = which(d == paste(dir, boundFileName, ".tif", sep=""))
  #  fin.loc = which(d == paste(dir, fname.root, ".tif", sep=""))
  #}

  # if tifs where created continue 
  if( (length(dot.loc) != 0) | (length(fin.loc) != 0) ){

    require("rtiff")
    # reads tiff files 
    tif.dot = readTiff(paste(dir,boundFileName, ".tif", sep=""))
    tif.fin = readTiff(paste(dir,fname.root, ".tif", sep=""))

    # can compare based on three different channels, RGBobject 
    channels = c("blue", "red", "green")
    # keeps track of which channel attempting 
    idx = 1
    # keeps track if mapped correctly without errors 
    mapDif = FALSE
    # while not mapped correctly and more channels to try 
    while( (idx <= length(channels)) & (mapDif == FALSE) ){
      # attempt to get bounding limits 
      bounds = try(getBounds(channels[idx], tif.fin, tif.dot, automap.method=automap.method), silent=TRUE)
      idx = idx + 1
      # if bounds were retrived, mapping is correct
      if(class(bounds) != "try-error") mapDif=TRUE
      
    }
    # return bounding locations
    return(bounds)
    
  # if tifs where not created correctly 
  }else{ 
    tif.dot = paste(dir,boundFileName, ".tif", sep="")
    tif.fin = paste(dir,fname.root, ".tif", sep="")
    
    cat(paste("ERROR: could not map. \n       one of the required tif files is not found \n       Missing either ", tif.dot, " or ", tif.fin, "\n"))
    return(NA)
    
  }
  
}




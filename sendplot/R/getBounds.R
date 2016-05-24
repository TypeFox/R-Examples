
#
# function to try different channel comparisons to find difference
#    from original file to file with bounding points added to figure 
#

getBounds <- function(channelClr, tif.fin, tif.dot, automap.method="mode"){

  # blue channel
  if(channelClr == "blue"){

    # find where tifs differ
    temp.loc = which(tif.dot@blue != tif.fin@blue)
    temp = matrix(NA, nrow = dim(tif.dot@blue)[1], ncol =dim(tif.dot@blue)[2] )
    # convert to numeric matrix
    temp[which(tif.dot@blue == tif.fin@blue)] = 0
    temp[which(tif.dot@blue != tif.fin@blue)] = 1

  }

  # red channel
  if(channelClr == "red"){

    # find where tifs differ
    temp.loc = which(tif.dot@red != tif.fin@red)
    temp = matrix(NA, nrow = dim(tif.dot@red)[1], ncol =dim(tif.dot@red)[2] )
    # convert to numeric matrix
    temp[which(tif.dot@red == tif.fin@red)] = 0
    temp[which(tif.dot@red != tif.fin@red)] = 1

  }
  
  # green channel
  if(channelClr == "green"){

    # find where tifs differ
    temp.loc = which(tif.dot@green != tif.fin@green)
    temp = matrix(NA, nrow = dim(tif.dot@green)[1], ncol =dim(tif.dot@green)[2] )
    # convert to numeric matrix
    temp[which(tif.dot@green == tif.fin@green)] = 0
    temp[which(tif.dot@green != tif.fin@green)] = 1

  }
    
  # try mapping based on given indicies to find xlim/ylim for figure 
  bounds = try(mapMethod(automap.method, temp),silent=TRUE)

  # return bounding locations 
  return(bounds)
  
  
}

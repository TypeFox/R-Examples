# this function preprocesses panel data, difining time ranges, entities 
# creates maxtrix for x and y in wide format
# computes number of time (number of repeated measurements)
# computes number of entities (number of individual observations)

preprocess_data <- function(indnr, xv, yv, zv, vv) 
{  
  if (indnr == 2)
  {
    # arranging variable x and y in wide format (matrix, with rows representing entities and columns representing times)
    # computing number of measurement times and number of observation entities (e.g. number of countries)
    xwide <- as.matrix(xv)
    ywide <- as.matrix(yv)
    numTimes <- length(xwide[1, ]) 
    numEntities <- length(xwide[, 1]) 
    
    # reshaping variable x and variable y into a single column (long format ordered by time)
    tmpxw <- xwide[, -(ncol(xwide))]
    tmpyw <- ywide[, -(ncol(ywide))]
    allXt <- as.vector(tmpxw)
    allYt <- as.vector(tmpyw)
    
    # creating a long format vector for changes in variable x and variable y at each measurement point
    tmpchx <- rowDiffs(xwide)
    tmpchy <- rowDiffs(ywide)
    tiChXt <- as.vector(tmpchx)
    tiChYt <- as.vector(tmpchy)
    
    # computing means for the variables x and y 
    mx <- mean(allXt, na.rm=TRUE)
    my <- mean(allYt, na.rm=TRUE)
    
    # handling missing data, NAs are removed from the vectors
    idx <- which(is.na(allXt + allYt + tiChXt + tiChYt))   
    allX <- allXt[-idx]
    allY <- allYt[-idx]
    tiChX <- tiChXt[-idx]
    tiChY <- tiChYt[-idx]
    
    # scaling allX, allY, tiChX and tiChY
    xs <- allX/mx
    ys <- allY/my
    chXs <- tiChX/mx
    chYs <- tiChY/my
    
    # saving preprocessed data, if user wants to visually check them, uncomment please
    # save(xwide, ywide, numTimes, numEntities, allX, allY, tiChX, tiChY, mx, my, xs, ys, chXs, chYs, file = "var.RData")
    
    # defining what the function returns and asigning the returns to objects to be called
    return(list(xwide=xwide, ywide=ywide, numTimes=numTimes, numEntities=numEntities, allX=allX, 
         allY=allY, tiChX=tiChX, tiChY=tiChY, mx=mx, my=my, xs=xs, ys=ys, chXs=chXs, chYs=chYs))
  }
  
  if (indnr == 3)
  {
    # arranging variables x, y, z in wide format (matrix, with rows representing entities and columns representing times)
    # computing number of measurement times and number of observation entities (e.g. number of countries)
    xwide <- as.matrix(xv)
    ywide <- as.matrix(yv)
    zwide <- as.matrix(zv)
    numTimes <- length(xwide[1, ]) 
    numEntities <- length(xwide[, 1]) 
    
    # reshaping variables x, y, z into a single column (long format ordered by time)
    tmpxw <- xwide[, -(ncol(xwide))]
    tmpyw <- ywide[, -(ncol(ywide))]
    tmpzw <- zwide[, -(ncol(zwide))]
    allXt <- as.vector(tmpxw)
    allYt <- as.vector(tmpyw)
    allZt <- as.vector(tmpzw)
    
    # creating a long format vector for changes in variables x, y, z at each measurement point
    tmpchx <- rowDiffs(xwide)
    tmpchy <- rowDiffs(ywide)
    tmpchz <- rowDiffs(zwide)
    tiChXt <- as.vector(tmpchx)
    tiChYt <- as.vector(tmpchy)
    tiChZt <- as.vector(tmpchz)
    
    # computing means for the variables x, y, z 
    mx <- mean(allXt, na.rm=TRUE)
    my <- mean(allYt, na.rm=TRUE)
    mz <- mean(allZt, na.rm=TRUE)
    
    # handling missing data, NAs are removed from the vectors
    idx <- which(is.na(allXt + allYt + allZt + tiChXt + tiChYt + tiChZt))   
    allX <- allXt[-idx]
    allY <- allYt[-idx]
    allZ <- allZt[-idx]
    tiChX <- tiChXt[-idx]
    tiChY <- tiChYt[-idx]
    tiChZ <- tiChZt[-idx]
    
    # scaling allX, allY, allZ, tiChX, tiChY, tiChZ
    xs <- allX/mx
    ys <- allY/my
    zs <- allZ/mz
    chXs <- tiChX/mx
    chYs <- tiChY/my
    chZs <- tiChZ/mz
    
    # saving preprocessed data, if user wants to visually check them, uncomment please
    # save(xwide, ywide, zwide, numTimes, numEntities, allX, allY, allZ, tiChX, tiChY, tiChZ, mx, my, mz, xs, ys, zs, chXs, chYs, chZs, file = "var.RData")
    
    # defining what the function returns and asigning the returns to objects to be called
    return(list(xwide=xwide, ywide=ywide, zwide=zwide, numTimes=numTimes, numEntities=numEntities, allX=allX, 
         allY=allY, allZ=allZ, tiChX=tiChX, tiChY=tiChY, tiChZ=tiChZ, mx=mx, my=my, mz=mz, xs=xs, 
         ys=ys, zs=zs, chXs=chXs, chYs=chYs, chZs=chZs))
  }
  
  if (indnr == 4)
  {
    # arranging variables x, y, z in wide format (matrix, with rows representing entities and columns representing times)
    # computing number of measurement times and number of observation entities (e.g. number of countries)
    xwide <- as.matrix(xv)
    ywide <- as.matrix(yv)
    zwide <- as.matrix(zv)
    vwide <- as.matrix(vv)
    numTimes <- length(xwide[1, ]) 
    numEntities <- length(xwide[, 1]) 
    
    # reshaping variables x, y, z into a single column (long format ordered by time)
    tmpxw <- xwide[, -(ncol(xwide))]
    tmpyw <- ywide[, -(ncol(ywide))]
    tmpzw <- zwide[, -(ncol(zwide))]
    tmpvw <- vwide[, -(ncol(vwide))]
    allXt <- as.vector(tmpxw)
    allYt <- as.vector(tmpyw)
    allZt <- as.vector(tmpzw)
    allVt <- as.vector(tmpvw)
    
    # creating a long format vector for changes in variables x, y, z at each measurement point
    tmpchx <- rowDiffs(xwide)
    tmpchy <- rowDiffs(ywide)
    tmpchz <- rowDiffs(zwide)
    tmpchv <- rowDiffs(vwide)
    tiChXt <- as.vector(tmpchx)
    tiChYt <- as.vector(tmpchy)
    tiChZt <- as.vector(tmpchz)
    tiChVt <- as.vector(tmpchv)
    
    # computing means for the variables x, y, z 
    mx <- mean(allXt, na.rm=TRUE)
    my <- mean(allYt, na.rm=TRUE)
    mz <- mean(allZt, na.rm=TRUE)
    mv <- mean(allVt, na.rm=TRUE)
    
    # handling missing data, NAs are removed from the vectors
    idx <- which(is.na(allXt + allYt + allZt + allVt + tiChXt + tiChYt + tiChZt + tiChVt))   
    allX <- allXt[-idx]
    allY <- allYt[-idx]
    allZ <- allZt[-idx]
    allV <- allVt[-idx]
    tiChX <- tiChXt[-idx]
    tiChY <- tiChYt[-idx]
    tiChZ <- tiChZt[-idx]
    tiChV <- tiChVt[-idx]
    
    # scaling allX, allY, allZ, tiChX, tiChY, tiChZ
    xs <- allX/mx
    ys <- allY/my
    zs <- allZ/mz
    vs <- allV/mv
    chXs <- tiChX/mx
    chYs <- tiChY/my
    chZs <- tiChZ/mz
    chVs <- tiChV/mv
    
    # saving preprocessed data, if user wants to visually check them, uncomment please
    # save(xwide, ywide, zwide, vwide, numTimes, numEntities, allX, allY, allZ, allV, tiChX, tiChY, tiChZ, tiChV, mx, my, mz, mv, xs, ys, zs, vs, chXs, chYs, chZs, chVs, file = "var.RData")
    
    # defining what the function returns and asigning the returns to objects to be called
    return(list(xwide=xwide, ywide=ywide, zwide=zwide, vwide=vwide, numTimes=numTimes, numEntities=numEntities, 
         allX=allX, allY=allY, allZ=allZ, allV=allV, tiChX=tiChX, tiChY=tiChY, tiChZ=tiChZ, 
         tiChV=tiChV, mx=mx, my=my, mz=mz, mv=mv, xs=xs, ys=ys, zs=zs, vs=vs, chXs=chXs, chYs=chYs, 
         chZs=chZs, chVs=chVs))
  }  
}
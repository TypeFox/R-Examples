roctable <- function(indval, boolval, maskval=NULL, nthres=NULL, thres=NULL, NAval=0, progress=FALSE,
                     ones.bool=NULL, zeros.bool=NULL){ 
  
  # mask out nodata cells in the index and boolean vectors if a mask vector is given
  if(!is.null(maskval)){
    maskval[maskval == NAval] <- NA
    indval <- indval*maskval
    boolval <- boolval*maskval
  }
  
  if(is.null(ones.bool) | is.null(zeros.bool)){   
# extract total number of cells with ones and zeros in the boolean vector
boolvals <- boolval[!is.na(boolval)]
ones.bool <- sum(as.bit(boolvals))
zeros.bool <- length(boolvals) - ones.bool
  }

# generic function for crosstabing two boolean vectors 
func_logical3   <- function(v1,v2){
  r1  <- sum(v1 & v2)
  r2  <- sum(v1 & !v2)
  return(c(r1, r2))
}

# extract only no NA values from the boolean and index vectors
# length of (not NA) boolean and index vectors must be equal as the mask was applied previously
if(any(ls()=='boolvals')){
  boolval <- boolvals
} else{ 
  boolval <- boolval[!is.na(boolval)]
}

indval <- indval[!is.na(indval)]

if (length(boolval)!=length(indval)) stop('different NA values in input maps')

zeroIndVal <- indval*0
maxInd <- max(indval, na.rm=TRUE)

# define the thresholds vector
ifelse(!is.null(thres), minInd <- min(thres), minInd <- min(indval, na.rm=TRUE))

ifelse(!is.null(nthres), newThres <- (maxInd - minInd)/(nthres-2)*(0:(nthres-2)) + minInd, ifelse(!is.null(thres), 
                                                                                                  newThres <- thres, newThres <- unique(indval)))
newThres <- sort(newThres, decreasing=TRUE)

# create results data.frame
res <- cbind(newThres, "Hits"=0, "HitsRate"=0, "falseAlarms"=0, "falseAlarmsRate"=0)

# loop for reclassifying the index vector for each new threshold
# and then perform the crosstab with the boolean vector

for (j in 2:(nrow(res))){
  i <- newThres[j]
  zeroIndVal[which(indval > i)]  <- 1
  
  xb <- as.bit(zeroIndVal)
  yb <- as.bit(boolval)
  crsstb <- func_logical3(xb,yb)
  
  res[j,"Hits"] <- crsstb[1]
  res[j,"HitsRate"] <- crsstb[1]/ones.bool*100
  res[j,"falseAlarms"] <- crsstb[2]
  res[j,"falseAlarmsRate"] <- crsstb[2]/zeros.bool*100
  zeroIndVal <- indval*0
  
  if(progress){
    pb <- txtProgressBar(min = 0, max = 100, style = 3)
    setTxtProgressBar(pb, round(j/nrow(res)*100, 0))
  }
}


tocd <- as.data.frame(rbind(res, c(NA, ones.bool, 100, zeros.bool, 100)))

names(tocd) <- c("Threshold", "A", "HitsRate", "B", "falseAlarmsRate")

# calculate ROC table as shown in TOCfigure1.xlsx created by Pontius
tocd$Model1 <- tocd$HitsRate/100
tocd$falseAlarms1 <- tocd$falseAlarmsRate/100
tocd$Uniform <- tocd$falseAlarms1

return(list(tocd=tocd, minInd=minInd))
}
#' Function to impute the missing values in time series data
#'
#' @param dataIn as input time series data with missing values (NAs)
#' @param patch as missing data patches as returned with "missing_patch()" function
#' @return returns imputed data
#' @import PSF
#' @importFrom graphics plot
#' @importFrom graphics points
#' @importFrom stats na.omit
#' @export

#========================================================
# Function "heal_data()" starts here-----------------
#========================================================
heal_data <- function(dataIn, patch)
{
  blue <- NULL
  patch1 <- patch
  ###dataIn[is.na(dataIn)] <- 0
  #if(!(is.vector(dataIn)))
  #{
  #  dataIn <- dataIn[,1]
  #}
  #datax <- dataIn
  len <- length(dataIn)
  dataIn <- as.numeric(na.omit(dataIn))
  length(dataIn)
  datax <- dataIn
  while(nrow(patch) != 0)
  {
    max_val_patch <- which.max(patch$z)
    f_first <- patch[max_val_patch,][1]
    if(!is.vector(f_first))
    {
      f_first <- f_first[, 1]
    }

    f_last <- patch[max_val_patch,][2]
    if(!is.vector(f_last))
    {
      f_last <- f_last[, 1]
    }

    nextVal <- patch[max_val_patch,][3]
    if(!is.vector(nextVal))
    {
      nextVal <- nextVal[, 1]
    }

    f_mid <- (f_first + f_last)/2

    lens <- length(dataIn)

    if(f_first < 0.3*lens)
    #if(f_mid <= (len/3))
    {
      dataIn <- as.numeric(na.omit(dataIn))
      dataIn1 <- rev(dataIn)
      len1 <- length(dataIn1)
      f_first_rev <- len1 - f_last
      kOpt <- optimum_k(dataIn1[1:(f_first_rev)])
      wOpt <- optimum_w(dataIn1[1:(f_first_rev)], 1)$Optimum_W
      x2 <- pred_for_w(dataIn1[1:(f_first_rev)], wOpt, kOpt, nextVal)
      x3 <- rev(x2)
    }

    if((f_first >=  0.3*lens) && (f_last <= 0.7*lens))
    #if((f_mid > (len/3)) && (f_mid < (2*len/3)))
    {
      #x1 <- AUTO_PSF(dataIn[1:f_first-1],nextVal)$Predicted_Values
      dataIn <- as.numeric(na.omit(dataIn))
      kOpt <- optimum_k(dataIn[1:(f_first-1)])
      wOpt <- optimum_w(dataIn[1:(f_first-1)], 1)$Optimum_W
      x1 <- pred_for_w(dataIn[1:(f_first-1)], wOpt, kOpt, nextVal)

      dataIn1 <- rev(dataIn)
      dataIn1 <- as.numeric(na.omit(dataIn1))
      #x2 <- AUTO_PSF(dataIn1[1:f_first-1],nextVal)$Predicted_Values
      len1 <- length(dataIn1)
      f_first_rev <- len1 - f_last
      kOpt <- optimum_k(dataIn1[1:(f_first_rev)])
      wOpt <- optimum_w(dataIn1[1:(f_first_rev)], 1)$Optimum_W
      x2 <- pred_for_w(dataIn1[1:(f_first_rev)], wOpt, kOpt, nextVal)
      x2 <- rev(x2)

      x3 <- (x1+x2)/2

    }

    if(f_last > 0.7*lens)
    #if(f_mid >= (2*len/3))
    {
      dataIn <- as.numeric(na.omit(dataIn))
      kOpt <- optimum_k(dataIn[1:(f_first-1)])
      wOpt <- optimum_w(dataIn[1:(f_first-1)], 1)$Optimum_W
      x3 <- pred_for_w(dataIn[1:(f_first-1)], wOpt, kOpt, nextVal)
    }
  j <- 1
    ###for(i in f_first:f_last)
    ###{
      ###datax[i] <- x3[j]
      ###j <- j + 1
    ###}

    datax <- insert_patch(datax, (f_first - 1), x3)

    dataIn <- datax
    blue <- append(blue,f_first)
    blue <- append(blue,f_last)
    limit1 <- length(dataIn)
    plot(dataIn[1:f_first],type = "o", xlim = c(0,length(dataIn)),col=c("red"), xlab = "Time Series Data", ylab=NA)
    points(f_first:f_last+1, dataIn[f_first:f_last+1],col="blue", type ="o",pch=22)
    points(f_last+1:limit1,dataIn[f_last+1:limit1],col="red", type ="o")
    patch <- patch[-which.max(patch$z),]

    dataIn1 <- rev(dataIn)

  }
  blue <- sort(blue)
  blue <- append(blue, length(dataIn))
  plot(dataIn[1:blue[1]],type = "o", xlim = c(0,length(dataIn)),col=c("red"), xlab = "Time Series Data", ylab=NA)
  i=1
  #for(j in 1:length(blue)-2)

  while(nrow(patch1) != 0)

  {
    points(blue[i]:blue[i+1], dataIn[blue[i]:blue[i+1]],col="blue", type ="o",pch=22)
    points(blue[i+1]:blue[i+2], dataIn[blue[i+1]:blue[i+2]],col="red", type ="o")
    i <- i+2
    patch1 <- patch1[-1,]
  }

  #output <- list("Imputed Data" = dataIn, "Plot" = jk)
  #return(output)
  return(dataIn)
}
#========================================================
# Function "heal_data()" ends here-----------------
#========================================================




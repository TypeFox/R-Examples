#'Test the correlation between two climate variables
#'
#'Test the correlation between two climate variables across all considered climate
#'windows.
#'@param xvar The first climate variable of interest. Please specify the parent 
#'  environment and variable name (e.g. Climate$Temp).
#'@param xvar2 The second climate variable of interest. Please specify the parent 
#'  environment and variable name (e.g. Climate$Temp).
#'@param cdate The climate date variable (dd/mm/yyyy). Please specify the parent
#'  environment and variable name (e.g. Climate$Date).
#'@param bdate The biological date variable (dd/mm/yyyy). Please specify the 
#'  parent environment and variable name (e.g. Biol$Date).
#'@param furthest The furthest number of time intervals (set by cinterval) back
#'  from the cutoff date or biological record that will be included in the
#'  climate window search.
#'@param closest The closest number of time intervals (set by cinterval) back 
#'  from the cutoff date or biological record that will be included in the 
#'  climate window search.
#'@param stat The aggregate statistic used to analyse the climate data. Can 
#'  currently use basic R statistics (e.g. mean, min), as well as slope. 
#'  Additional aggregate statistics can be created using the format function(x) 
#'  (...). See FUN in \code{\link{apply}} for more detail.
#'@param stat2 Second aggregate statistic used to analyse climate data (xvar2). Can 
#'  currently use basic R statistics (e.g. mean, min), as well as slope. 
#'  Additional aggregate statistics can be created using the format function(x) 
#'  (...). See FUN in \code{\link{apply}} for more detail.
#'@param type fixed or variable, whether you wish the climate window to be variable
#'  (i.e. the number of days before each biological record is measured) or fixed
#'  (i.e. number of days before a set point in time).
#'@param cutoff.day,cutoff.month If type is "fixed", the day and month of the year
#'  from which the fixed window analysis will start.
#'@param cmissing TRUE or FALSE, determines what should be done if there are 
#'  missing climate data. If FALSE, the function will not run if missing climate
#'  data is encountered. If TRUE, any records affected by missing climate data 
#'  will be removed from climate window analysis.
#'@param cinterval The resolution at which climate window analysis will be 
#'  conducted. May be days ("day"), weeks ("week"), or months ("month"). Note the units 
#'  of parameters 'furthest' and 'closest' will differ depending on the choice 
#'  of cinterval
#'@return Will return a dataframe containing the correlation between the two
#'  climate variables.
#'@author Liam D. Bailey and Martijn van de Pol
#'@examples
#'\dontrun{
#'# Test correlation between temperature and rainfall in the MassClimate dataset.
#' 
#'data(Mass)
#'data(MassClimate)
#'
#'cross <- crosswin(xvar = list(Temp = MassClimate$Temp), 
#'                  xvar2 = list(Rain = MassClimate$Rain), 
#'                  cdate = MassClimate$Date, bdate = Mass$Date, 
#'                  furthest = 365, closest = 0,
#'                  stat = "mean", stat2 = "mean", type = "variable",
#'                  cmissing = FALSE, cinterval = "day")
#'                 
#'# View the output
#'head(cross)
#' 
#'# Plot the output
#'plotcor(cross, type = "C")
#' 
#'}
#' 
#'@export

#LAST EDITED: 18/02/2015
#EDITED BY: LIAM
#NOTES: Tidy up code

crosswin <- function(xvar, xvar2, cdate, bdate, furthest, closest, 
                     stat, stat2, type, cutoff.day, cutoff.month,
                     cinterval = "day", cmissing = FALSE){
  
  print("Initialising, please wait...")
  
  xvar  <- xvar[[1]]
  xvar2 <- xvar2[[1]]
  
  duration <- (furthest - closest) + 1
  maxmodno <- (duration * (duration + 1))/2 
  cont     <- convertdate(bdate = bdate, cdate = cdate, xvar = xvar, xvar2 = xvar2, 
                            cinterval = cinterval, type = type, 
                            cutoff.day = cutoff.day, cutoff.month = cutoff.month, cross = TRUE)   # create new climate dataframe with continuous daynumbers, leap days are not a problem
  modno    <- 1  #Create a model number variable that will count up during the loop#
  modlist  <- list()   # dataframes to store ouput
  cmatrix1 <- matrix(ncol = (duration), nrow = length(bdate))  # matrix that stores the weather data for variable or fixed windows
  cmatrix2 <- matrix(ncol = (duration), nrow = length(bdate))  # matrix that stores the weather data for variable or fixed windows
  
  for (i in 1:length(bdate)){
    for (j in closest:furthest){
      k <- j - closest + 1
      cmatrix1[i, k] <- cont$xvar[which(cont$cintn == cont$bintno[i] - j)]  #Create a matrix which contains the climate data from furthest to furthest from each biological record#
      cmatrix2[i, k] <- cont$xvar2[which(cont$cintno == cont$bintno[i] - j)]
    }
  }
  
  if (cmissing == FALSE && length(which(is.na(cmatrix1))) > 0){
    if(cinterval == "day"){
      .GlobalEnv$missing <- as.Date(cont$cintno[is.na(cont$xvar)], origin = min(as.Date(cdate, format = "%d/%m/%Y")) - 1)
    }
    if(cinterval == "month"){
      .GlobalEnv$missing <- c(paste("Month:", month(as.Date(cont$cintno[is.na(cont$xvar)], origin = min(as.Date(cdate, format = "%d/%m/%Y")) - 1)),
                                    "Year:", year(as.Date(cont$cintno[is.na(cont$xvar)], origin = min(as.Date(cdate, format = "%d/%m/%Y")) - 1))))
    }
    if(cinterval == "week"){
      .GlobalEnv$missing <- c(paste("Week:", month(as.Date(cont$cintno[is.na(cont$xvar)], origin = min(as.Date(cdate, format = "%d/%m/%Y")) - 1)),
                                    "Year:", year(as.Date(cont$cintno[is.na(cont$xvar)], origin = min(as.Date(cdate, format = "%d/%m/%Y")) - 1))))
    }
    stop(c("Climate data xvar should not contain NA values: ", length(.GlobalEnv$missing),
           " NA value(s) found. Please add missing climate data or set cmissing=TRUE.
           See object missing for all missing climate data"))
  }  
  
  if (cmissing == FALSE && length(which(is.na(cmatrix2))) > 0){
    if(cinterval == "day"){
      .GlobalEnv$missing <- as.Date(cont$cintno[is.na(cont$xvar2)], origin = min(as.Date(cdate, format = "%d/%m/%Y")) - 1)
    }
    if(cinterval == "month"){
      .GlobalEnv$missing <- c(paste("Month:", month(as.Date(cont$cintno[is.na(cont$xvar2)], origin = min(as.Date(cdate, format = "%d/%m/%Y")) - 1)),
                                    "Year:", year(as.Date(cont$cintno[is.na(cont$xvar2)], origin = min(as.Date(cdate, format = "%d/%m/%Y")) - 1))))
    }
    if(cinterval == "week"){
      .GlobalEnv$missing <- c(paste("Week:", month(as.Date(cont$cintno[is.na(cont$xvar2)], origin = min(as.Date(cdate, format = "%d/%m/%Y")) - 1)),
                                    "Year:", year(as.Date(cont$cintno[is.na(cont$xvar2)], origin = min(as.Date(cdate, format = "%d/%m/%Y")) - 1))))
    }
    stop(c("Climate data xvar2 should not contain NA values: ", length(.GlobalEnv$missing),
           " NA value(s) found. Please add missing climate data or set cmissing=TRUE.
           See object missing for all missing climate data"))
  }
  
  if (cmissing == TRUE){ 
    cmatrix1  <- cmatrix1[complete.cases(cmatrix1), ]
    cmatrix2  <- cmatrix2[complete.cases(cmatrix2), ]
  }
  
  climate1 <- matrix(ncol = 1, nrow = nrow(cmatrix1), 1)
  climate2 <- matrix(ncol = 1, nrow = nrow(cmatrix2), 1)
  
  pb <- txtProgressBar(min = 0, max = maxmodno, style = 3, char = "|")
  
  for (m in closest:furthest){
    for (n in 1:duration){
      if ( (m - n) >= (closest - 1)){  # do not use windows that overshoot the closest possible day in window   
        if (stat != "slope" || stat2 != "slope" || n > 1){
          windowopen  <- m - closest + 1
          windowclose <- windowopen - n + 1
          if (stat == "slope"){ 
            time <- seq(1, n, 1)
            climate1 <- apply(cmatrix1[, windowclose:windowopen], 1, FUN = function(x) coef(lm(x ~ time))[2])
            climate2 <- apply(cmatrix2[, windowclose:windowopen], 1, FUN = function(x) coef(lm(x ~ time))[2])
          } else { 
            if (n == 1) {
              climate1 <- cmatrix1[, windowclose:windowopen]
              climate2 <- cmatrix2[, windowclose:windowopen]
            } else {
              climate1 <- apply(cmatrix1[, windowclose:windowopen], 1, FUN = stat) 
              climate2 <- apply(cmatrix2[, windowclose:windowopen], 1, FUN = stat)
            }
            if (stat2 == "slope"){ 
              time     <- seq(1, n, 1)
              climate2 <- apply(cmatrix2[, windowclose:windowopen], 1, FUN = function(x) coef(lm(x ~ time))[2])
            } else { 
              if (n == 1) {
                climate2 <- cmatrix2[, windowclose:windowopen]
              } else {
                climate2 <- apply(cmatrix2[, windowclose:windowopen], 1, FUN = stat)
              }
            }
          }
          # Run the model
          modeloutput <- cor(climate1, climate2)
          # Add model parameters to list#
          modlist$cor[[modno]]         <- modeloutput
          modlist$WindowOpen[[modno]]  <- m
          modlist$WindowClose[[modno]] <- m - n + 1
          modno                        <- modno + 1 # Increase modno#
        }
      }
    }  
    #Fill progress bar
    setTxtProgressBar(pb, modno - 1)
  }
  modlist$Furthest    <- furthest
  modlist$Closest     <- closest
  modlist$Statistics  <- stat
  modlist$Statistics2 <- stat2
  modlist$Type        <- type
  
  if (type == "fixed"){
    modlist$Cutoff.day   <- cutoff.day
    modlist$Cutoff.month <- cutoff.month
  }
  return(as.data.frame(modlist))
}
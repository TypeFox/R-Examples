#'Fit a single climate window
#'
#'Fit a single climate window with a known start and end time.
#'@param xvar A list object containing all climate variables of interest. 
#'  Please specify the parent environment and variable name (e.g. Climate$Temp).
#'@param cdate The climate date variable (dd/mm/yyyy). Please specify the parent
#'  environment and variable name (e.g. Climate$Date).
#'@param bdate The biological date variable (dd/mm/yyyy). Please specify the 
#'  parent environment and variable name (e.g. Biol$Date).
#'@param baseline The baseline model structure used for testing correlation. 
#'  Currently known to support lm, glm, lmer and glmer objects.
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
#'@param func The functions used to fit the climate variable. Can be linear 
#'  ("lin"), quadratic ("quad"), cubic ("cub"), inverse ("inv") or log ("log").
#'@param type fixed or variable, whether you wish the climate window to be 
#'  variable (i.e. the number of days before each biological record is 
#'  measured) or fixed (i.e. number of days before a set point in time).
#'@param cutoff.day,cutoff.month If type is "fixed", the day and month of the 
#'  year from which the fixed window analysis will start.
#'@param cmissing TRUE or FALSE, determines what should be done if there are 
#'  missing climate data. If FALSE, the function will not run if missing 
#'  climate data is encountered. If TRUE, any records affected by missing 
#'  climate data will be removed from climate window analysis.
#'@param cinterval The resolution at which climate window analysis will be 
#'  conducted. May be days ("day"), weeks ("week"), or months ("month"). Note the units
#'  of parameters 'furthest' and 'closest' will differ depending on the choice
#'  of cinterval.
#'@param upper Cut-off values used to determine growing degree days or positive 
#'  climate thresholds (depending on parameter thresh). Note that when values
#'  of lower and upper are both provided, climatewin will instead calculate an 
#'  optimal climate zone.
#'@param lower Cut-off values used to determine chill days or negative 
#'  climate thresholds (depending on parameter thresh). Note that when values
#'  of lower and upper are both provided, climatewin will instead calculate an 
#'  optimal climate zone.
#'@param thresh TRUE or FALSE. Determines whether to use values of upper and
#'  lower to calculate binary climate data (thresh = TRUE), or to use for
#'  growing degree days (thresh = FALSE).
#'@param centre Variable used for mean centring (e.g. Year, Site, Individual).
#'  Please specify the parent environment and variable name (e.g. Biol$Year).
#'@return Will return a list containing two objects:
#'  
#'  \itemize{
#'  \item BestModel, a model object of the fitted climate window
#'  model.
#'  
#'  \item BestModelData, a dataframe with the biological and climate data
#'  used to fit the climate window model.}
#'  
#'@author Liam D. Bailey and Martijn van de Pol
#'@examples
#'\dontrun{
#'# Fit a known climate window to the datasets Mass and MassClimate
#'
#'data(Mass)
#'data(MassClimate)
#'
#'# Test for a fixed climate window, starting from 20th May
#'# Fit a climate window starting 72 days ago and ending 15 days ago
#'# Fit a linear term for the mean climate
#'# Fit climate windows at the resolution of days
#'
#'single <- singlewin(xvar = list(Temp = MassClimate$Temp), 
#'                    cdate = MassClimate$Date, bdate = Mass$Date,
#'                    baseline = lm(Mass ~ 1, data = Mass), 
#'                    furthest = 72, closest = 15,
#'                    stat = "mean", func = "lin",
#'                    type = "fixed", cutoff.day = 20, cutoff.month = 5,
#'                    cmissing = FALSE, cinterval = "day")
#'                
#'##View data##
#'single$BestModel
#'head(single$BestModelData)
#'}
#'
#'@importFrom MuMIn AICc
#'@importFrom lubridate year
#'@importFrom lubridate month
#'@export

singlewin <- function(xvar, cdate, bdate, baseline, 
                      furthest, closest, stat, func, 
                      type, cutoff.day, cutoff.month, 
                      cmissing = FALSE, cinterval = "day",
                      upper = NA, lower = NA, thresh = FALSE,
                      centre = NULL){
  
  xvar = xvar[[1]]
  
  if(stat == "slope" & func == "log" || stat == "slope" & func == "inv"){
    stop("stat = slope cannot be used with func = LOG or I as negative values may be present")
  }
  
  duration  <- (furthest - closest) + 1
  
  bdate  <- as.Date(bdate, format = "%d/%m/%Y") # Convert the date variables into the R date format
  cdate2 <- seq(min(as.Date(cdate, format = "%d/%m/%Y")), max(as.Date(cdate, format = "%d/%m/%Y")), "days") # Convert the date variables into the R date format
  cdate  <- as.Date(cdate, format = "%d/%m/%Y")
  
  if(min(cdate) > min(bdate)){
    stop("Climate data does not cover all years of biological data. Please increase range of climate data")
  }
  
  xvar <- xvar[match(cdate2, cdate)]
  
  cintno     <- as.numeric(cdate2) - min(as.numeric(cdate2)) + 1   # atrribute daynumbers for both datafiles with first date in CLimateData set to cintno 1
  realbintno <- as.numeric(bdate) - min(as.numeric(cdate2)) + 1
  
  if (length(cintno) != length(unique(cintno))){
    stop ("There are duplicate dayrecords in climate data")
  }
  
  if (cinterval != "day" && cinterval != "week" && cinterval != "month"){
    stop("cinterval should be either day, week or month")
  }
  
  if(cinterval == "day"){  
    if(type == "fixed"){   
      bintno            <- as.numeric(as.Date(paste(cutoff.day, cutoff.month, year(bdate), sep = "-"), format = "%d-%m-%Y")) - min(as.numeric(cdate2)) + 1 
      wrongyear         <- which(bintno < realbintno)
      bintno[wrongyear] <- (as.numeric(as.Date(paste(cutoff.day, cutoff.month, (year(bdate[wrongyear]) + 1), sep = "-"), format = "%d-%m-%Y")) - min(as.numeric(cdate2)) + 1)
    } else {
      bintno <- realbintno
    }
  } else if (cinterval == "week"){
    cintno     <- ceiling((as.numeric(cdate2) - min(as.numeric(cdate2)) + 1) / 7)   # atrribute weeknumbers for both datafiles with first week in CLimateData set to cintno 1
    realbintno <- ceiling((as.numeric(bdate) - min(as.numeric(cdate2)) + 1) / 7)
    newclim    <- data.frame("cintno" = cintno, "xvar" = xvar)
    newclim2   <- melt(newclim, id = "cintno")
    newclim3   <- cast(newclim2, cintno~variable, mean)
    cintno     <- newclim3$cintno
    xvar       <- newclim3$xvar
    if (type == "fixed"){ 
      bintno            <- ceiling((as.numeric(as.Date(paste(cutoff.day, cutoff.month, year(bdate), sep = "-"), format = "%d-%m-%Y")) - min(as.numeric(cdate2)) + 1) / 7) 
      wrongyear         <- which(bintno < realbintno)
      bintno[wrongyear] <- ceiling((as.numeric(as.Date(paste(cutoff.day, cutoff.month, (year(bdate[wrongyear]) + 1), sep = "-"), format = "%d-%m-%Y")) - min(as.numeric(cdate2)) + 1) / 7)
    } else {
      bintno <- realbintno
    }
  } else if (cinterval == "month"){ 
    cmonth     <- month(cdate2)
    cyear      <- year(cdate2) - min(year(cdate2))
    cintno     <- cmonth + 12 * cyear
    realbintno <- month(bdate) + 12 * (year(bdate) - min(year(cdate2)))
    newclim    <- data.frame("cintno" = cintno, "xvar" = xvar)
    newclim2   <- melt(newclim, id = "cintno")
    newclim3   <- cast(newclim2, cintno ~ variable, mean)
    cintno     <- newclim3$cintno
    xvar       <- newclim3$xvar
    if (type == "fixed"){ 
      bintno            <- cutoff.month + 12 * (year(bdate) - min(year(cdate2)))
      wrongyear         <- which(bintno < realbintno)
      bintno[wrongyear] <- cutoff.month + 12 * (year(bdate[wrongyear]) + 1 - min(year(cdate2)))
    } else {
      bintno <- realbintno
    }
  }
  
  if(cinterval == "day"){
    if((min(bintno) - furthest) < min(cintno)){
      stop("You do not have enough climate data to search that far back. Please adjust the value of furthest or add additional climate data.")
    }
  }

  if(cinterval == "week"){
    if((min(bintno) - furthest * 7) < min(cintno)){
      stop("You do not have enough climate data to search that far back. Please adjust the value of furthest or add additional climate data.")
    }
  }
  
  if(cinterval == "month"){
    if((as.numeric(min(as.Date(bdate, format = "%d/%m/%Y")) - months(furthest)) - (as.numeric(min(as.Date(cdate, format = "%d/%m/%Y"))))) <= 0){
      stop("You do not have enough climate data to search that far back. Please adjust the value of furthest or add additional climate data.")
    }
  }
  
  if(max(bintno) > max(cintno)){
    if(type == "fixed"){
      stop("You need more recent biological data. This error may be caused by your choice of cutoff.day/cutoff.month")
    } else {
      stop("You need more recent biological data")
    }
  }
  
  baseline  <- update(baseline, .~.)
  nullmodel <- AICc(baseline)  
  modlist   <- list()   # dataframes to store ouput
  cmatrix   <- matrix(ncol = (duration), nrow = length(bdate))
  
  modeldat      <- model.frame(baseline)
  modeldat$yvar <- modeldat[, 1]
  
  if(is.null(centre) == FALSE){
    func = "centre"
  }
  
  if(length(modeldat$yvar) != length(bdate)){
    stop("NA values present in biological response. Please remove NA values")
  }
  
  if(is.na(upper) == FALSE && is.na(lower) == TRUE){
    if(thresh == TRUE){
      xvar <- ifelse(xvar > upper, 1, 0)
    } else {
      xvar <- ifelse(xvar > upper, xvar, 0)
    }
  }
  
  
  if(is.na(lower) == FALSE && is.na(upper) == TRUE){
    if(thresh == TRUE){
      xvar <- ifelse(xvar < lower, 1, 0)
    } else {
      xvar <- ifelse(xvar < lower, xvar, 0)
    }
  }
  
  if(is.na(lower) == FALSE && is.na(upper) == FALSE){
    if(thresh == TRUE){
      xvar <- ifelse(xvar > lower & xvar < upper, 1, 0)
    } else {
      xvar <- ifelse(xvar > lower & xvar < upper, xvar - lower, 0)
    } 
  }  
  
  for (i in 1:length(bdate)){
    for (j in closest:furthest){
      k <- j - closest + 1
      cmatrix[i, k] <- xvar[which(cintno == bintno[i] - j)]   #Create a matrix which contains the climate data from furthest to furthest from each biological record#    
    }
  }
  
  if (cmissing == FALSE && length(which(is.na(cmatrix))) > 0){
    if(cinterval == "day"){
      .GlobalEnv$missing <- as.Date(cintno[is.na(xvar)], origin = min(as.Date(cdate, format = "%d/%m/%Y")) - 1)
    }
    if(cinterval == "month"){
      .GlobalEnv$missing <- c(paste("Month:", month(as.Date(cintno[is.na(xvar)], origin = min(as.Date(cdate, format = "%d/%m/%Y")) - 1)),
                                    "Year:", year(as.Date(cintno[is.na(xvar)], origin = min(as.Date(cdate, format = "%d/%m/%Y")) - 1))))
    }
    if(cinterval == "week"){
      .GlobalEnv$missing <- c(paste("Week:", month(as.Date(cintno[is.na(xvar)], origin = min(as.Date(cdate, format = "%d/%m/%Y")) - 1)),
                                    "Year:", year(as.Date(cintno[is.na(xvar)], origin = min(as.Date(cdate, format = "%d/%m/%Y")) - 1))))
    }
    stop(c("Climate data should not contain NA values: ", length(.GlobalEnv$missing),
           " NA value(s) found. Please add missing climate data or set cmissing = TRUE.
           See object missing for all missing climate data"))
  }
  
  if (cmissing == TRUE && length(which(is.na(cmatrix))) > 0){
    modeldat      <- modeldat[complete.cases(cmatrix), ]
    baseline      <- update(baseline, yvar~., data = modeldat)
    cmatrix       <- cmatrix[complete.cases(cmatrix), ]
  }
  
  modeldat           <- model.frame(baseline)
  modeldat$yvar      <- modeldat[, 1]
  modeldat$climate   <- matrix(ncol = 1, nrow = nrow(cmatrix), seq(from = 1, to = nrow(cmatrix), by = 1))
  
  if(is.null(weights(baseline)) == FALSE){
    if(class(baseline)[1] == "glm" & sum(weights(baseline)) == nrow(model.frame(baseline)) || class(baseline)[1] == "lmerMod" & sum(weights(baseline)) == nrow(model.frame(baseline))){
    } else {
      modeldat$modweights <- weights(baseline)
      baseline <- update(baseline, .~., weights = modeldat$modweights, data = modeldat)
    }
  }
  
  if (func == "lin"){
    modeloutput <- update(baseline, yvar~. + climate, data = modeldat)
  } else if (func == "quad") {
    modeloutput <- update(baseline, yvar~. + climate + I(climate ^ 2), data = modeldat)
  } else if (func == "cub") {
    modeloutput <- update(baseline, yvar~. + climate + I(climate ^ 2) + I(climate ^ 3), data = modeldat)
  } else if (func == "log") {
    modeloutput <- update(baseline, yvar~. + log(climate), data = modeldat)
  } else if (func == "inv") {
    modeloutput <- update (baseline, yvar~. + I(climate ^ -1), data = modeldat)
  } else if (func == "centre"){
    modeldat$wgdev  <- matrix(ncol = 1, nrow = nrow(cmatrix), seq(from = 1, to = nrow(cmatrix), by = 1))
    modeldat$wgmean <- matrix(ncol = 1, nrow = nrow(cmatrix), seq(from = 1, to = nrow(cmatrix), by = 1))
    modeloutput <- update (baseline, yvar ~. + wgdev + wgmean, data = modeldat)
  } else {
    print("Define func")
  }
  
  #CREATE A FOR LOOP TO FIT DIFFERENT CLIMATE WINDOWS#
  m     <- closest
  n     <- duration

  #Save the best model output
  if (stat == "slope"){
    time             <- seq(1, n, 1)
    modeldat$climate <- apply(cmatrix, 1, FUN = function(x) coef(lm(x ~ time))[2])
  } else {
    ifelse(n == 1, modeldat$climate <- cmatrix, modeldat$climate <- apply(cmatrix, 1, FUN = stat))
  }
  if(is.null(centre) == FALSE){
    modeldat$WGdev  <- wgdev(modeldat$climate, centre)
    modeldat$WGmean <- wgmean(modeldat$climate, centre)
    LocalBestModel  <- update(modeloutput, .~., data = modeldat)
  } else {
    LocalBestModel <- update(modeloutput, .~.)
  }
  LocalData           <- model.frame(LocalBestModel)
  return(list(BestModel = LocalBestModel, BestModelData = LocalData))
}
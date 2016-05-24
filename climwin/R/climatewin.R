#'Test for a climate windows in data.
#'
#'Finds the time period when a biological variable is most strongly affected 
#'by climate. Note that climate data and biological data should be loaded as 
#'two seperate objects. Both objects should contain a date column to designate
#'when the data were recorded (dd/mm/yyyy) and a response variable column.
#'
#'Note that climatewin allows you to test multiple possible parameters with the
#'same code (e.g. func, stat, xvar). See examples for more detail.
#'@param xvar A list object containing all climate variables of interest. 
#'  Please specify the parent environment and variable name (e.g. climate$Temp).
#'@param cdate The climate date variable (dd/mm/yyyy). Please specify the 
#'  parent environment and variable name (e.g. climate$Date).
#'@param bdate The biological date variable (dd/mm/yyyy). Please specify the 
#'  parent environment and variable name (e.g. Biol$Date).
#'@param baseline The baseline model structure used for model testing. 
#'  Currently known to support lm, glm, lmer and glmer objects.
#'@param furthest The furthest number of time intervals (set by cinterval) 
#'  back from the cutoff date or biological record that will be included in
#'  the climate window search.
#'@param closest The closest number of time intervals (set by cinterval) back 
#'  from the cutoff date or biological record that will be included in the
#'  climate window search.
#'@param stat The aggregate statistics used to analyse the climate data. Can 
#'  currently use basic R statistics (e.g. mean, min), as well as slope. 
#'  Additional aggregate statistics can be created using the format 
#'  function(x) (...). See FUN in \code{\link{apply}} for more detail.
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
#'@param cvk The number of folds used for k-fold cross validation. By default
#'  this value is set to 0, so no cross validation occurs. Value should be a
#'  minimum of 2 for cross validation to occur.
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
#'@return Will return a list with an output for each tested set of climate
#'  window parameters. Each list item contains three objects:
#'  
#'  \itemize{
#'  \item BestModel, a model object. The strongest climate window model based on AICc. 
#'  \item BestModelData, a dataframe with the data used to fit the strongest 
#'  climate window model.
#'  \item Dataset, a dataframe with information on all fitted climate windows. 
#'  Ordered using deltaAICc, with most negative deltaAICc values first. 
#'  See \code{\link{MassOutput}} as an example.}
#'  
#'  In addition, the returned list includes an object 'combos', a summary of all
#'  tested sets of climate window parameters. 
#'@author Liam D. Bailey and Martijn van de Pol
#'@importFrom MuMIn AICc
#'@import lme4
#'@import stats
#'@import utils
#'@import graphics
#'@importFrom lubridate weeks  
#'@examples
#'\dontrun{
#'##EXAMPLE 1## 
#'  
#'# Test both a linear and quadratic variable climate window using datasets "Offspring"
#'# and "OffspringClimate".
#'
#'# Load data.
#'
#'data(Offspring) 
#'data(OffspringClimate)
#'
#'# Test both linear and quadratic functions with climate variable temperature
#'
#'OffspringWin <- climatewin(xvar = list(Temp = OffspringClimate$Temperature), 
#'                           cdate = OffspringClimate$Date, 
#'                           bdate = Offspring$Date, 
#'                           baseline = glm(Offspring ~ 1, data = Offspring, family = poisson),
#'                           furthest = 150, closest = 0, 
#'                           type = "variables", stat = "mean", 
#'                           func = c("lin", "quad"), cmissing = FALSE, cinterval = "day")
#'
#'# Examine tested combinations
#'  
#'OffspringWin$combos
#'      
#'# View output for func = "lin"
#'  
#'head(OffspringWin[[1]]$Dataset) 
#'summary(OffspringWin[[1]]$BestModel)
#'  
#'# View output for func = "quad"
#'  
#'head(OffspringWin[[2]]$Dataset)
#'summary(OffspringWin[[2]]$BestModel)
#'  
#'##EXAMPLE 2##
#'  
#'# Test for a fixed climate window with both 'mean' and 'max' aggregate statistics
#'# using datasets 'Mass' and 'MassClimate'.
#'  
#'# Load data.
#'  
#'data(Mass)
#'data(MassClimate)
#'  
#'# Test a fixed window, starting 20 May (cutoff.month = 5, cutoff.day = 20)
#'# Test for climate windows between 100 and 0 days ago (furthest = 100, closest = 0)
#'# Test both mean and max aggregate statistics (stat = c("mean", "max"))
#'# Fit a linear term (func = "lin")
#'# Test at the resolution of days (cinterval = "day")
#'  
#'MassWin <- climatewin(xvar = list(Temp = MassClimate$Temp), cdate = MassClimate$Date, 
#'                      bdate = Mass$Date, baseline = lm(Mass ~ 1, data = Mass),
#'                      furthest = 100, closest = 0,
#'                      stat = c("mean", "max"), func = "lin",
#'                      type = "fixed", cutoff.day = 20, cutoff.month = 5,
#'                      cmissing = FALSE, cinterval = "day")
#'                        
#'# Examine tested combinations
#'  
#'MassWin$combos                      
#'  
#'# View output for mean temperature
#'  
#'head(MassWin[[1]]$Dataset)
#'summary(MassWin[[1]]$BestModel)
#'  
#'# View output for max temperature
#'  
#'head(MassWin[[2]]$Dataset)
#'summary(MassWin[[2]]$BestModel)
#'  
#'}
#'  
#'@export

climatewin <- function(xvar, cdate, bdate, baseline, furthest, closest, 
                       type, cutoff.day, cutoff.month, stat = "mean", func = "lin",
                       cmissing = FALSE, cinterval = "day", cvk = 0,
                       upper = NA, lower = NA, thresh = FALSE, centre = NULL){
  
  #Make xvar a list where the name of list object is the climate variable (e.g. Rain, Temp)
  if (is.list(xvar) == FALSE){
    stop("xvar should be an object of type list")
  }

  if (is.null(names(xvar)) == TRUE){
    numbers <- seq(1, length(xvar), 1)
    for (xname in 1:length(xvar)){
      names(xvar)[xname] = paste("climate", numbers[xname])
    }
  }
  
  if (is.na(upper) == FALSE && is.na(lower) == FALSE){
    combos       <- expand.grid(list(upper = upper, lower = lower))
    combos       <- combos[which(combos$upper >= combos$lower), ]
    allcombos    <- expand.grid(list(climate = names(xvar), type = type, stat = stat, func = func, gg = c(1:nrow(combos)), thresh = thresh))
    allcombos    <- cbind(allcombos, combos[allcombos$gg, ], deparse.level = 2)
    threshlevel  <- "two"
    allcombos$gg <- NULL
  } else if (is.na(upper) == FALSE && is.na(lower) == TRUE){
    allcombos   <- expand.grid(list(climate = names(xvar), type = type, stat = stat, func = func, upper = upper, thresh = thresh))
    threshlevel <- "upper"
  } else if (is.na(upper) == TRUE && is.na(lower) == FALSE){
    allcombos   <- expand.grid(list(climate = names(xvar), type = type, stat = stat, func = func, lower = lower, thresh = thresh))
    threshlevel <- "lower"
  } else if (is.na(upper) == TRUE && is.na(lower) == TRUE){
    allcombos   <- expand.grid(list(climate = names(xvar), type = type, stat = stat, func = func))
    threshlevel <- "none"
  }
  
  rownames(allcombos) <- seq(1, nrow(allcombos), 1)
  print("All combinations to be tested...")
  print(allcombos)
  
  combined <- list()
  for (combo in 1:nrow(allcombos)){
    runs <- basewin(xvar = xvar[[paste(allcombos[combo, 1])]], cdate = cdate, bdate = bdate, baseline = baseline,
                    furthest = furthest, closest = closest, type = paste(allcombos[combo, 2]), cutoff.day = cutoff.day,
                    cutoff.month = cutoff.month, stat = paste(allcombos[combo, 3]), func = paste(allcombos[combo, 4]),
                    cmissing = cmissing, cinterval = cinterval, cvk = cvk, 
                    upper = ifelse(threshlevel == "two" || threshlevel == "upper", allcombos$upper[combo], NA),
                    lower = ifelse(threshlevel == "two" || threshlevel == "lower", allcombos$lower[combo], NA),
                    thresh = paste(allcombos$thresh[combo]), centre = centre)
    combined[[combo]] <- runs
  }
  combined <- c(combined, combos = list(allcombos))
  return(combined)
}

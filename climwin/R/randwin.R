#'Climate window analysis for randomised data
#'
#'Will randomise biological data and carry out a climate window analysis. Used
#'to help determine the chance of obtaining an observed result at random.
#'@param repeats The number of times that data will be randomised and analysed 
#'  for climate windows.
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
#'@return Will return a dataframe containing information on all fitted climate
#'  windows. See \code{\link{MassRand}} as an example.
#'@author Liam D. Bailey and Martijn van de Pol
#' @examples
#' \dontrun{
#'# Test climate windows for random data using Mass dataset
#' 
#'data(Mass)
#'data(MassClimate)
#' 
#'# Randomise data twice
#'# Note all other parameters are fitted in the same way as the climatewin function.
#' 
#'rand <- randwin(repeats = 2, xvar = list(Temp = MassClimate$Temp), 
#'                cdate = MassClimate$Date, bdate = Mass$Date,
#'                baseline = lm(Mass ~ 1, data = Mass), 
#'                furthest = 100, closest = 0,
#'                stat = "mean", func = "lin", type = "fixed", 
#'                cutoff.day = 20, cutoff.month = 5,
#'                cmissing = FALSE, cinterval = "day")
#'                 
#'# View output #
#' 
#'head(rand)    
#'        }
#'
#'@export

randwin <- function(repeats = 1, xvar, cdate, bdate, baseline, 
                    furthest, closest, stat,  
                    func, type, cutoff.day, cutoff.month,
                    cmissing = FALSE, cinterval = "day",
                    upper = NA, lower = NA, thresh = FALSE, centre = NULL){
  
  xvar = xvar[[1]]
  
  for (r in 1:repeats){
    print (c("randomization number ", r))
    bdateNew        <- sample(bdate)
    outputrep <- basewin(xvar = xvar, cdate = cdate, bdate = bdateNew, 
                         baseline = baseline, furthest = furthest,
                         closest = closest, stat = stat, 
                         func = func, type = type,
                         cutoff.day = cutoff.day, cutoff.month = cutoff.month,
                         nrandom = repeats, cmissing = cmissing, cinterval = cinterval,
                         upper = upper, lower = lower, thresh = thresh, centre = centre)
    
    outputrep$Repeat <- r
    
    if(r == 1){ 
      outputrand <- outputrep 
    } else { 
      outputrand <- rbind(outputrand, outputrep)
    }
    rm(outputrep)
  }
  return(as.data.frame(outputrand))
} 
#'Find a weighted climate window
#'
#'Finds the best weighted average of a weather variable over a period that 
#'correlates most strongly with a biological variable. Uses weibull or 
#'Generalised Extreme Value (GEV) distribution. See references for a full 
#'description.
#'
#'@param xvar A list object containing all climate variables of interest. 
#'  Please specify the parent environment and variable name (e.g. Climate$Temp).
#'@param cdate The climate date variable. Please specify the parent environment 
#'  and variable name (e.g. Climate$Date).
#'@param bdate The biological date variable. Please specify the parent 
#'  environment and variable name (e.g. Biol$Date).
#'@param baseline The baseline model structure used for testing correlation. 
#'  Currently known to support lm, lme, glm and glmer objects.
#'@param furthest The furthest number of time intervals (set by cinterval) 
#'  back from the cutoff date or biological record that will be included in
#'  the climate window search.
#'@param closest The closest number of time intervals (set by cinterval) back 
#'  from the cutoff date or biological record that will be included in the
#'  climate window search.
#'@param func The function used to fit the climate variable in the model. Can be
#'  linear ("lin"), quadratic ("quad"), cubic ("cub"), inverse ("inv") or log ("log").
#'@param type fixed or variable, whether you wish the climate window to be variable
#'  (i.e. the number of days before each biological record is measured) or fixed
#'  (i.e. number of days before a set point in time).
#'@param cutoff.day,cutoff.month If type is "fixed", the day and month of the year
#'  from which the fixed window analysis will start.
#'@param weightfunc The distribution to be used for optimisation. Can be 
#'  either a Weibull ("W") or Generalised Extreme Value distribution ("G").
#'@param cinterval The resolution at which the climate window analysis will be 
#'  conducted. May be days ("day"), weeks ("week"), or months ("month"). Note the units 
#'  of parameters 'furthest' and 'closest' will differ depending on the choice 
#'  of cinterval.
#'@param par Shape, scale and location parameters of the Weibull of GEV weight 
#'  function used as start weight function. For Weibull : Shape and scale 
#'  parameters must be greater than 0, while location parameter must be less 
#'  than or equal to 0. For GEV : Scale parameter must be greater than 0.
#'@param control Parameters used to determine step size for the optimisation 
#'  function. Please see \code{\link{optim}} for more detail.
#'@param method The method used for the optimisation function. Please see 
#'  \code{\link{optim}} for more detail.
#'@references van de Pol & Cockburn 2011 Am Nat 177(5):698-707 (doi: 
#'  10.1086/659101) "Identifying the critical climatic time window that affects 
#'  trait expression"
#'@return Produces a constantly updating grid of plots as the optimisation 
#'  function is running. 
#'  \itemize{ 
#'  \item Right panel from top to bottom: The
#'  three parameters (shape, scale and location) determining the weight 
#'  function.
#'  
#'  \item Left top panel: The resulting weight function.
#'  
#'  \item Right middle panel: The delta AICc compared to the baseline model.
#'  
#'  \item Right bottom panel: The weighted mean of climate for the current
#'  weight function. }
#'  
#'  Also returns a list containing three objects: \itemize{ 
#'  \item BestModel, a model object. The best weighted window model determined
#'  by AICc.
#'  
#'  \item BestModelData, a dataframe. Biological and climate data used to fit
#'  the best weighted window model.
#'  
#'  \item WeightedOutput, a list. Parameter values for the best weighted window.
#'  }
#'@author Martijn van de Pol and Liam D. Bailey
#'  @examples
#'  \dontrun{
#'  
#'# Test for a weighted average over a fixed climate window 
#'# using datasets 'Offspring' and 'OffspringClimate'
#'  
#'# N.B. THIS EXAMPLE MAY TAKE A MOMENT TO CONVERGE ON THE BEST MODEL.
#'  
#'# Load data
#'  
#'data(Offspring)
#'data(OffspringClimate)
#'  
#'# Test for climate windows between 365 and 0 days ago (furthest=365, closest=0)
#'# Fit a quadratic term for the mean weighted climate (func="quad")
#'# in a Poisson regression (offspring number ranges 0-3)
#'# Test a variable window (type = "fixed")
#'# Test at the resolution of days (cinterval="day")
#'# Uses a Weibull weight function (weightfunc="week")
#'  
#'weight <- weightwin(xvar = list(Temp = OffspringClimate$Temperature), 
#'                    cdate = OffspringClimate$Date, 
#'                    bdate = Offspring$Date, 
#'                    baseline = glm(Offspring ~ 1, family = poisson, data = Offspring), 
#'                    furthest = 365, closest = 0, func = "quad", 
#'                    type = "variable", weightfunc = "W", cinterval = "day", 
#'                    par = c(3, 0.2, 0), control = list(ndeps = c(0.01, 0.01, 0.01)), 
#'                    method = "L-BFGS-B") 
#'  
#'# View output
#'  
#'head(weight[[3]])
#'summary(weight[[1]])
#'head(weight[[2]])
#'  }
#'
#'@importFrom evd dgev
#'@export

weightwin <- function(xvar, cdate, bdate, baseline, furthest, closest, 
                      func = "lin", type = "fixed", cutoff.day, cutoff.month, 
                      weightfunc = "W", cinterval = "day",
                      par = c(3, 0.2, 0), control = list(ndeps = c(0.01, 0.01, 0.01)), 
                      method = "L-BFGS-B"){
  
  xvar = xvar[[1]]
  
  funcenv                 <- environment()
  cont                    <- convertdate(bdate = bdate, cdate = cdate, xvar = xvar, 
                                         cinterval = cinterval, type = type, 
                                         cutoff.day = cutoff.day, cutoff.month = cutoff.month )   # create new climate dataframe with continuous daynumbers, leap days are not a problem 

  modno        <- 1
  DAICc        <- list()
  par_shape    <- list()
  par_scale    <- list()
  par_location <- list()
  duration     <- (furthest - closest) + 1
  cmatrix      <- matrix(ncol = (duration), nrow = length(bdate))
  baseline     <- update(baseline, .~.)
  nullmodel    <- AICc(baseline)
  
  for (i in 1:length(bdate)){
    for (j in closest:furthest){
      k <- j - closest + 1
      cmatrix[i, k] <- xvar[match(cont$bintno[i] - j, cont$cintno)]   #Create a matrix which contains the climate data from furthest to furthest from each biological record#    
    }
  }
  
  funcenv$modeldat         <- model.frame(baseline)
  funcenv$modeldat$climate <- matrix(ncol = 1, nrow = nrow(cmatrix), seq(from = 1, to = nrow(cmatrix), by = 1))
  
  if (func == "lin"){
    modeloutput <- update(baseline, .~. + climate, data = modeldat)
  } else if (func == "quad") {
    modeloutput <- update(baseline, .~. + climate + I(climate ^ 2), data = modeldat)
  } else if (func == "cub") {
    modeloutput <- update(baseline, .~. + climate + I(climate ^ 2) + I(climate ^ 3), data = modeldat)
  } else if (func == "log") {
    modeloutput <- update(baseline, .~. + log(climate), data = modeldat)
  } else if (func == "inv") {
    modeloutput <- update (baseline, .~. + I(climate ^ -1), data = modeldat)
  } else {
    print("Define func")
  }
  
  # now run one of two optimization functions
  if (weightfunc == "W"){
    if (par[1] <= 0){
      stop("Weibull shape parameter should be >0")
    }
    if (par[2] <= 0){
      stop("Weibull scale parameter should be >0")
    }
    if (par[3] > 0){
      stop("Weibull location parameter should be <=0")
    }
    j      <- seq(1:duration) / duration
    result <- optim(par = par, fn = modloglik_W, control = control, 
                    method = method, lower = c(0.0001, 0.0001, -Inf), 
                    upper = c(Inf, Inf, 0), duration = duration, 
                    modeloutput = modeloutput, funcenv = funcenv,  
                    cmatrix = cmatrix, nullmodel = nullmodel)  
    
  } else if (weightfunc == "G"){
    if (par[2] <= 0){
      stop("GEV scale parameter should be >0")
    }
    j      <- seq(-10, 10, by = (2 * 10 / duration))
    result <- optim(par = par, fn = modloglik_G, control = control, 
                    method = method, lower = c(-Inf, 0.0001, -Inf), 
                    upper = c(Inf, Inf, Inf), duration = duration, 
                    modeloutput = modeloutput, funcenv = funcenv,
                    cmatrix = cmatrix, nullmodel = nullmodel)
  } else {
    stop("Please choose Method to equal either W or G")
  } 
  bestmodel                     <- which(as.numeric(funcenv$DAICc) == min(as.numeric(funcenv$DAICc)))[1] # sometimes there are several bestmodels with similar DAICc, in which case we just pick one as they are all very similar
  WeightedOutput                <- list()   # prepare output of best model
  WeightedOutput$DeltaAICc      <- funcenv$DAICc[bestmodel]
  WeightedOutput$par_shape      <- funcenv$par_shape[bestmodel]
  WeightedOutput$par_scale      <- funcenv$par_scale[bestmodel]
  WeightedOutput$par_loc        <- funcenv$par_location[bestmodel]
  WeightedOutput$Function       <- func
  WeightedOutput$weightfunc     <- weightfunc
  
  ifelse (weightfunc == "W", 
          weight <- weibull3(x = j[1:duration], shape = as.numeric(funcenv$par_shape[bestmodel]), 
                                                scale = as.numeric(funcenv$par_scale[bestmodel]), 
                                                location = as.numeric(funcenv$par_location[bestmodel])),
          weight <- dgev(j[1:duration], loc = as.numeric(funcenv$par_location[bestmodel]), 
                         scale = as.numeric(funcenv$par_scale[bestmodel]), 
                         shape = as.numeric(funcenv$par_shape[bestmodel]), 
                         log = FALSE))
  
  weight[is.na(weight)] <- 0
  if (sum(weight) == 0){
    weight <- weight + 1
  }
  
  weight                <- weight / sum(weight) 
  modeldat$climate      <- apply(cmatrix, 1, FUN = function(x) {sum(x * weight)})
  LocalModel            <- update(modeloutput, .~., data = modeldat)
  WeightedOutput$Weight <- weight

  return(list(BestModel = LocalModel, BestModelData = model.frame(LocalModel), WeightedOutput = WeightedOutput))  
}
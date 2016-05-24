#' Error function
#'
#' @param x x
#' @return Error function result
#' @references https://stat.ethz.ch/R-manual/R-devel/library/stats/html/Normal.html
erf <- function(x) {
  2 * pnorm(x * sqrt(2)) - 1
}

#' Helper function hantushS
#'
#' @param x  distance between 0 and half lenght of recharge basin
#' @param alpha alpha
#' @param beta beta
#' @return Hantush star
#' @references p.22, http://pubs.usgs.gov/sir/2010/5102/support/sir2010-5102.pdf
hantushS <- function(x, alpha, beta) {
  xRoot <- sqrt(x)
  erf(alpha/xRoot) * erf(beta/xRoot)
}


#' Hantush function Sstar
#' @param alpha alpha
#' @param beta beta
#' @param dbg If True additional messages on integration quality of function
#' hantushSstar are printed on screen
#' @return Hantush Sstar result
#' @references p.22, http://pubs.usgs.gov/sir/2010/5102/support/sir2010-5102.pdf

hantushSstar <- function(alpha, beta, dbg) {
  
  x <- integrate(f = hantushS, lower = 0, upper = 1,
                 alpha = alpha, beta = beta # arguments passed to hansuhS
  )
  
  if (dbg) {
    msg <- sprintf("Sstar: %3.5f with an absolute error < %f", x$value, x$abs.error)
    message(msg)
  }
  return(x$value)
}


#' Hantush equation base properties
#'
#' @param time time elapsed since recharge began (T), (Default: 1.5 days)
#' @param basinWidth half width of the recharge basin (L), (Default: 10 m)
#' @param basinLength half length of the recharge basin (L), (Default: 10 m)
#' @param infiltrationRate recharge (infiltration) rate (L/T), (Default: 0.5 m/d)
#' @param horizConductivity horizontal hydraulic conductivity (L/T), (Default: 10 m/d)
#' @param iniHead initial head (height of the water table above the base of the
#' aquifer);(L), (Default: 10)
#' @param specificYield specific yield (Default: 0.2)
#' @param numberTimeSteps number of time steps to be used for average aquifer
#' thickness calculation (Default: 150)
#' @return Base properties for Hantush equation
#' @references p.22, http://pubs.usgs.gov/sir/2010/5102/support/sir2010-5102.pdf
baseProperties <- function(time = 10,
                            basinWidth = 10, ### m
                            basinLength = 10, ### m
                            infiltrationRate = 0.5 , ### 0.5 m/day
                            horizConductivity = 10, ### m/d
                            iniHead = 10, ### meter
                            specificYield = 0.2,
                            numberTimeSteps = 150) {
  
  x <- list(time = time,
            basinWidth = basinWidth,
            basinLength = basinLength,
            infiltrationRate = infiltrationRate,
            horizConductivity = horizConductivity,
            iniHead = iniHead,
            specificYield = specificYield,
            numberTimeSteps = numberTimeSteps)
  return(x)
  
}

#' Hantush equation
#'
#' @param x distance from the center of the recharge basin in the x direction (L)
#' @param y distance from the center of the recharge basin in the y direction (L)
#' @param baseProps basic model properties as retrieved by baseProperties()
#' @param dbg If True additional messages on integration quality of function
#' hantushSstar are printed on screen
#' @return Head at a given time after recharge begins
#' @references p.22, http://pubs.usgs.gov/sir/2010/5102/support/sir2010-5102.pdf
#' @seealso \code{\link{baseProperties}} for basic model properties

hantush <- function(x = 0,
                     y = 0,
                     baseProps = baseProperties(),
                     dbg = TRUE) {
  indices <- seq_len(baseProps$numberTimeSteps)
  
  timeSteps <- baseProps$time/baseProps$numberTimeSteps * indices
  
  h <- vector(length = baseProps$numberTimeSteps)
  
  L <- baseProps$basinLength
  sy <- baseProps$specificYield
  W <- baseProps$basinWidth
  
  for (i in indices)
  {
    if (i == 1) {
      avgHead <- baseProps$iniHead
    } else {
      avgHead <- 0.5 * (baseProps$iniHead + h[i - 1])
    }
    
    if (dbg) {
      cat(sprintf("Elapsed time: %4.5f ;  Havg: %4.5f Hi-1: %4.5f\n\n",
                  timeSteps[i], avgHead, h[i - 1] ))
    }
    
    diffusivity <- baseProps$horizConductivity * avgHead / sy
    
    factor1 <- (baseProps$infiltrationRate * avgHead * timeSteps[i]) / (2*sy)
    
    divisor <- sqrt( 4 * timeSteps[i] * baseProps$horizConductivity *  avgHead / baseProps$specificYield)
    
    alphaPlus <- (L + x) / divisor
    alphaMinus <- (L - x) / divisor
    
    betaPlus <- (W + y) / divisor
    betaMinus <- (W - y) / divisor
    
    factor2 <-
      hantushSstar(alphaPlus, betaPlus, dbg) +
      hantushSstar(alphaPlus, betaMinus, dbg) +
      hantushSstar(alphaMinus, betaPlus, dbg) +
      hantushSstar(alphaMinus, betaMinus, dbg)
    
    h[i] <- sqrt(factor1 * factor2 + baseProps$iniHead ^ 2)
    
  }
  
  res <- data.frame(timeSteps = timeSteps,
                    x = rep(x, baseProps$numberTimeSteps),
                    y = rep(y, baseProps$numberTimeSteps),
                    head = h)
  return(list(res = res,
              baseProps = baseProps))
}

#' Hantush distance: for for multiple coordinate
#'
#' @param x vector with distances from the center of the recharge basin in the x
#' direction (L), (Default: each meter between 0-100m)
#' @param y vector with distances from the center of the recharge basin in the y
#' direction (L), (Default: 0 times length of x)
#' @param config function as retrieved by hantush()
#' @param dbg If True additional debug messages are printed on screen
#' @return Head at a given time after recharge begins
#' @seealso \code{\link{hantush}} for parameterising the Hantush equation
hantushDistances <- function(x = 0:10,
                              y = rep(0, length(x)),
                              config = hantush,
                              dbg = TRUE) {
  out <- NULL
  baseProps <- list()
  for (i in 1:length(x))
  {
    if (i == 1) {
      tmp <- config(x = x[i], y = y[i], dbg = dbg)
      baseProps <- tmp$baseProps
      tmp <- tmp$res
    } else {
      tmp <- config(x = x[i], y = y[i], dbg = dbg)$res
    }
    out <- rbind(out, tmp)
  }
  
  if (dbg)
  {
    cat(sprintf("Head of %4.5f m at position %5.2f m (x), %5.2f (y)\n", out$head[i],x[i], y[i]))
  }
  
  
  maxSimTime <- max(out$timeSteps)
  
  simTime <- out[out$timeSteps == maxSimTime,c("x", "y","head") ]
  

  out <- list(timeSteps = out,
              simTime = simTime,
              baseProps = baseProps
  )
  
  return(out)
}

#' Hantush distances & base properties: allows input of vector of x,y coordinates and
#' also a vector for one of the base properties
#'
#' @param x vector with distances from the center of the recharge basin in the x
#' direction (L), (Default: every 5 meter between 0-200m)
#' @param y vector with distances from the center of the recharge basin in the y
#' direction (L), (Default: 0 times length of x)
#' @param baseProps as retrieved by baseProperties(), but one property is allowed
#' to be a vector (e.g. infiltrationRate = c(1,2,4))
#' @param dbg If True additional debug messages are printed on screen
#' @return List with sublists "dat" (x,y,head & WLincrease), "changedBaseProp.Name"
#' (name of base property with multiple values) and "baseProps" (complete base
#' parameterisation)
#' @seealso \code{\link{baseProperties}} for basic model properties
#' @examples
#'  baseProps <- baseProperties( time = 2^(0:6),
#'                              infiltrationRate = 1,
#'                              basinWidth = 10,
#'                              basinLength = 50,
#'                              horizConductivity = 10,
#'                              iniHead = 10,
#'                              specificYield = 0.2,
#'                              numberTimeSteps = 15)
#' res <- hantushDistancesBaseProps(baseProps = baseProps)
#' cols <- length(unique(res$dat[[res$changedBaseProp.Name]]))
#' mainTxt <- sprintf("Changed baseProperty: %s", res$changedBaseProp.Name)
#' xyplot(WLincrease ~ x,
#'        groups=res$dat[[res$changedBaseProp.Name]],
#'        data=res$dat,
#'        type="b",
#'        auto.key=list(columns=cols),
#'        main=mainTxt)

hantushDistancesBaseProps <- function(x = seq(0,200, 5),
                                      y = rep(0, length(x)),
                                      baseProps = baseProperties( time = 2 ^ (0:6),
                                                                  infiltrationRate = 1,
                                                                  basinWidth = 10,
                                                                  basinLength = 50,
                                                                  horizConductivity = 10,
                                                                  iniHead = 10,
                                                                  specificYield = 0.2),
                                      dbg = FALSE) {
  newRes <- NULL
  
  multiplePars <- unlist(lapply(X = baseProps, FUN = length))
  multipleParVals <- NULL
  changedBaseProp.Name <- ""
  
  if (any(multiplePars > 1)) {
    selMultiplePars <- which(multiplePars > 1)
    if (length(selMultiplePars) == 1)
    {
      multipleParVals <- baseProps[selMultiplePars]
      changedBaseProp.Name <- names(multipleParVals)
    } else {
      variedPars <- paste(names(selMultiplePars), collapse = ", ")
      msg <- sprintf("Only one base property can be a vector!
                     Now you defined multiple values the following base properties:\n%s", variedPars)
      stop(msg)
    }
  } else {
    multipleParVals <- baseProps$time ### select a random variable from baseProp
  }
  
  baseProps.sel <- baseProps
  
  for (parVal in as.vector(multipleParVals[[1]]))
  {
    
    if (any(multiplePars > 1))
    {
      cat(sprintf("Calculating results for changed baseProperty '%s': %3.5f  ... ",
                  changedBaseProp.Name,
                  parVal))
    }
    config <- function(x,y, dbg) {
      
      baseProps.sel[[ changedBaseProp.Name]] <- parVal
      
      hantush(x = x,
              y = y,
              baseProps = baseProps.sel,
              dbg = dbg)
    }
    
    res <- hantushDistances(x = x,
                            y = y,
                            config =  config,
                            dbg = dbg)
    
    
    res$simTime$WLincrease <- res$simTime$head -  res$baseProps$iniHead
    
    if (any(multiplePars > 1)) {
      changedBaseProp <- data.frame(rep(parVal, nrow(res$simTime)))
      names(changedBaseProp) <- names(multipleParVals)
      tmpRes <- cbind(changedBaseProp, res$simTime)
    } else {
      tmpRes <- res$simTime
    }
    
    newRes <- rbind(newRes, tmpRes)
    
    if (any(multiplePars > 1)) {
      cat("Done!\n") }
  }
  
  return(list(dat = newRes,
              changedBaseProp.Name = changedBaseProp.Name,
              baseProps = baseProps)
  )
}

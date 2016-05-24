#' Compute tide tables
#' 
#' @description Takes a data frame as input with date time high water and height information and returns a tide table
#' @references Horn, W. (1960) Some Recent Approaches to Tidal Problems. Int. Hydrogr. Rev. 37(2), 65-84
#' @references  Godin, Gabriel (1972) The Analysis of Tides. Toronto, 264pp
#' @references \url{http://tidesandcurrents.noaa.gov/publications/glossary2.pdf}
#' @references \url{http://www.bsh.de/de/Produkte/Buecher/Berichte_/Bericht50/BSH-Bericht50.pdf}
#' @param dataInput A data frame with the columns observation_date, observation_time, high_or_low_water and height. See attached data for correct formats.
#' @param otz The time zone of the observations
#' @param hwi The average of all intervals between the Moon's transit (upper or lower) over the Greenwich meridian and the following           high or low waters for all phases of the Moon is known as mean high water lunitidal interval and is abbreviated to high              water interval (hwi). Please only supply a value, when you are sure. Otherwise leave the default value "99:99" untouched.            hwi is then computed for you.
#' @param asdate A string indication the date you want the analysis to start with. Format: "yyyy/mm/dd".
#' @param astime A string indicating the time you want the analysis to start with. Format: "hh:mm:ss"
#' @param aedate A string indication the date you want the analysis to start with. Format: "yyyy/mm/dd".
#' @param aetime A string indicating the time you want the analysis to start with. Format: "hh:mm:ss"
#' @param ssdate Synthesis start date. This indicates the date you want your tide table to start with. Format: See above
#' @param sstime Synthesis start time. The starting time for your tide table. Format: See above
#' @param sedate Synthesis end date. Format: See above
#' @param setime Synthesis end time. Format: See above
#' @param stz Dummy for later extension to modify target time zone.
#' @return Returns a list with elements of the analysis, fitting and the tide table for given data 
#' \item{c.table}{The complete synthesis data as a data.table object}
#' \item{tide.table}{The tide table as a data.table object}
#' \item{lm.coeff}{Coefficients for the eight fitted linear models used in the synthesis}
#' \item{diff.analyse}{Time in days spanning the analysis}
#' \item{i.analyse}{How many different cases where used in the analysis}
#' @examples 
#' TideTable(dataInput=observation, asdate=observation$observation_date[1], 
#' astime=observation$observation_time[1], 
#' aedate="1991/05/01", aetime="21:00:00", ssdate="1995/01/01", 
#' sstime="00:00:00", sedate="1995/01/31", setime="21:00:00")
#' 
#' @import data.table
#' @import chron
#' @importFrom stats coef
#' @importFrom stats coefficients
#' @importFrom stats lm.fit
#' @importFrom stats sd
#' @export
TideTable <- function(dataInput, otz = 1, hwi = "99:99", asdate, astime, aedate, aetime, ssdate, sstime, sedate, setime, stz = 1) {
  
  chron.origin <- chron(dates. = "1900/01/01", 
                        times. = "00:00:00",
                        format = c(dates = "y/m/d", times = "h:m:s"),
                        out.format = c(dates = "y/m/d", times = "h:m:s"))
  tperiode.m2  <- 360 / 28.9841042
  tmean.moon   <- tperiode.m2 * 2
  tm24         <- tmean.moon / 24
  tmoon.0      <- chron(dates. = "1949/12/31", 
                        times. = "21:08:00",
                        format = c(dates = "y/m/d", times = "h:m:s"),
                        out.format = c(dates = "y/m/d", times = "h:m:s")) - chron.origin
                        
  tplus        <- tmoon.0 + 24.2491 / 1440.00
  
  if (unlist(strsplit(hwi, ":"))[1] == "99") {
    tmhwi <- 0
  } else {
    tmhwi <- as.numeric(unlist(strsplit(hwi, ":"))[1]) / 24 + as.numeric(unlist(strsplit(hwi, ":"))[2]) / 1440
  }

  chron.beob      <- chron(dates. = as.character(dataInput$observation_date),
                           times. = as.character(dataInput$observation_time),
                           format = c(dates = "y/m/d", times = "h:m:s"),
                           out.format = c(dates = "y/m/d", times = "h:m:s"))
                        
  diff.days       <- (chron.beob - chron.origin) - otz / 24 
  high.low        <- dataInput$high_or_low_water  
  nummk4          <- NumCulm(t = diff.days, tmhwi = tmhwi)
  tmmt.numm       <- nummk4$numm * tm24 + tplus
  stunden.transit <- (diff.days - tmmt.numm) * 24

  design.frame  <- data.table(numm            = nummk4$numm,
                              k               = nummk4$k4,
                              stunden.transit = stunden.transit,
                              height          = dataInput$height,
                              high.low        = high.low)
  
  #removing rows where k+high.low is odd
  design.frame <- design.frame[((k + high.low) %% 2) != 1]
  #Dropping the high.low column
  design.frame[, high.low := NULL]
  
  if (mean(design.frame[k == 1, (stunden.transit)]) / 24 < 0) {
    tmhwi           <- mean(design.frame[k == 1, (stunden.transit)]) / 24 + tm24 / 2 
    nummk4          <- NumCulm(t = diff.days, tmhwi = tmhwi)
    tmmt.numm       <- nummk4$numm * tm24 + tplus
    stunden.transit <- (diff.days - tmmt.numm) * 24
    numm            <- NULL
    design.frame    <- data.table(numm            = nummk4$numm,
                                  k               = nummk4$k4,
                                  stunden.transit = stunden.transit,
                                  height          = dataInput$height)
  }    

  asdate.time <- chron(dates. = asdate,
                       times. = astime,
                       format = c(dates = "y/m/d", times = "h:m:s"),
                       out.format = c(dates = "y/m/d", times = "h:m:s")) - chron.origin
  
  aedate.time <- chron(dates. = aedate,
                       times. = aetime,
                       format = c(dates = "y/m/d", times = "h:m:s"),
                       out.format = c(dates = "y/m/d", times = "h:m:s")) - chron.origin
  
  #Computing NumCulm
  astart.nummculm  <- NumCulm(t = asdate.time, tmhwi = tmhwi)
  
  aend.nummculm    <- NumCulm(t = aedate.time, tmhwi = tmhwi)
  
  tdiff.analyse    <- -astart.nummculm$numm+aend.nummculm$numm + 1
  
  #Subsetting design.frame
  
  design.frame     <- design.frame[(numm >= astart.nummculm$numm) & (numm <= aend.nummculm$numm)]
  
  lm.fitting          <- list()
  lm.fits             <- list()
  fitting.coef        <- list()
  i.analyse           <- matrix(nrow = 4, ncol = 2)
  colnames(i.analyse) <- c("stunden.transit", "height")
  rownames(i.analyse) <- 1 : 4
  temp.design         <- list()
  design.list         <- list()
  
  for(k4 in 1 : 4) {    
    predictor     <- design.frame[k == k4, c("numm", "stunden.transit", "height"), with = FALSE]$numm
    design.matrix <- matrix(nrow = length(predictor), ncol = 89)
    
    for(i in 1 : nrow(design.matrix)) {
      design.matrix[i, ] <- Funcs(xi = predictor[i], tdiff = tdiff.analyse)
    }
    
    for(l in c("stunden.transit", "height")) {
      col <- as.symbol(l)  
      predictant       <- design.frame[k == k4, l, with = FALSE]
      temp.design[[l]] <- data.table(design.matrix, predictant, predictor)    
      temp.design[[l]] <- temp.design[[l]][(eval(col) >= (mean(eval(col)) - 3 * sd(eval(col)))) & (eval(col) <= (mean(eval(col)) +                             3 * sd(eval(col))))]
      i.analyse[k4, l]  <- nrow(temp.design[[l]])
      lm.fitting[[l]]   <- lm.fit(x = as.matrix(temp.design[[l]][, !c(l,"predictor"), with = FALSE]),
                                 y = temp.design[[l]][, eval(col)]) 
    }
    lm.fits[[k4]]      <- lm.fitting
    design.list[[k4]]  <- temp.design
    fitting.coef[[k4]] <- lapply(lm.fitting, coef)
  }
  
  #Synthesis period
  ssdate.time <- chron(dates. = ssdate,
                       times. = sstime,
                       format = c(dates = "y/m/d", times = "h:m:s"),
                       out.format = c(dates = "y/m/d", times = "h:m:s")) - chron.origin
  sedate.time <- chron(dates. = sedate,
                       times. = setime,
                       format = c(dates = "y/m/d", times = "h:m:s"),
                       out.format = c(dates = "y/m/d", times = "h:m:s")) - chron.origin
  
  start.nummculm  <- NumCulm(t = ssdate.time, tmhwi = tmhwi)
  end.nummculm    <- NumCulm(t = sedate.time, tmhwi = tmhwi)
  #Synthesis
  time1        <- vector()
  height       <- vector()
  afunc        <- vector()
  coeff        <- vector()
  st.transit   <- vector()
  time.height  <- data.table(matrix(0.0, ncol = 6,
                                        nrow = ((end.nummculm$numm - start.nummculm$numm + 1) * 4)))
  m <- 0
  for (i in start.nummculm$numm : end.nummculm$numm) {
    afunc <- Funcs(i, tdiff = tdiff.analyse)
    for (k in 1 : 4) {
      m <- m + 1
      for (l in c("stunden.transit", "height")) {        
        coeff <- coefficients(lm.fits[[k]][[l]])
        summe <- coeff[1]
        for (h in 1 : 44){
          if (is.na(coeff[h * 2])) next
          summe <- summe + coeff[h * 2] * (afunc[h * 2]) + coeff[h * 2 + 1] * (afunc[h * 2 + 1]) 
        }
        if (l == "stunden.transit") {
          st.transit[m] <- summe
          tmmt.numm     <- i * tm24 + tplus
          time1[m]      <- (tmmt.numm + summe / 24 + otz / 24)
          
        }
        else {
          height[m] <- summe
        }   
      }
      if (k == 1 | k == 3){
        ihn <- 1
      } else {
        ihn <- 0
      }
      if (k == 1 | k == 2) {
        trans <- 1
      } else {
        trans <- 0
      }
      time.height[m,] <- as.list(c(chron(dates. = time1[m],
                                        origin. = c(month = 1, day = 1, year = 1900)),
                                  ihn,
                                  trans,
                                  height[m],
                                  st.transit[m],
                                  i))  
    }
  }
  date_time       <- NULL
  prediction_date <- NULL
  prediction_time <- NULL
  V1              <- NULL
  time.height[, date_time := chron(dates.  = V1,
                                  origin. = c(month = 1, day = 1, year = 1900))]
  time.height[, prediction_date := dates(date_time)]
  time.height[, prediction_time := paste(hours(date_time), 
                                         minutes(date_time), 
                                         seconds(date_time), 
                                         sep = ":")]
  setnames(time.height, "V6", "i")
  setnames(time.height, "V5", "st.transit")
  setnames(time.height, "V4", "height")
  setnames(time.height, "V3", "upper_or_lower_transit")
  setnames(time.height, "V2", "high_or_low_water")
  setnames(time.height, "V1", "transit")
  time.height[, date_time := NULL]
  setcolorder(time.height, c("transit",
                             "prediction_date",
                             "prediction_time",
                             "high_or_low_water",
                             "upper_or_lower_transit",
                             "height",
                             "st.transit",
                             "i"))
  tide.table <- time.height[,c("prediction_date",
                             "prediction_time",
                              "high_or_low_water",
                               "height"), with=FALSE]
  #we return a list called report containing the tide table, diff.analyse, i.analyse and lm.coeff
  report              <- list()
  report$c.table      <- time.height
  report$tide.table   <- tide.table
  report$diff.analyse <- tdiff.analyse
  report$i.analyse    <- i.analyse
  report$lm.coeff     <- fitting.coef
  
  return(report)
}


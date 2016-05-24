#' Computes tide curves
#'
#' @description Takes a data frame as input with three columns (see example dataset) and returns a tide curve.
#' Internally the analysis is carried out in lunar days.
#' One mean lunar day lasts 1.0350501 mean solar days.
#' Therefore the analysis time period
#' should start one lunar day after the first observation and end one lunar day before the last observation.
#' @param dataInput A data frame with the columns observation_date, observation_time and height. See attached data for correct    formats.
#' @param otz The time zone of the observations
#' @param km The number of nodes between two consecutive mean moon transits. Shall be less or equal to: round(1440 [min] / time step [min])
#' Example: Time step 5 min: Use km = 288 or even smaller. Leave on default (km = -1) and supply mindt, when unsure.
#' @param mindt Observation time step in [min]
#' @param asdate A string indication the date you want the analysis to start with. Format: "yyyy/mm/dd".
#' @param astime A string indicating the time you want the analysis to start with. Format: "hh:mm:ss"
#' @param aedate A string indication the date you want the analysis to start with. Format: "yyyy/mm/dd".
#' @param aetime A string indicating the time you want the analysis to start with. Format: "hh:mm:ss"
#' @param ssdate Synthesis start date. This indicates the date you want your tide curve to start with. Format: See above
#' @param sstime Synthesis start time. The starting time for your tide table. Format: See above
#' @param sedate Synthesis end date. Format: See above
#' @param setime Synthesis end time. Format: See above
#' @return Returns a list with elements of the analysis, fitting and the tidal curve for given data
#' \item{synthesis.lunar}{The lunar synthesis data as a data.table object}
#' \item{data.matrix}{The data needed for analysis}
#' \item{tide.curve}{The solar tide curve as a data.table object}
#' \item{lm.coeff}{Coefficients for the km fitted linear models used in the synthesis}
#' \item{diff.analyse}{Time in days spanning the analysis}
#' \item{i.analyse}{How many different cases where used in the analysis}
#' @references  Godin, Gabriel (1972) The Analysis of Tides. Toronto, 264pp
#' @references \url{http://tidesandcurrents.noaa.gov/publications/glossary2.pdf}
#' @references \url{http://www.bsh.de/de/Produkte/Buecher/Berichte_/Bericht50/BSH-Bericht50.pdf}
#' @examples
#' TideCurve(dataInput = tideObservation, asdate = "2015/12/06",
#'              astime = "00:00:00",      aedate = "2015/12/31",
#'              aetime = "23:30:00",      ssdate = "2015/12/17",
#'              sstime = "00:00:00",      sedate = "2015/12/31",
#'              setime = "23:30:00")
#'
#' @import data.table
#' @import chron
#' @importFrom stats coef
#' @importFrom stats coefficients
#' @importFrom stats lm.fit
#' @importFrom stats sd
#' @importFrom fields splint
#' @export
TideCurve <- function(dataInput, otz = 1, km = -1, mindt = 30, asdate, astime, aedate, aetime, ssdate, sstime, sedate, setime) {

  if(km == -1){
    km <- round(1440 / mindt)
  }

  height  <- dataInput$height
  nspline <- 7

  chron.origin <- chron(dates. = "1900/01/01",
                        times. = "00:00:00",
                        format = c(dates = "y/m/d", times = "h:m:s"),
                        out.format = c(dates = "y/m/d", times = "h:m:s"))
  tperiode.m2  <- 360 / 28.9841042373
  tmean.moon   <- tperiode.m2 * 2
  tm24         <- tmean.moon / 24
  tmoon.0      <- as.numeric(chron(dates. = "1949/12/31",
                                   times. = "21:08:00",
                                   format = c(dates = "y/m/d", times = "h:m:s"),
                                   out.format = c(dates = "y/m/d", times = "h:m:s")) - chron.origin)
  tmdm         <- 24.2325 / 1440
  tplus        <- tmoon.0 + tmdm
  tmmh         <- tm24 / km
  tdtobs       <- mindt / 1440
  tdtobsn      <- tdtobs * (nspline - 1)
  tdtobsn.105  <- tdtobsn + 10^-5
  tdtobs.105   <- tdtobs + 10^-5
  tmondkm      <- numeric()

  for (i in 1:km) {
    tmondkm[i] <- tmmh * (i - 0.5)
  }

  chron.beob      <- chron(dates. = as.character(dataInput$observation_date),
                           times. = as.character(dataInput$observation_time),
                           format = c(dates = "y/m/d", times = "h:m:s"),
                           out.format = c(dates = "y/m/d", times = "h:m:s"))

  diff.days       <- as.numeric((chron.beob - chron.origin) - otz / 24)
  length.diffdays <- length(diff.days)
  moona           <- as.numeric(floor((diff.days[1] - tplus) / tm24))
  moone           <- as.numeric(floor((diff.days[length.diffdays] - tplus) / tm24))
  tmmt.numm       <- numeric()

  for(i in moona:moone){
    tmmt.numm[i] <- i * tm24 + tplus
  }

  asdate.time <- chron(dates. = asdate,
                       times. = astime,
                       format = c(dates = "y/m/d", times = "h:m:s"),
                       out.format = c(dates = "y/m/d", times = "h:m:s")) - chron.origin - (otz / 24)

  aedate.time <- chron(dates. = aedate,
                       times. = aetime,
                       format = c(dates = "y/m/d", times = "h:m:s"),
                       out.format = c(dates = "y/m/d", times = "h:m:s")) - chron.origin - (otz / 24)

  numma   <- as.numeric(floor((asdate.time - tplus) / tm24))
  numme   <- as.numeric(floor((aedate.time - tplus) / tm24))
  xa      <- numeric()
  ya      <- numeric()
  ty      <- numeric()

  data.matrix <- matrix(0.0, nrow = length.diffdays, ncol = 4)
  colnames(data.matrix) <- c("numm","imm", "tmmttmond", "height")
  numm <- NULL
  data.matrix <- data.table(data.matrix)

  #can put the first floor call outside the loop?
  floored   <- floor((diff.days - tplus) / tm24)
  tdtobs.2  <- tdtobs / 2
  imm       <- numeric()
  tmmttmond <- numeric()
  tmd       <- numeric()
  dx        <- numeric()

  for (ii in 4 : (length.diffdays - 3)) {
    #ik        <- floor((diff.days[i] - tplus) / tm24)
    ik <- floored[ii]
    if((ik < numma) || (ik > numme)) next
    imm       <- floor((diff.days[ii] - tmmt.numm[ik]) / tmmh) + 1
    tmmttmond <- tmmt.numm[ik] + tmondkm [imm]
    tmd       <- diff.days[ii] - tmmttmond
    if (abs(tmd) > tdtobs.2) next
    insp <- 0
    for (j in (ii - 3) : (ii + 3)){
      insp     <- insp + 1
      xa[insp] <- diff.days[j]
      ya[insp] <- height[j]

    }
    dx       <- xa[insp] - xa[1]
    if(dx > tdtobsn.105) next
    ty                 <- splint(xa, ya, tmmttmond)

    set(data.matrix,i = ii, j = 1L, value = ik)
    set(data.matrix,i = ii, j = 2L, value = imm)
    set(data.matrix,i = ii, j = 3L, value = tmmttmond)
    set(data.matrix,i = ii, j = 4L, value = ty)

  }

  tdiff.analyse    <- numme - numma + 1
  design.frame     <- data.matrix[(numm >= numma) & (numm <= numme)]

  lm.fits             <- list()
  fitting.coef        <- list()
  i.analyse           <- matrix(nrow = km, ncol = 1)
  rownames(i.analyse) <- 1 : km
  temp.design         <- list()
  design.list         <- list()

  for(k in 1 : km) {
    predictor     <- design.frame[imm == k, "numm", with = FALSE]$numm
    design.matrix <- matrix(nrow = length(predictor), ncol = 89)

    for(i in 1 : nrow(design.matrix)) {
      design.matrix[i, ] <- Funcs(xi = predictor[i], tdiff = tdiff.analyse)
    }

    predictant    <- design.frame[imm == k, "height",  with = FALSE]
    temp.design   <- data.table(design.matrix, predictant, predictor)
    temp.design   <- temp.design[(height >= (mean(height) - 3 * sd(height))) & (height <= (mean(height) +  3 * sd(height)))]
    i.analyse[k,] <- nrow(temp.design)
    lm.fitting    <- lm.fit(x = as.matrix(temp.design[, !c("predictor","height"), with = FALSE]),
                            y = as.matrix(temp.design[,"height", with = FALSE]))

    lm.fits[[k]]      <- lm.fitting
    design.list[[k]]  <- temp.design
    fitting.coef[[k]] <- coef(lm.fitting)
  }

  #Synthesis
  ssdate.time <- chron(dates. = ssdate,
                       times. = sstime,
                       format = c(dates = "y/m/d", times = "h:m:s"),
                       out.format = c(dates = "y/m/d", times = "h:m:s")) - chron.origin - otz / 24
  sedate.time <- chron(dates. = sedate,
                       times. = setime,
                       format = c(dates = "y/m/d", times = "h:m:s"),
                       out.format = c(dates = "y/m/d", times = "h:m:s")) - chron.origin - otz / 24

  nummsa  <- as.numeric(floor((ssdate.time - tplus) / tm24))
  nummse  <- as.numeric(floor((sedate.time - tplus) / tm24))
  m.length     <- (nummse - nummsa + 1) * km
  time1        <- numeric(length = m.length)
  height       <- numeric(length = m.length)
  tobstime     <- numeric(length = m.length)
  afunc        <- vector()
  coeff        <- vector()
  time.height           <- matrix(1.0, nrow = (m.length), ncol = 4)
  colnames(time.height) <- c("height", "i", "k", "time1")
  time.height           <- data.table(time.height)
  options(chron.origin = c(month = 1, day = 1, year = 1900))

  m <- 0L
  for (ii in nummsa : nummse) {
    afunc <- Funcs(ii, tdiff = tdiff.analyse)
    for (k in 1 : km) {
      m <- m + 1L
      coeff <- coefficients(lm.fits[[k]])
      summe <- coeff[1]
      for (h in 1 : 44){
        if (is.na(coeff[h * 2])) next
        summe <- summe + coeff[h * 2] * (afunc[h * 2]) + coeff[h * 2 + 1] * (afunc[h * 2 + 1])
      }
      time1[m]      <- ii * tm24 + tplus + (k - 0.5) * tmmh
      height[m]     <- summe
      tobstime[m]   <- round(time1[m] / tdtobs) * tdtobs

      set(time.height,i = m, j = 1L, value = height[m])
      set(time.height,i = m, j = 2L, value = ii)
      set(time.height,i = m, j = 3L, value = k)
      set(time.height,i = m, j = 4L, value = time1[m])
    }
  }
  time.height[,date_time := chron(dates. = time1)]
  setcolorder(time.height, c("date_time", "height", "i", "k", "time1"))

  l <- 1L

  tsyntstd <- numeric()
  ty       <- numeric()

  for(j in 4 : (m - 3)){
    tsyn  <- tobstime[j]
    tsynp <- tsyn + tdtobs.2
    tsynm <- tsyn - tdtobs.2
    if((tsynp < ssdate.time) || (tsynm > sedate.time)) next

    insp <- 0
    for(is in (j - 3) : (j + 3)){
      insp     <- insp + 1
      xa[insp] <- time1[is]
      ya[insp] <- height[is]
    }

    tdtt <- tobstime[j] - tobstime [j - 1]

    if(tdtt > (tdtobs.105)){

      tsynn    <- tobstime[j] - tdtobs
      ty[l]       <- splint(xa, ya, tsynn)
      tsyntstd[l] <- tsynn + otz / 24
      l <- l + 1L
    }
    ty[l]      <- splint(xa, ya, tsyn)
    tsyntstd[l] <- tsyn + otz / 24
    l <- l + 1L
  }
  #Add date/time columns to synthesis
  prediction_date <- NULL
  prediction_time <- NULL
  date_time   <- NULL
  tidal.curve <- data.table(date_time = chron(dates. = (tsyntstd + 1 / 864000)), height = ty)


  tidal.curve[, prediction_date := strftime(dates(date_time), format = "%Y/%m/%d")]
  tidal.curve[, prediction_time := paste(hours(date_time),
                                           minutes(date_time),
                                           seconds(date_time),
                                           sep = ":")]
  time.height[, prediction_date := strftime(dates(date_time), format = "%Y/%m/%d")]
  time.height[, prediction_time := paste(hours(date_time),
                                         minutes(date_time),
                                         seconds(date_time),
                                         sep = ":")]

  #we return a list called report containing the tide curve (lunar and solar), diff.analyse, i.analyse and lm.coeff
  report                 <- list()
  report$data.matrix     <- data.matrix
  report$synthesis.lunar <- time.height
  report$tide.curve      <- tidal.curve
  report$lm.coeff        <- fitting.coef
  report$i.analyse       <- diff.days
  report$diff.analyse    <- tdiff.analyse
  return(report)
}



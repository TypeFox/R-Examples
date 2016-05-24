

InitLogEntry <- function(dateIn, ticker, capital, timeIn = NA, execVol = 0,
                         execQuant = 0, avgPrice = NA, depthIn = NA) {
  ## formalize argument ##
  if (capital < 0) {
    stop('capital must be non-negative')
  }  
  dateIn <- as.numeric(dateIn)
  timeIn <- as.numeric(timeIn)
  if (class(ticker) == 'numeric') {
    ticker <- sprintf('%06d', ticker)
  }
  
  ## calculate returning entries ##
  capLeft   = capital - execQuant
  execRatio = execQuant / capital
  
  out <- data.frame(DATEIN      = dateIn, 
                    TIMEIN      = timeIn,
                    TICKER      = ticker, 
                    AVGPRICEIN  = avgPrice,
                    VOLUMEIN    = execVol,
                    QUANTIN     = execQuant,
                    CAPALLOW    = capital,
                    CAPLEFT     = capLeft,
                    EXECRATIO   = execRatio,
                    DEPTHIN     = depthIn,
                    DATEOUT     = NA,
                    TIMEOUT     = NA,
                    AVGPRICEOUT = NA,
                    VOLUMEOUT   = 0,
                    QUANTOUT    = 0,
                    DEPTHOUT    = NA,
                    SELLATTEMPT = 0,
                    VOLUMEHOLD  = execVol)
  
  return (out)
}
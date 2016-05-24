apcal <-
  function(lat,days,rad_mea,extraT=NULL, DL=NULL, SSD) {
    i <- dayOfYear(days)
    
    if (is.null(extraT) | is.null(DL))
    {
      ex <- extrat(lat=radians(lat),i)
      DL <- ex$DayLength   #[hours]
      extraT <- ex$ExtraTerrestrialSolarRadiationDaily  # [MJ]
    }
    Y <- rad_mea/extraT      
    X <- SSD/DL
    m <- lm(Y ~ X)
    rval <- c(m$coefficients,summary(m)$r.squared)
    names(rval) <- c("APa","APb","APr2")
    rval
  }


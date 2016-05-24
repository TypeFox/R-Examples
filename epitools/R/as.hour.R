"as.hour" <-
function(x, mindt, maxdt, half.hour = FALSE){
  if(half.hour==TRUE){
    unit <- 3600/2
  } else {
    unit <- 3600
  }
  ct <- as.POSIXct(x)
  lt <- as.POSIXlt(ct)
  mindt <- as.POSIXct(mindt)
  maxdt <- as.POSIXct(maxdt)
  ct.sec <- as.numeric(ct)
  ct.hour <- ct.sec%/%unit
  ct.hour12 <- as.numeric(format(ct, format = "%I"))
  ct.wkday <- format(ct, format = "%a")
  ct.weekday <- format(ct, format = "%A")
  cct <- seq(mindt, maxdt, 1)
  cct.sec <- as.numeric(cct)
  cct.hour <- cct.sec%/%unit
  cct.hour.tab <- as.numeric(names(table(cct.sec%/%unit)))
  cct.tab <- cct[!duplicated(cct.hour)]
  cct.hour12 <- as.numeric(format(cct.tab, format = "%I"))
  cct.ampm <- cct.ampm2 <- format(cct.tab, format = "%p")
  cct.ampm2[cct.ampm=="AM" & !is.na(cct.ampm)] <- "am"
  cct.ampm2[cct.ampm=="PM" & !is.na(cct.ampm)] <- "pm"
  cct.weekday <- format(cct.tab, format = "%A")
  cct.wkday <- format(cct.tab, format = "%a")
  cct.month <- format(cct.tab, format = "%B")
  cct.mon <- format(cct.tab, format = "%b")  
  cct.year <- format(cct.tab, format = "%Y")  
  names(cct.tab) <- cct.hour.tab
  ct.hour.factor <- factor(ct.hour, levels = cct.hour.tab)
  ct.stratum <- cct.tab[as.character(ct.hour)]
  ct.stratum.factor <- factor(unname(ct.stratum),
                       levels = unname(cct.tab))
  clt <- as.POSIXlt(cct.tab)
  list(ct = ct,
       sec = lt$sec,
       min = lt$min,
       hour = lt$hour,
       hour12 = ct.hour12,
       stratum = unname(ct.hour.factor),
       stratum2 = unname(ct.stratum),
       stratum3 = ct.stratum.factor,
       cstratum = unname(cct.hour.tab),
       cstratum2 = unname(cct.tab),
       csec = clt$sec,
       cmin = clt$min,
       chour = clt$hour,
       chour12 = cct.hour12,
       campm = cct.ampm,
       campm2 = cct.ampm2,
       cweekday = cct.weekday,
       cwkday = cct.wkday,
       cmday = clt$mday,
       cmonth = cct.month,
       cmon = cct.mon,
       cyear = cct.year,
       half.hour = half.hour
       )
}

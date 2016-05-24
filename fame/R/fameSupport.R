fameDateString <- function(xTi){
  tif <- tif(xTi)
  freq <- frequency(xTi)
  ans <- character(length(xTi))
  intraday <- isIntradayTif(tif)

  if(any(intraday)){
    intradayTif <- tif[intraday]
    firstTif <- intradayTif[1]
    if(length(intradayTif) > 1 && any(intradayTif != firstTif))
      stop("intraday args must all have same frequency")
    fmt <-paste("%d%b%Y:",
                c("%H","%H:%M","%H:%M:%S")[(firstTif %/% 1000) - 1],
                sep = "")
    ans[intraday] <- format(xTi[intraday], fmt)
  }
  if(any(!intraday)){
    z <- xTi[!intraday]
    ans[!intraday] <- paste(year(z), cycle(z), sep = ":")
  }
  ans
}

convertFreqCode <- function(tif = NULL, fame = NULL){
  ## Vectors of Fame and ti (time index) frequency codes
  tiFreqs   <- 1000 + (1:50)[-11]  ## Fame has no "reserves" frequency
  fameFreqs <- c(8:9, 16:22, 32, 64:77, 128:129, 144:145, 160:162, 192:209)
  if(!missing(tif)){ ## converting R tif to Fame freq
    tif[tif == 1011] <- 1015
    z <- numeric(length(tif))
    intraday <- isIntradayTif(tif)
    if(any(intraday)){
      ## fame uses 228 for hourly, 227 for minutely, 226 for secondly
      ## R uses 2000 for hourly, 3000 for minutely, 4000 for secondly
      base <- 230 - (tif[intraday] %/% 1000)
      nUnits <- tif[intraday] %% 1000
      z[intraday] <- base + 65536*nUnits
    }
    if(any(!intraday)){
      z[!intraday] <- fameFreqs[match(tif[!intraday], tiFreqs, nomatch = NA)]
    }
  }
  if(!missing(fame)){ ## converting Fame freq to R tif
    z <- numeric(length(fame))
    base <- fame %% 65536
    intraday <- between(base, 226, 228)
    if(any(intraday)){
      nUnits <- fame[intraday] %/% 65536
      z[intraday] <- (230 - base[intraday])*1000 + nUnits
    }
    if(any(!intraday))
      z[!intraday] <- tiFreqs[match(fame[!intraday], fameFreqs, nomatch = NA)]
  }
  z
}

tifToFame <- function(tif){
  nTif <- tif(tif)
  convertFreqCode(tif = nTif)
}

fameToTif <- function(fameFreq)
  convertFreqCode(fame = fameFreq)

tifToFameName <- function(tif){
  nTif <- tif(tif)
  cTif <- tifName(tif)
  intraday <- isIntradayTif(nTif)
  if(any(intraday)){
    base <- c("hourly", "minutely", "secondly")[(nTif[intraday] %/% 1000) -1]
    nUnits <- nTif[intraday] %% 1000
    nUnits[nUnits == 0] <- 1
    cTif[intraday] <- paste(base, "(", nUnits, ")", sep = "")
  }
  if(any(!intraday))
    cTif[!intraday] <- c("daily",
                         "business",
                         "weekly(sunday)",
                         "weekly(monday)",
                         "weekly(tuesday)",
                         "weekly(wednesday)",
                         "weekly(thursday)",
                         "weekly(friday)",
                         "weekly(saturday)",
                         "tenday",
                         "biweekly(awednesday)",
                         "biweekly(asunday)",
                         "biweekly(amonday)",
                         "biweekly(atuesday)",
                         "biweekly(awednesday)",
                         "biweekly(athursday)",
                         "biweekly(afriday)",
                         "biweekly(asaturday)",
                         "biweekly(bsunday)",
                         "biweekly(bmonday)",
                         "biweekly(btuesday)",
                         "biweekly(bwednesday)",
                         "biweekly(bthursday)",
                         "biweekly(bfriday)",
                         "biweekly(bsaturday)",
                         "twicemonthly",
                         "monthly",
                         "bimonthly(november)",
                         "bimonthly(december)",
                         "quarterly(october)",
                         "quarterly(november)",
                         "quarterly(december)",
                         "annual(january)",
                         "annual(february)",
                         "annual(march)",
                         "annual(april)",
                         "annual(may)",
                         "annual(june)",
                         "annual(july)",
                         "annual(august)",
                         "annual(september)",
                         "annual(october)",
                         "annual(november)",
                         "annual(december)",
                         "semiannual(july)",
                         "semiannual(august)",
                         "semiannual(september)",
                         "semiannual(october)",
                         "semiannual(november)",
                         "semiannual(december)")[nTif[!intraday] - 1000]
  return(cTif)
}

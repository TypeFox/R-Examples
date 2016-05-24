hffmcRaster <- function(weatherstream, ffmc_old = 85, time.step = 1, 
                        hourlyFWI = FALSE) {
  #############################################################################
  # Description: Raster-based Diurnal (Hourly) Fine Fuel Moisture Code 
  #              Calculation. Most of  the equations in this code refer to the 
  #              Van Wagner (1977), withsome equations contained in Van Wagner 
  #              & Pickett (1985). Additionally, some modifications were made 
  #              for precision.
  #
  #              Van Wagner, C.E. 1977. A method of computing fine fuel moisture
  #              content throughout the diurnal cycle. Environment Canada, 
  #              Canadian Forestry Service, Petawawa Forest Experiment Station, 
  #              Chalk River, Ontario. Information Report PS-X-69. 
  #              http://cfs.nrcan.gc.ca/pubwarehouse/pdfs/25591.pdf
  # 
  #              Equations and FORTRAN program for the Canadian Forest Fire 
  #              Weather Index System. 1985. Van Wagner, C.E.; Pickett, T.L. 
  #              Canadian Forestry Service, Petawawa National Forestry 
  #              Institute, Chalk River, Ontario. Forestry Technical Report 33. 
  #              18 p.
  #
  #              Additional reference on FWI system
  #
  #              Development and structure of the Canadian Forest Fire Weather 
  #              Index System. 1987. Van Wagner, C.E. Canadian Forestry Service,
  #              Headquarters, Ottawa. Forestry Technical Report 35. 35 p.
  #  
  #
  # Args: weatherstream:   Input weather stream data.frame which includes
  #                        temperature, relative humidity, wind speed, 
  #                        precipitation, hourly value, and bui. More specific
  #                        info can be found in the hffmc.Rd help file.
  #            ffmc_old:   ffmc from previous timestep
  #           time.step:   The time (hours) between previous FFMC and current
  #                        time.
  #           hourlyFWI:   Can calculated hourly ISI & FWI as well 
  #                        (TRUE/FALSE, default=FALSE)
  #
  # Returns: A single or multiple hourly ffmc raster
  #############################################################################
  names(weatherstream) <- tolower(names(weatherstream))
  #local scope variables
  Tp <- weatherstream$temp
  H  <- weatherstream$rh
  W  <- weatherstream$ws
  ro <- weatherstream$prec
  #Check that the parameters are correct
  if (is.numeric(ffmc_old)){
    Fo <- Tp
    values(Fo)<-ffmc_old
  }else{Fo<-ffmc_old}
  if (!exists("Tp") | is.null(Tp)) 
    warning("temperature (temp) is missing!")
  if (!exists("ro") | is.null(ro)) 
    warning("precipitation (prec) is missing!")
  if (!exists("W") | is.null(W)) 
    warning("wind speed (ws) is missing!")
  if (!exists("H") | is.null(H)) 
    warning("relative humidity (rh) is missing!")
  #Eq. 1 (with a more precise multiplier than the daily)
  mo <- 147.27723 * (101 - Fo) / (59.5 + Fo)
  
  mr1 <- mo
  mr1[mr1 > 150] <- NA
  #masking values
  rf1 <- mask(ro, mr1)
  #Eqs. 3a (Van Wagner & Pickett 1985)
  mr1 <- mr1 + 42.5 * rf1 * exp(-100 / (251 - mr1)) * (1 - exp(-6.93 / rf1))
  mr2 <- mo
  mr2[mr2 <= 150] <- NA
  rf2 <- mask(ro, mr2)
  #Eqs. 3b (Van Wagner & Pickett 1985)
  mr2 <- mr2 + 42.5 * rf2 * exp(-100 / (251 - mr2)) *(1 - exp(-6.93 / rf2)) + 
         0.0015 * ((mr2 - 150)^2) * (rf2^0.5)
  mr3 <- cover(mr1,mr2)
  #The real moisture content of pine litter ranges up to about 250 percent,
  # so we cap it at 250
  mr3[mr3 > 250] <- 250
  #raster manipulation to speed up processing
  r1 <- ro
  r1[r1 <= 0] <- NA
  mr<-mask(mr3, r1)
  r1 <- ro
  r1[r1 > 0] <- NA
  mo1 <- mask(mo, r1)
  mo <- cover(mo1, mr)
  #Eq. 2a Equilibrium moisture content from drying
  Ed <- 0.942 * (H^0.679) + 11 * exp((H - 100)/10) + 0.18 * (21.1 - Tp) * 
        (1 - exp(-0.115 * H))
  #Eq. 3a Log drying rate at the normal temperature of 21.1C
  ko <- 0.424 * (1 - (H / 100)^1.7) + 0.0694 * (W^0.5) * (1 - (H / 100)^8)
  #Eq. 3b
  kd <- ko * 0.0579 * exp(0.0365 * Tp)
  #Eq. 8 (Van Wagner & Pickett 1985)
  md <- Ed + (mo - Ed) * (10^(-1*kd*time.step))
  #Eq. 2b Equilibrium moisture content from wetting  
  Ew <- 0.618 * (H^0.753) + 10 * exp((H - 100) / 10) + 0.18 * (21.1 - Tp) * 
        (1 - exp(-0.115 * H))
  #Eq. 7a Log wetting rate at the normal temperature of 21.1 C 
  k1 <- 0.424 * (1 - ((100 - H)/100)^1.7) + 0.0694 * (W^0.5) * 
        (1 - ((100 - H) / 100)^8)
  #Eq. 4b
  kw <- k1 * 0.0579 * exp(0.0365 * Tp)
  #Eq. 8 (Van Wagner & Pickett 1985)
  mw <- Ew - (Ew - mo) * (10^(-1*kw*time.step)) 
  #Constraints using raster manipulation
  m0 <- overlay(mo, Ed, fun = function(a, b){ return(a > b) })
  md[m0 == 0] <- NA
  mw[m0 == 1] <- NA
  m <- cover(md, mw)

  m1 <- overlay(Ed, mo, Ew, fun = function(a, b, c) return(a >= b & b >= c))
  mo[m1 == 0] <- NA
  m[m1 == 1] <- NA
  m <- cover(mo, m)
  #Eq. 6 - Final hffmc calculation (modified 3rd constant to 147.27723)
  fo <- 59.5 * (250 - m) / (147.27723 + m)
  fo[fo <= 0] <- 0
  #Calculate hourly isi and fwi
  if (hourlyFWI){
    if ("bui" %in% names(weatherstream)){
      bui <- weatherstream$bui

      #Calculate ISI
      fW <- exp(0.05039 * W)
      fm <- 147.2 * (101 - fo) / (59.5 + fo)
      fF <- 91.9 * exp(-0.1386 * fm) * (1 + (fm^5.31) / 49300000)
      isi <- 0.208 * fW * fF
      
      #Calculate FWI
      bui1 <- bui
      bui1[bui1 <= 80] <- NA
      bui1 <- 0.1 * isi * (1000 / (25 + 108.64 / exp(0.023 * bui1)))
      
      bui2 <- bui
      bui2[bui1 > 80] <- NA
      bui2 <- 0.1 * isi * (0.626 * (bui2^0.809) + 2)
      bb <- cover(bui1, bui2)
      
      bb1 <- bb
      bb1[bb > 1] <- NA
      
      bb2 <- bb
      bb2[bb <= 1] <- NA
      bb2 <- exp(2.72 * ((0.434 * log(bb2))^0.647))
      fwi <- cover(bb1, bb2)
      #Calculate DSR
      dsr <- 0.0272 * (fwi^1.77)
      #Create Raster Stack for the ouput
      output <- stack(fo, isi, fwi, dsr)
      names(output) <- c("hffmc", "hisi", "hfwi", "hdsr")
      return(output)
    } else {
      warning("Daily BUI is required to calculate hourly FWI")
    }
  } else {
    names(fo) <- "hffmc"
    return(fo)
  }
}


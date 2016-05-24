fwiRaster <- function(input, init = c(ffmc = 85, dmc = 6, dc = 15), mon = 7,
                      out = "all", lat.adjust = TRUE, uppercase = TRUE) {
  #############################################################################
  # Description: Raster-based implementation of the Canadian Forest Fire Weather 
  #              Index Calculations. All code is based on a C code library that 
  #              was written by Canadian Forest Service Employees, which was 
  #              originally based on the Fortran code listed in the reference 
  #              below. All equations in this code refer to that document, 
  #              unless otherwise noted.
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
  # Args: init:     Initializing moisture values
  #                 ffmc:     Fine Fuel Moisture Code (default 85)
  #                 dmc:      Duff Moisture Code (default 6)
  #                 dc:       Drought Code (default 15)
  #                 lat:      Latitude (decimal degrees, default 55)
  #       batch:    Function can be run in a batch mode, where multiple 
  #                 weather stations or points can be calculated at once. 
  #                 (TRUE/FALSE, default TRUE)
  #       out:      Display the calculated FWI values, with or without the 
  #                 inputs. (all/fwi, default all)
  #       lat.adjust: Option to adjust day length in the calculations 
  #                   (TRUE/FALSE, default TRUE)
  #       uppercase:  Output names in upper or lower case - a commonly 
  #                   asked for feature, as dataset naming conventions vary 
  #                   considerably. (TRUE/FALSE, default TRUE)
  #       
  #
  # Returns: A data.frame of the calculated FWI values with or without
  #          the input data attached to it.
  #
  #############################################################################
  #Reference latitude for DMC day length adjustment
  #46N: Canadian standard, latitude >= 30N   (Van Wagner 1987)
  ell01 <- c(6.5, 7.5, 9, 12.8, 13.9, 13.9, 12.4, 10.9, 9.4, 8, 7, 6)
  #20N: For 30 > latitude >= 10
  ell02 <- c(7.9, 8.4, 8.9, 9.5, 9.9, 10.2, 10.1, 9.7, 9.1,8.6, 8.1, 7.8)
  #20S: For -10 > latitude >= -30  
  ell03 <- c(10.1, 9.6, 9.1, 8.5, 8.1, 7.8, 7.9, 8.3, 8.9, 9.4, 9.9, 10.2)
  #40S: For -30 > latitude
  ell04 <- c(11.5, 10.5, 9.2, 7.9, 6.8, 6.2, 6.5, 7.4, 8.7, 10, 11.2, 11.8)
  #For latitude near the equator, we simple use a factor of 9 for all months
  
  #Day length factor for DC Calculations
  #20N: North of 20 degrees N
  fl01 <- c(-1.6, -1.6, -1.6, 0.9, 3.8, 5.8, 6.4, 5, 2.4, 0.4, -1.6, -1.6)
  #20S: South of 20 degrees S
  fl02 <- c(6.4, 5, 2.4, 0.4, -1.6, -1.6, -1.6, -1.6, -1.6, 0.9, 3.8, 5.8)
  #Near the equator, we just use 1.4 for all months.
  
  #Quite often users will have a data frame called "input" already attached
  #  to the workspace. To mitigate this, we remove that if it exists, and warn
  #  the user of this case.
  if (!is.na(charmatch("input", search()))) {
    detach(input)
  }
  names(input) <- tolower(names(input))
  temp <- input$temp
  prec <- input$prec
  ws <- input$ws
  rh <- input$rh
  if ("lat" %in% names(input)) {
    lat <- input$lat
  }else {
    lat <- temp
    values(lat) <- 55
  }
  
  if (!exists("temp") | is.null(temp)) 
    warning("temperature (temp) is missing!")
  if (!exists("prec") | is.null(prec)) 
    warning("precipitation (prec) is missing!")
  if (!exists("ws") | is.null(ws)) 
    warning("wind speed (ws) is missing!")
  if (!exists("rh") | is.null(rh)) 
    warning("relative humidity (rh) is missing!")

  names(init) <- tolower(names(init))

  #Assign values for initializing variables
  if (is.numeric(init)){
    if (is.null(names(init))){
      names(init)<-c('ffmc', 'dmc', 'dc')
    }
    ffmc_yda <- dmc_yda <- dc_yda <- temp
    values(ffmc_yda) <- init[['ffmc']]
    values(dmc_yda) <- init[['dmc']]
    values(dc_yda) <- init[['dc']]
  } else {
    ffmc_yda <- init$ffmc
    dmc_yda  <- init$dmc
    dc_yda   <- init$dc
  }
  #constrain relative humidity
  rh[rh>=100]<- 99.9999
  ###########################################################################
  #                    Fine Fuel Moisture Code (FFMC)
  ###########################################################################
  #Eq. 1
  wmo <- 147.2 * (101 - ffmc_yda)/(59.5 + ffmc_yda)
  #Eq. 2 Rain reduction to allow for loss in overhead canopy
  ra1 <- prec
  ra1[ra1 <= 0.5] <- NA
  ra1 <- ra1-0.5
  ra2 <- prec
  ra2[ra2 > 0.5] <- NA
  ra <- cover(ra1, ra2)
  #masking values
  wmo1 <- mask(wmo, ra1)
  wmo2 <- mask(wmo,ra2)
  wmo11 <- wmo1
  wmo11[wmo11 <= 150] <- NA
  ra11 <- ra1
  ra11[wmo1 <= 150] <- NA
  #Eqs. 3a & 3b
  wmo11 <- wmo11 + 0.0015 * (wmo11 - 150) * (wmo11 - 150) * sqrt(ra11) + 42.5 * 
           ra11 * exp(-100 / (251 - wmo11)) * (1 - exp(-6.93 / ra11))
  wmo12 <- wmo1
  wmo12[wmo12 > 150] <- NA
  ra12 <- ra1
  ra12[wmo1 > 150] <- NA
  wmo12 <- wmo12 + 42.5 * ra12 * exp(-100 / (251 - wmo12)) * 
          (1 - exp(-6.93 / ra12))
  wmo1 <- cover(wmo11, wmo12)
  wmo <- cover(wmo1, wmo2)
  #The real moisture content of pine litter ranges up to about 250 percent,
  # so we cap it at 250
  wmo[wmo > 250] <- 250
  #cleanup intermediate values
  rm(ra1, ra11, ra12, ra2, wmo1, wmo2, wmo11, wmo12)
  #Eq. 4 Equilibrium moisture content from drying
  ed <- 0.942 * (rh^0.679) + (11 * exp((rh - 100)/10)) + 0.18 * 
    (21.1 - temp) * (1 - 1/exp(rh * 0.115))
  #Eq. 5 Equilibrium moisture content from wetting
  ew <- 0.618 * (rh^0.753) + (10 * exp((rh - 100)/10)) + 0.18 * 
    (21.1 - temp) * (1 - 1/exp(rh * 0.115))
  #Create a new raster object based on wmo, ed, and ew
  z0 <- overlay(wmo, ed, ew, fun = function(a, b, c){ return(a < b & a < c) })
  #Create new rasters and mask out missing values
  z0[z0 == 0] <- NA
  rh0 <- mask(rh, z0)
  ws0 <- mask(ws, z0)
  #Eq. 6a (ko) Log drying rate at the normal termperature of 21.1 C
  z <- 0.424 * (1 - (((100 - rh0)/100)^1.7)) + 0.0694 * sqrt(ws0) * (1 - ((100 - rh0)/100)^8)
  # Assigning to 0 instead of NA, as 0 makes more sense
  z[is.na(z)] <- 0
  # Mask missing temp values
  z <- mask(z, temp)
  rm(rh0, ws0, z0)
  #Eq. 6b Affect of temperature on  drying rate
  x <- z * 0.581 * exp(0.0365 * temp)
  #Create a new raster object based on wmo, ed, and ew
  z0 <- overlay(wmo, ed, ew, fun = function(a, b, c){ return(a < b & a < c) })
  #Create new rasters and mask out missing values
  z0[z0 == 0] <- NA
  ew0 <- mask(ew, z0)
  x0 <- mask(x, z0)
  wmo0 <- mask(wmo, z0)
  #Eq. 8
  wmo1 <- ew0 - (ew0 - wmo0) / (10^x0)
  wmo2 <- wmo
  wmo2[!is.na(wmo0)] <- NA
  wm <- cover(wmo1, wmo2)
  rm(z0, ew0, x0, wmo0, wmo1, wmo2)
  #Create a new raster object based on wmo, and ed
  z0 <- overlay(wmo, ed, fun = function(a, b){ return(a > b) }) 
  #Create new rasters and mask out missing values
  z0[z0 == 0] <- NA
  rh0 <- mask(rh, z0)
  ws0 <- mask(ws, z0)
  #Eq. 7a (ko) Log wetting rate at the normal termperature of 21.1 C   
  z0 <- 0.424 * (1 - (rh0 / 100)^1.7) + 0.0694 * sqrt(ws0) * (1 - (rh0 / 100)^8)
  z1 <- z
  z1[!is.na(z0)] <- NA
  z <- cover(z0, z1)
  rm(rh0, ws0)
  #Eq. 7b Affect of temperature on  wetting rate
  x <- z * 0.581 * exp(0.0365 * temp)
  ed0 <- mask(ed,z0)
  wmo0 <- mask(wmo,z0)
  x0 <- mask(x,z0)
  #Eq. 9
  wm0 <- ed0 + (wmo0 - ed0)/(10^x0)
  wm1 <- mask(wm, z1)
  wm <- cover(wm0, wm1)
  rm(ed0, x0, wm0, wm1, wmo0)
  #Eq. 10 Final FFMC calculation
  ffmc <- (59.5 * (250 - wm))/(147.2 + wm)
  #Constraints
  ffmc[ffmc>101] <- 101 
  ffmc[ffmc<0] <- 0

  
  ###########################################################################
  #                        Duff Moisture Code (DMC)
  ###########################################################################
  t0 <- temp
  #constrain low end of temperature
  t0[t0 < -1.1] <- -1.1
  #Eq. 16 - The log drying rate
  rk <- 1.894 * (t0 + 1.1) * (100 - rh) * ell01[mon] * 1e-04
  #Adjust the day length  and thus the drying r, based on latitude and month
  if (lat.adjust) {
    rk[lat <= 30 & lat > 10] <- 1.894 * (t0[lat <= 30 & lat > 10] + 1.1) * 
      (100 - rh[lat <= 30 & lat > 10]) * ell02[mon] * 1e-04
    rk[lat <= -10 & lat > -30] <- 1.894 * (t0[lat <= -10 & lat > -30] + 1.1) * 
      (100 - rh[lat <= -10 & lat > -30]) * ell03[mon] * 1e-04
    rk[lat <= -30 & lat >= -90] <- 1.894 * (t0[lat <= -30 & lat >= -90] + 1.1) *
      (100 - rh[lat <= -30 & lat >= -90]) * ell04[mon] * 1e-04
    rk[lat <= 10 & lat > -10] <- 1.894 * (t0[lat <= 10 & lat > -10] + 1.1) * 
      (100 - rh[lat <= 10 & lat > -10]) * 9 * 1e-04
  }
  ra <- prec
  #Eq. 11 - Net rain amount
  rw <- 0.92 * ra - 1.27
  #Alteration to Eq. 12 to calculate more accurately
  wmi <- 20 + 280 / exp(0.023 * dmc_yda)
  #Eqs. 13a, 13b, 13c
  b <- dmc_yda
  b[dmc_yda <= 33] <- 100 / (0.5 + 0.3 * dmc_yda[dmc_yda <= 33])
  if (!is.null(dmc_yda[dmc_yda > 33 & dmc_yda <= 65])){
    b[dmc_yda > 33 & dmc_yda <= 65] <- 
      14 - 1.3 * log(dmc_yda[dmc_yda > 33 & dmc_yda <= 65])
  }
  if(!is.null(dmc_yda[dmc_yda > 65])){
    b[dmc_yda > 65] <- 
      6.2 * log(dmc_yda[dmc_yda > 65]) - 17.2
  }
  #Eq. 14 - Moisture content after rain
  wmr <- wmi + 1000 * rw / (48.77 + b * rw)
  op <- options(warn = (-1))
  #Alteration to Eq. 15 to calculate more accurately
  pr0 <- 43.43 * (5.6348 - log(wmr - 20))
  options(op)

  pr<-pr0
  #Constrain P
  pr[prec <= 1.5] <-dmc_yda[prec <= 1.5]
  pr[pr < 0] <- 0
  #Calculate final P (DMC)
  dmc <- pr + rk
  dmc[dmc < 0] <- 0 
  ###########################################################################
  #                             Drought Code (DC)
  ###########################################################################
  #Constrain temperature
  t0[temp< (-2.8)] <- -2.8
  #Eq. 22 - Potential Evapotranspiration
  pe <- (0.36 * (t0 + 2.8) + fl01[mon])/2
  #Daylength factor adjustment by latitude for Potential Evapotranspiration
  if (lat.adjust) {
    pe[lat <= -10] <- (0.36 * (t0[lat <= -10] + 2.8) + fl02[mon]) / 2
    pe[lat > -10 & lat <= 10] <- (0.36 * (t0[lat > -10 & lat <= 10] + 2.8) + 1.4) / 2
  }
  ra <- prec
  #Eq. 18 - Effective Rainfall
  rw <- 0.83 * ra - 1.27
  #Eq. 19
  smi <- 800 * exp(-1 * dc_yda / 400)
  #Alteration to Eq. 21
  dr0 <- dc_yda - 400 * log(1 + 3.937 * rw / smi)
  dr0[dr0 < 0] <- 0 
  dr <- dr0
  #if precip is less than 2.8 then use yesterday's DC
  dr[prec <= 2.8] <- dc_yda[prec <= 2.8]
  #Alteration to Eq. 23
  dc <- dr + pe
  dc[dc < 0] <- 0
  
  ###########################################################################
  #                    Initial Spread Index (ISI)
  ###########################################################################
  #Eq. 24 - Wind Effect
  fW <- exp(0.05039 * ws)
  #Eq. 10 - Moisture content
  fm <- 147.2 * (101 - ffmc) / (59.5 + ffmc)
  #Eq. 25 - Fine Fuel Moisture
  fF <- 91.9 * exp(-0.1386 * fm) * (1 + (fm^5.31) / 49300000)
  #Eq. 26 - Spread Index Equation
  isi <- 0.208 * fW * fF
  
  ###########################################################################
  #                       Buildup Index (BUI)
  ###########################################################################
  #Eq. 27a
  bui <- 0.8 * dc * dmc/(dmc + 0.4 * dc)
  bui[dmc == 0 & dc == 0] <- 0
  #Eq. 27b - next 4 lines
  p <- (dmc - bui)/dmc
  p[dmc == 0] <- 0
  cc <- 0.92 + ((0.0114 * dmc)^1.7)
  bui0 <- dmc - cc * p
  #Constraints
  bui0[bui0 < 0] <- 0
  bui[bui < dmc] <- bui0[bui < dmc]
  
  ###########################################################################
  #                     Fire Weather Index (FWI)
  ###########################################################################
  #Eqs. 28b, 28a, 29
  bb <-0.1 * isi * (0.626 * (bui^0.809) + 2)
  bb[bui > 80] <- 0.1 * isi[bui > 80] * (1000 / (25 + 108.64 / exp(0.023 * bui[bui > 80])))
  #Eqs. 30b, 30a
  fwi <-exp(2.72 * ((0.434 * log(bb))^0.647))
  #Constraint
  fwi[bb <= 1] <- bb[bb <= 1]
  ###########################################################################
  #                   Daily Severity Rating (DSR)
  ###########################################################################
  #Eq. 31
  dsr <- 0.0272 * (fwi^1.77)

  #If output specified is "fwi", then return only the FWI variables
  if (out == "fwi") {
    #Creating a raster stack of FWI variables to return
    new_FWI <- stack(ffmc, dmc, dc, isi, bui, fwi, dsr)
    names(new_FWI) <- c("ffmc", "dmc", "dc", "isi", "bui", "fwi", "dsr")
    if (uppercase){
      names(new_FWI) <- toupper(names(new_FWI))
    }
    #If output specified is "all", then return both FWI and input weather vars
  } else {
    if (out == "all") {
      #Create a raster stack of input and FWI variables
      new_FWI <- stack(input, ffmc, dmc, dc, isi, bui, fwi, dsr)
      names(new_FWI) <- c(names(input),"ffmc", "dmc", "dc", "isi", "bui", "fwi", "dsr")
      if (uppercase){
        names(new_FWI) <- toupper(names(new_FWI))
      }
    }
  }
  return(new_FWI)
}


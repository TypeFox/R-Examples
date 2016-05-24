.dmcCalc <- function(dmc_yda, temp, rh, prec, lat, mon, lat.adjust=TRUE) {
  #############################################################################
  # Description: Duff Moisture Code Calculation. All code
  #              is based on a C code library that was written by Canadian
  #              Forest Service Employees, which was originally based on
  #              the Fortran code listed in the reference below. All equations
  #              in this code refer to that document.
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
  # Args:  dmc_yda:   The Duff Moisture Code from previous iteration
  #           temp:   Temperature (centigrade)
  #             rh:   Relative Humidity (%)
  #           prec:   Precipitation(mm)
  #            lat:   Latitude (decimal degrees)
  #            mon:   Month (1-12)
  #     lat.adjust:   Latitude adjustment (TRUE, FALSE, default=TRUE)
  #       
  #
  # Returns: A single dmc value
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
  
  #constrain low end of temperature
  temp <- ifelse(temp < (-1.1), -1.1, temp)
  #Eq. 16 - The log drying rate
  rk <- 1.894 * (temp + 1.1) * (100 - rh) * ell01[mon] * 1e-04
  #Adjust the day length  and thus the drying r, based on latitude and month
  if (lat.adjust) {
    rk <- ifelse(lat <= 30 & lat > 10, 1.894 * (temp + 1.1) * 
                   (100 - rh) * ell02[mon] * 1e-04, rk)
    rk <- ifelse(lat <= -10 & lat > -30, 1.894 * (temp + 1.1) * 
                   (100 - rh) * ell03[mon] * 1e-04, rk)
    rk <- ifelse(lat <= -30 & lat >= -90, 1.894 * (temp + 1.1) * 
                   (100 - rh) * ell04[mon] * 1e-04, rk)
    rk <- ifelse(lat <= 10 & lat > -10, 1.894 * (temp + 1.1) * 
                   (100 - rh) * 9 * 1e-04, rk)
  }
  ra <- prec
  #Eq. 11 - Net rain amount
  rw <- 0.92 * ra - 1.27
  #Alteration to Eq. 12 to calculate more accurately
  wmi <- 20 + 280/exp(0.023 * dmc_yda)
  #Eqs. 13a, 13b, 13c
  b <- ifelse(dmc_yda <= 33, 
              100/(0.5 + 0.3 * dmc_yda), 
              ifelse(dmc_yda <= 65, 
                     14 - 1.3 * log(dmc_yda), 
                     6.2 * log(dmc_yda) - 17.2))
  #Eq. 14 - Moisture content after rain
  wmr <- wmi + 1000 * rw/(48.77 + b * rw)
  op <- options(warn = (-1))
  #Alteration to Eq. 15 to calculate more accurately
  pr0 <- 43.43 * (5.6348 - log(wmr - 20))
  options(op)
  #Constrain P
  pr <- ifelse(prec <= 1.5, dmc_yda, pr0)
  pr <- ifelse(pr < 0, 0, pr)
  #Calculate final P (DMC)
  dmc1 <- pr + rk
  dmc1 <- ifelse(dmc1 < 0, 0, dmc1)
  return(dmc1)
}
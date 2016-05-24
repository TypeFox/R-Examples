gfmc <- function(input, GFMCold = 85, batch = TRUE, time.step = 1, roFL = 0.3,
                 out = "GFMCandMC") {
  #############################################################################
  # Description: Calculation of the Grass Fuel Moisture Code. This calculates
  #              the moisture content of both the surface of a fully cured
  #              matted grass layer and also an equivalent Grass Fuel Moisture
  #              Code. All equations come from Wotton (2009) as cited below
  #              unless otherwise specified.
  # 
  #              Wotton, B.M. 2009. A grass moisture model for the Canadian 
  #              Forest Fire Danger Rating System. In: Proceedings 8th Fire and
  #              Forest Meteorology Symposium, Kalispell, MT Oct 13-15, 2009. 
  #              Paper 3-2. https://ams.confex.com/ams/pdfpapers/155930.pdf
  #           
  # Args: input (data.frame):
  #         temp (required)	    Temperature (centigrade)
  #         rh	 (required)	    Relative humidity (%)
  #         ws	 (required)	    10-m height wind speed (km/h)
  #         prec (required)	    1-hour rainfall (mm)
  #         isol (required)	    Solar radiation (kW/m^2)
  #         mon	 (recommended)	Month of the year (integer 1-12)
  #         day	 (optional)	    Day of the month (integer)
  #       GFMCold:    GFMC from yesterday (double, default=85)
  #       batch:      Compute multiple locations (TRUE/FALSE, default=TRUE)
  #       time.step:  The hourly time steps (integer hour, default=1)
  #       roFL:       Nominal fuel load of the fine fuel layer 
  #                   (kg/m^2 double, default=0.3)
  #       out:        Output format (GFMCandMC/MC/GFMC/ALL, default=GFMCandMC)
  #       
  # Returns: Returns a data.frame of either MC, GMFC, All, or GFMCandMC
  #
  #############################################################################
  t0 <- time.step
  names(input) <- tolower(names(input))
  #Quite often users will have a data frame called "input" already attached
  #  to the workspace. To mitigate this, we remove that if it exists, and warn
  #  the user of this case.
  if (!is.na(charmatch("input", search()))) {
    warning("Attached dataset 'input' is being detached to use fbp() function.")
    detach(input)
  }
  #set local scope variables
  temp <- input$temp
  prec <- input$prec
  ws <- input$ws
  rh <- input$rh
  isol <- input$isol

  #show warnings when inputs are missing
  if (!exists("temp") | is.null(temp)) 
    warning("temperature (temp) is missing!")
  if (!exists("prec") | is.null(prec)) 
    warning("precipitation (prec) is missing!")
  if (!exists("ws") | is.null(ws)) 
    warning("wind speed (ws) is missing!")
  if (!exists("rh") | is.null(rh)) 
    warning("relative humidity (rh) is missing!")
  if (!exists("isol") | is.null(isol)) 
    warning("ISOL is missing!")  

  #check for issues with batching the function
  if (batch){
    if ("id" %in% names(input)) {
      n <- length(unique(input$id))
      if(length(unique(input[1:n, "id"])) != n){
        stop("Multiple stations have to start and end at the same dates/time, 
             and input data must be sorted by date/time and id")
      }
    } else {
      n <- 1
    }
  } else {n <- nrow(input)}
  
  if (length(temp)%%n != 0)
    warning("Input data do not match with number of weather stations")

  if (length(GFMCold) != n & length(GFMCold) == 1){
    warning("One GFMCold value for multiple weather stations") 
    GFMCold <- rep(GFMCold, n)
  }
  
  if (length(GFMCold) != n & length(GFMCold) > 1)
    stop("Number of GFMCold doesn't match number of wx stations")
  validOutTypes = c("GFMCandMC", "MC", "GFMC", "ALL")
  if(!(out %in% validOutTypes)){
    stop(paste("'",out, "' is an invalid 'out' type.", sep=""))
  }
  
  #get the length of the data stream
  n0 <- length(temp)%/%n
  GFMC <- NULL
  MC <- NULL
  #iterate through timesteps
  for (i in 1:n0){
    #k is the data for all stations by time step
    k <- (n * (i - 1) + 1):(n * i)
    #Eq. 13 - Calculate previous moisture code
    MCold <- 147.2772 * ((101 - GFMCold) / (59.5 + GFMCold))
    #Eq. 11 - Calculate the moisture content of the layer in % after rainfall
    MCr <- ifelse(prec[k] > 0, MCold + 100 * (prec[k] / roFL), MCold)
    #Constrain to 250
    MCr <- ifelse(MCr > 250, 250, MCr)
    MCold <- MCr
    #Eq. 2 - Calculate Fuel temperature
    Tf <- temp[k] + 35.07 * isol[k] * exp(-0.06215 * ws[k])
    #Eq. 3 - Calculate Saturation Vapour Pressure (Baumgartner et a. 1982)
    eS.T <- 6.107 * 10^(7.5 * temp[k] / (237 + temp[k]))
    #Eq. 3 for Fuel temperature
    eS.Tf <- 6.107 * 10^(7.5 * Tf / (237 + Tf))
    #Eq. 4 - Calculate Fuel Level Relative Humidity
    RH.f <- rh[k] * (eS.T / eS.Tf)
    #Eq. 7 - Calculate Equilibrium Moisture Content for Drying phase
    EMC.D <- (1.62 * RH.f^0.532 + 13.7 * exp((RH.f - 100) / 13.0)) + 
              0.27 * (26.7 - Tf) * (1 - exp(-0.115 * RH.f))
    #Eq. 7 - Calculate Equilibrium Moisture Content for Wetting phase
    EMC.W <- (1.42 * RH.f^0.512 + 12.0 * exp((RH.f - 100) / 18.0)) + 
              0.27 * (26.7 - Tf) * (1 - exp(-0.115 * RH.f))
    #RH in terms of RH/100 for desorption
    Rf <- ifelse(MCold > EMC.D, RH.f / 100, rh)
    #RH in terms of 1-RH/100 for absorption
    Rf <- ifelse(MCold < EMC.W, (100 - RH.f) / 100, Rf)
    #Eq. 10 - Calculate Inverse Response time of grass (hours)
    K.GRASS <- 0.389633 * exp(0.0365 * Tf) * (0.424 * (1 - Rf^1.7) + 0.0694 * 
                                              sqrt(ws[k]) * (1 - Rf^8))
    #Fuel is drying, calculate Moisture Content
    MC0 <- ifelse(MCold > EMC.D, EMC.D + (MCold - EMC.D) * 
                  exp(-1.0 * log(10.0) * K.GRASS * t0), MCold)
    #Fuel is wetting, calculate moisture content
    MC0 <- ifelse(MCold < EMC.W, EMC.W + (MCold - EMC.W) * 
                  exp(-1.0 * log(10.0) * K.GRASS * t0), MC0)

    #Eq. 12 - Calculate GFMC
    GFMC0 <- 59.5 * ((250 - MC0) / (147.2772 + MC0))
    #Keep current and old GFMC
    GFMC <- c(GFMC, GFMC0)
    #Same for moisture content
    MC <- c(MC, MC0)
    #Reset vars
    GFMCold <- GFMC0
    MCold <- MC0
  }
  #Return requested 'out' type
  if (out=="ALL"){
    return(as.data.frame(cbind(input, GFMC, MC)))
  } else if(out == "GFMC"){
    return(GFMC)
  } else if (out == "MC"){
    return(MC)
  } else { #GFMCandMC
    return(data.frame(GFMC = GFMC, MC = MC))
  }
}


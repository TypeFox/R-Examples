fwi <- function(input, init = data.frame(ffmc = 85, dmc = 6, dc = 15, lat = 55),
              batch = TRUE, out = "all", lat.adjust = TRUE, uppercase = TRUE) {
  #############################################################################
  # Description: Canadian Forest Fire Weather Index Calculations. All code
  #              is based on a C code library that was written by Canadian
  #              Forest Service Employees, which was originally based on
  #              the Fortran code listed in the reference below. All equations
  #              in this code refer to that document, unless otherwise noted.
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
  #Args:  input:    View Documentation (fwi.Rd) for full description
  #                 of input data frame
  #       init:     Initializing moisture values
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
  
  #Quite often users will have a data frame called "input" already attached
  #  to the workspace. To mitigate this, we remove that if it exists, and warn
  #  the user of this case.
  if (!is.na(charmatch("input", search()))) {
    detach(input)
  }
  names(input) <- tolower(names(input))
  
  #convert vector to data.frame to ensure consitency
  if (is.vector(init)){
    init <- as.data.frame(t(init))
  }
  names(init) <- tolower(names(init))
  #resolve missing names of the initializing variables if necessary
  if(substr(names(init), 1, 1)[1] == "x" | substr(names(init), 1, 1)[1] == "v"){
    if (ncol(init) == 3){
      names(init) <- c("ffmc", "dmc", "dc")
      init$lat <- 55
    }else if(ncol(init) == 4){
      names(init) <- c("ffmc", "dmc", "dc", "lat")
    }
  }
    
  #############################################################################
  #                                 
  # Set local variables and display warnings to user if default is being used
  #############################################################################
  ffmc_yda <- init$ffmc
  dmc_yda  <- init$dmc
  dc_yda   <- init$dc

  if ("lat" %in% names(input)) {
    lat <- input$lat
  }
  else {
    warning("latitude was not provided, assign default value 55")
    lat <- rep(55, nrow(input))
  }
  if ("long" %in% names(input)) {
    long <- input$long
  }
  else {
    warning("long was not provided, assign a default number -120")
    long <- rep(-120, nrow(input))
  }
  if ("yr" %in% names(input)) {
    yr <- input$yr
  }
  else {
    warning("Year was not provided, assigned default number 5000")
    yr <- rep(5000, nrow(input))
  }
  if ("mon" %in% names(input)) {
    mon <- input$mon
  }
  else {
    warning("Month was not provided, assigned the default value, July")
    mon <- rep(7, nrow(input))
  }
  if ("day" %in% names(input)) {
    day <- input$day
  }
  else {
    warning("day was not provided, assigned default number -99")
    day <- rep(-99, nrow(input))
  }

  #If batch selected, then sort the data by Date and id and determine the 
  # length of each run.
  # Currently when running multiple stations, the stations much have the same
  # amount of data and same start/end dates
  #Function stops running if these requirements are not met
  if (batch){
    if ("id" %in% names(input)) {
      input <- input[with(input,order(yr,mon,day,id)),]
      #number of stations
      n <- length(unique(input$id))
      if(length(unique(input[1:n, "id"])) != n){
        stop("Multiple stations have to start and end at the same dates, and 
             input data must be sorted by date/time and id")
      }
    } else {
      n <- 1
    }
  }else{
    n <- nrow(input)
  }

  temp <- input$temp
  prec <- input$prec
  ws <- input$ws
  rh <- input$rh
  if (!exists("temp") | is.null(temp)) 
    warning("temperature (temp) is missing!")
  if (!exists("prec") | is.null(prec)) 
    warning("precipitation (prec) is missing!")
  if (!exists("ws") | is.null(ws)) 
    warning("wind speed (ws) is missing!")
  if (!exists("rh") | is.null(rh)) 
    warning("relative humidity (rh) is missing!")
  if (length(unique(!mon %in% 1:12))>1)
    warning("Month has to be between 1 and 12")
  if (length(unique(!day %in% 1:31))>1)
    warning("Day has to be between 1 and 31")
  if (length(unique(lat>90))>1|length(unique(lat< -90))>1)
    warning("Latitude has to be between -90 and 90")
  if (length(unique(long>180))>1|length(unique(long< -180))>1)
    warning("Longitude has to be between -180 and 180")
  #############################################################################
  #                                 END
  # Set local variables and display warnings to user if default is being used
  #############################################################################

  if (length(temp) %% n != 0)
    warning("Missing records may generate wrong outputs")
  if (nrow(init) == 1 & n > 1){
    warning("Same initial data were used for multiple weather stations")
    ffmc_yda <- rep(ffmc_yda, n)
    dmc_yda <- rep(dmc_yda, n)
    dc_yda <- rep(dc_yda, n)
  }
  #if the number of rows in the init file does not equal that of the number of
  # stations, then stop execution as we do not have a complete input set
  if(nrow(init) > 1 & nrow(init) != n) {
    stop("Number of initial values do not match with number of weather 
         stations")
  }
  
  #Length of weather run
  n0 <- length(temp) / n
  #Initialize variables
  ffmc <- dmc <- dc <- isi <- bui <- fwi <- dsr <- NULL
  #For each day in the run
  for (i in 1:n0){
    #k is the data for all stations by day
    k  <- ((i - 1) * n + 1):(i * n)
    #constrain relative humidity
    rh[k] <- ifelse(rh[k] >= 100, 99.9999, rh[k])
    ###########################################################################
    # Fine Fuel Moisture Code (FFMC)
    ###########################################################################
    ffmc1 = .ffmcCalc(ffmc_yda, temp[k], rh[k], ws[k], prec[k])
    
    ###########################################################################
    # Duff Moisture Code (DMC)
    ###########################################################################
    dmc1 = .dmcCalc(dmc_yda, temp[k], rh[k], prec[k], lat[k], mon[k], 
                    lat.adjust)
    
    ###########################################################################
    # Drought Code (DC)
    ###########################################################################
    dc1 <- .dcCalc(dc_yda, temp[k], rh[k], prec[k], lat[k], mon[k],
                   lat.adjust)
    
    ###########################################################################
    # Initial Spread Index (ISI)
    ###########################################################################
    isi1 <- .ISIcalc(ffmc1, ws[k], FALSE)
    
    ###########################################################################
    # Buildup Index (BUI)
    ###########################################################################
    bui1 <- .buiCalc(dmc1, dc1)
    
    ###########################################################################
    # Fire Weather Index (FWI)
    ###########################################################################
    fwi1 <- .fwiCalc(isi1, bui1)
    ###########################################################################
    #                   Daily Severity Rating (DSR)
    ###########################################################################
    #Eq. 31
    dsr1 <- 0.0272 * (fwi1^1.77)
    
    #Concatenate values
    ffmc<-c(ffmc,ffmc1)
    dmc<-c(dmc,dmc1)
    dc<-c(dc,dc1)
    isi<-c(isi,isi1)
    bui<-c(bui,bui1)
    fwi<-c(fwi,fwi1)
    dsr<-c(dsr,dsr1)
    ffmc_yda<-ffmc1
    dmc_yda<-dmc1
    dc_yda<-dc1
  } 
  
  #If output specified is "fwi", then return only the FWI variables
  if (out == "fwi") {
    new_FWI <- data.frame(ffmc = ffmc, dmc = dmc, dc = dc, isi = isi, 
                          bui = bui, fwi = fwi, dsr = dsr)
    if (uppercase){
      names(new_FWI) <- toupper(names(new_FWI))
    }
  }
  #If output specified is "all", then return both FWI and input weather vars
  else {
    if (out == "all") {
      new_FWI <- cbind(input, ffmc, dmc, dc, isi, bui, fwi, dsr)
      if (uppercase){
        names(new_FWI) <- toupper(names(new_FWI))
      }
    }
  }
  return(new_FWI)
}

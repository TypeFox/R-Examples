sdmc <- function(input, sdmc_old = NULL, batch = TRUE){
  #############################################################################
  # Description: Calculate the sheltered Duff Moisture Code (sDMC) based on
  #              noon weather observations of temp, rh, ws, 24-hour rain and 
  #              previous day's sDMC. Equations being referenced are from
  #              Wotton et. al. (2005) or Van Wagner & Pickett (1985).
  #
  #              Wotton, B.M., B.J. Stocks, and D.L. Martell. 2005. An index for
  #              tracking sheltered forest floor moisture within the Canadian 
  #              Forest Fire Weather Index System. International Journal of 
  #              Wildland Fire, 14, 169-182.
  #
  #              Equations and FORTRAN program for the Canadian Forest Fire 
  #              Weather Index System. 1985. Van Wagner, C.E.; Pickett, T.L. 
  #              Canadian Forestry Service, Petawawa National Forestry 
  #              Institute, Chalk River, Ontario. Forestry Technical Report 33. 
  #              18 p.
  #     
  #  
  # Args:  
  #       input:  View Documentation (sdmc.Rd) for full description of input
  #               data frame
  #    sdmc_old:  previous day's calculated sDMC value
  #       batch:  Function can be run in a batch mode, where multiple 
  #               weather stations or points can be calculated at once. 
  #               (TRUE/FALSE, default TRUE)
  #       
  #
  # Returns: sdmc single or vector of SDMC value(s)
  #
  #############################################################################
  #Quite often users will have a data frame called "input" already attached
  #  to the workspace. To mitigate this, we remove that if it exists, and warn
  #  the user of this case.
  if (!is.na(charmatch("input", search()))) {
    detach(input)
  }
  names(input) <- tolower(names(input))
  #order dataset
  input<-input[with(input,order(mon,day)),]
  #enable batch mode
  if (batch){
    #if id is set, then multiple stations is correct
    if ("id" %in% names(input)) {
      input <- input[with(input, order(mon, day, id)), ]
      n <- length(unique(input$id))
      if(length(unique(input[1:n, "id"])) != n){
        stop("Multiple stations have to start and end at the same dates,and 
             input data must be sorted by date/time and id")
      }
    } else {
      n <- 1
    }
  #not batch mode
  } else {
    n <- nrow(input)
  }
  #set local scope variables
  temp <- input$temp
  prec <- input$prec
  ws <- input$ws
  rh <- input$rh
  mon<-input$mon
  dmc<-input$dmc
  #set some warnings if data not setup appropriately
  if (!exists("dmc") | is.null(dmc)) 
    warning("dmc is missing!")
  if (!exists("temp") | is.null(temp)) 
    warning("temperature (temp) is missing!")
  if (!exists("prec") | is.null(prec)) 
    warning("precipitation (prec) is missing!")
  if (!exists("ws") | is.null(ws)) 
    warning("wind speed (ws) is missing!")
  if (!exists("rh") | is.null(rh)) 
    warning("relative humidity (rh) is missing!")
#   if (!exists("mon") | is.null(mon)) 
#     warning("month (mon) is missing!")
  if (length(temp)%%n != 0)
    warning("Input data do not match with number of weather stations")
  
  #Reference latitude for DMC day length adjustment
  #Using the Canadian reference only
  el <- c(6.5, 7.5, 9.0, 12.8, 13.9, 13.9, 12.4, 10.9, 9.4, 8.0, 7.0, 6.0)
  #Constrain rh, ws and precipitation
  rh <- ifelse(rh > 99.9, 99.9, rh)
  rh <- ifelse(rh < 0.0, 10.0, rh)
  ws <- ifelse(ws < 0.0, 0.0, ws)
  prec <- ifelse(prec < 0.0, 0.0, prec)
  
  #Length of weather run
  n0 <- length(temp) %/% n
  SDMC <- NULL
  #loop through all elements
  for (i in 1:n0){
    #k is the data for all stations by day
    k <- (n * (i - 1) + 1):(n * i)
    #initialize sdmc if it does not exist
    if (is.null(sdmc_old)){
      sdmc_old <- 2.6 + (1.7 * dmc[k]) 
      sdmc_old <- sdmc_old - 6.0
      sdmc_old <- ifelse(sdmc_old < 12, 12, sdmc_old)
    } 
    #Constrain temperature
    t0 <- ifelse(temp[k] < -1.1, -1.1, temp[k])     
    #This is a modification multpliper at front
    rk = 4.91 / 3.57 * 1.894 * (t0 + 1.1) * (100 - rh[k]) * el[mon[k]] * 0.0001
    #Eq.7 (Wotton et. al. 2005) calculates rain throughfall.
    rw <- ifelse(prec[k] < 7.69, 0.218 * prec[k] - 0.094, 0.83 * prec[k] - 4.8)
    #Alteration to Eq. 12 (Van Wagner & Pickett 1985)
    wmi <- 20.0 + 280.0 / exp(0.023 * sdmc_old)
    #Eqs. 13a, 13b, 13c (Van Wagner & Pickett 1985)
    b <- ifelse(sdmc_old <= 33, 100.0 / (0.5 + 0.3 * sdmc_old), 14.0 - 1.3 * 
                log(sdmc_old))
    b <- ifelse(sdmc_old > 65, 6.2 * log(sdmc_old) - 17.2, b)
    #Eq. 14 (Van Wagner & Pickett 1985) - Moisture content after rain
    wmr <- wmi + 1000.0 * rw / (48.77 + b * rw)
    #Alteration to Eq. 15 (Van Wagner & Pickett 1985)
    pr <- ifelse(prec[k] <= 0.44, sdmc_old, 43.43 * (5.6348 - log(wmr - 20)))
    #Constrain p
    pr<-ifelse(pr<0,0,pr)
    #Calculate final SDMC
    SDMC0 <- pr + rk
    #Constrain result
    SDMC0 <- ifelse(SDMC0 < 0, 0, SDMC0)
    SDMC <- c(SDMC, SDMC0)
    sdmc_old <- SDMC0
  }
  return(SDMC)
}  

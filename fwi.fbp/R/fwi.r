fwi <- function(input,yda.fwi=NULL,init=c(ffmc_yda=85,dmc_yda=6,dc_yda=15,lat=55),out="all",lat.adjust=TRUE)
{
  .Deprecated(new="fwi", 
              package="cffdrs", 
              msg="The 'fwi.fbp' package and contained functions are being deprecated and replaced by the 'cffdrs' package, please update your code to use the 'cffdrs' package.",
              old="fwi")
  
  
  ell01 <- c(6.5, 7.5, 9, 12.8, 13.9, 13.9, 12.4, 10.9, 9.4,8, 7, 6)
  ell02 <- c(7.9, 8.4, 8.9, 9.5, 9.9, 10.2, 10.1, 9.7, 9.1,8.6, 8.1, 7.8)
  ell03 <- c(10.1, 9.6, 9.1, 8.5, 8.1, 7.8, 7.9, 8.3, 8.9,9.4, 9.9, 10.2)
  ell04 <- c(11.5, 10.5, 9.2, 7.9, 6.8, 6.2, 6.5, 7.4, 8.7,10, 11.2, 11.8)
  fl01 <- c(-1.6, -1.6, -1.6, 0.9, 3.8, 5.8, 6.4, 5, 2.4, 0.4,-1.6, -1.6)
  fl02 <- c(6.4, 5, 2.4, 0.4, -1.6, -1.6, -1.6, -1.6, -1.6,0.9, 3.8, 5.8)
  names(input) <- tolower(names(input))
    ### detach the dataset when it is attached from last fail run
  if (!is.na(charmatch("input", search()))) {
    detach(input)
  }
  if (is.vector(init)) {init0 <- as.data.frame(t(init))} else {init0 <- init}
  if (!is.null(yda.fwi)) {init <- NULL}
  if (is.null(init)) {
    yda.fwi <- as.data.frame(yda.fwi)
    if ("lat" %in% names(yda.fwi)){old_lat<-yda.fwi$lat}else{old_lat<-rep(55, nrow(yda.fwi))}
    if("yr" %in% names(yda.fwi)){old_yr<-yda.fwi$yr}else{old_yr<-rep(-99, nrow(yda.fwi))}
    ffmc_yda <- yda.fwi$ffmc
    dmc_yda <- yda.fwi$dmc
    dc_yda <- yda.fwi$dc
    }  else {
    ffmc_yda <- init0[, 1]
    dmc_yda <- init0[, 2]
    dc_yda <- init0[, 3]
    old_yr <- rep(-99, nrow(init0))
    if(ncol(init0)==3){old_lat<-rep(55, nrow(init0))}else{old_lat<-init0[,4]}
  }
  if ("id" %in% names(input)){id<-input$id}else{id<-rep(-99, nrow(input))}
  if ("lat" %in% names(input)){lat<-input$lat}else{lat<-rep(55, nrow(input))}
  if ("long" %in% names(input)){long<-input$long}else{long<-rep(-120, nrow(input))}
  if ("yr" %in% names(input)){yr<-input$yr}else{yr<-rep(2011, nrow(input))}
  if ("mon" %in% names(input)){mon<-input$mon}else{mon<-rep(7, nrow(input))}
  if ("day" %in% names(input)){day<-input$day}else{day<-rep(-99, nrow(input))}
  temp <- input$temp
  prec <- input$prec
  ws <- input$ws
  rh <- input$rh
  ffmc_yda <- ifelse(lat > old_lat | old_yr > yr, init0[, 1],ffmc_yda)
  dmc_yda <- ifelse(lat > old_lat | old_yr > yr, init0[, 2],dmc_yda)
  dc_yda <- ifelse(lat > old_lat | old_yr > yr, init0[, 3],dc_yda)
  rh <- ifelse(rh >= 100, 99.9999, rh)
  # FWI computations.
    # the ffmc part
  wmo <- 147.2 * (101 - ffmc_yda)/(59.5 + ffmc_yda)
  ra <- ifelse(prec > 0.5, prec - 0.5, prec)
  wmo <- ifelse(prec > 0.5, ifelse(wmo > 150, wmo + 0.0015 * 
              (wmo - 150) * (wmo - 150) * sqrt(ra) + 42.5 * ra * exp(-100/(251 - 
              wmo)) * (1 - exp(-6.93/ra)), wmo + 42.5 * ra * exp(-100/(251 - 
              wmo)) * (1 - exp(-6.93/ra))), wmo)
  wmo <- ifelse(wmo > 250, 250, wmo)
  ed <- 0.942 * (rh^0.679) + (11 * exp((rh - 100)/10)) + 0.18 * (21.1 - temp) * (1 - 1/exp(rh * 0.115))
  ew <- 0.618 * (rh^0.753) + (10 * exp((rh - 100)/10)) + 0.18 * (21.1 - temp) * (1 - 1/exp(rh * 0.115))
  z <- ifelse(wmo < ed & wmo < ew, 0.424 * (1 - (((100 - rh)/100)^1.7)) + 0.0694 * sqrt(ws) * (1 - ((100 - rh)/100)^8), NA)
  x <- z * 0.581 * exp(0.0365 * temp)
  wm <- ifelse(wmo < ed & wmo < ew, ew - (ew - wmo)/(10^x), wmo)
  z <- ifelse(wmo > ed, 0.424 * (1 - (rh/100)^1.7) + 0.0694 * sqrt(ws) * (1 - (rh/100)^8), z)
  x <- z * 0.581 * exp(0.0365 * temp)
  wm <- ifelse(wmo > ed, ed + (wmo - ed)/(10^x), wm)
  ffmc <- (59.5 * (250 - wm))/(147.2 + wm)
  ffmc <- ifelse(ffmc > 101, 101, ffmc)
  ffmc <- ifelse(ffmc < 0, 0, ffmc)
      # the DMC part
  t0 <- temp
  t0 <- ifelse(t0 < (-1.1), -1.1, t0)
  rk <- 1.894 * (t0 + 1.1) * (100 - rh) * ell01[mon] * 1e-04
  if (lat.adjust) {
    rk <- ifelse(lat <= 30 & lat > 10, 1.894 * (t0 + 1.1) * (100 - rh) * ell02[mon] * 1e-04, rk)
    rk <- ifelse(lat <= -10 & lat > -30, 1.894 * (t0 + 1.1) * (100 - rh) * ell03[mon] * 1e-04, rk)
    rk <- ifelse(lat <= -30 & lat >= -90, 1.894 * (t0 + 1.1) * (100 - rh) * ell04[mon] * 1e-04, rk)
    rk <- ifelse(lat <= 10 & lat > -10, 1.894 * (t0 + 1.1) * (100 - rh) * 9 * 1e-04, rk)
  }
  ra <- prec
  rw <- 0.92 * ra - 1.27
  wmi <- 20 + 280/exp(0.023 * dmc_yda)
  b <- ifelse(dmc_yda <= 33, 100/(0.5 + 0.3 * dmc_yda), ifelse(dmc_yda <= 65, 14 - 1.3 * log(dmc_yda), 6.2 * log(dmc_yda) - 17.2))
  wmr <- wmi + 1000 * rw/(48.77 + b * rw)
  op <- options(warn = (-1))
  pr0 <- 43.43 * (5.6348 - log(wmr - 20))
  options(op)
  pr <- ifelse(prec <= 1.5, dmc_yda, pr0)
  pr <- ifelse(pr < 0, 0, pr)
  dmc <- pr + rk
  dmc <- ifelse(dmc < 0, 0, dmc)
    # the DC part
  t0 <- ifelse(temp < (-2.8), -2.8, t0)
  pe <- (0.36 * (t0 + 2.8) + fl01[mon])/2
  if (lat.adjust) {
    pe <- ifelse(lat <= -10, (0.36 * (t0 + 2.8) + fl02[mon])/2,pe)
    pe <- ifelse(lat > -10 & lat <= 10, (0.36 * (t0 + 2.8) + 1.4)/2, pe)
  }
  ra <- prec
  rw <- 0.83 * ra - 1.27
  smi <- 800 * exp(-1 * dc_yda/400)
  dr0 <- dc_yda - 400 * log(1 + 3.937 * rw/smi)
  dr0 <- ifelse(dr0 < 0, 0, dr0)
  dr <- ifelse(prec <= 2.8, dc_yda, dr0)
  dc <- dr + pe
  dc <- ifelse(dc < 0, 0, dc)

        #FWI equation that is consistent with C/C++
  fW <- exp(0.05039 * ws)
  fm <- 147.2 * (101 - ffmc)/(59.5 + ffmc)
  fF <- 91.9 * exp(-0.1386 * fm) * (1 + (fm^5.31)/49300000)
  isi <- 0.208 * fW * fF
  bui <- ifelse(dmc == 0 & dc == 0, 0, 0.8 * dc * dmc/(dmc + 0.4 * dc))
  p <- ifelse(dmc==0,0,(dmc-bui)/dmc)
    cc <- 0.92 + ((0.0114 * dmc)^1.7)
  bui0 <- dmc - cc * p
  bui0 <- ifelse(bui0 < 0, 0, bui0)
  bui <- ifelse(bui < dmc, bui0, bui)
  bb <- ifelse(bui > 80, 0.1 * isi * (1000/(25 + 108.64/exp(0.023 * bui))), 0.1 * isi * (0.626 * (bui^0.809) + 2))
  fwi <- ifelse(bb <= 1, bb, exp(2.72 * ((0.434 * log(bb))^0.647)))
  dsr <- 0.0272 * (fwi^1.77)
  if (out == "fwi") {
    new_FWI <- cbind(ffmc, dmc, dc, isi, bui, fwi, dsr)
    new_FWI
  } else {
    if (out == "all") {
      new_FWI <- cbind(input, ffmc, dmc, dc, isi, bui, fwi, dsr)
      new_FWI
    }
  }
}

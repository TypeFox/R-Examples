dead_reckoning <-
function(speed, heading, angle = "degree", ret = TRUE, depth = NULL,
                     pitch = NULL, endcoords=NULL, speedhorizontal = "corrected"){
  
  #exception handling
  if (angle != "degree" & angle != "radian") stop("heading must be in either degree or radian") 
  if (speedhorizontal != "corrected" & speedhorizontal != "pitch" & speedhorizontal != "depth") stop("speedhorizontal is mispecified") 
  if (length(speed) != length(heading)) stop("speed and heading must be the same length")
  
  dive_len <- length(speed)
  
  # Calculate horizontal speed values.
  if (speedhorizontal == "pitch" & is.null(pitch) == FALSE){
    speed <- cos(pitch)*speed
    print("horizontal speed is being calculated using cos(pitch)*speed")
  }
  if (speedhorizontal == "depth" & is.null(depth) == FALSE){
    speed <- cos(abs(pitch))*speed
    deltadepth <- c(0,abs(diff(depth)))
    speed <- ifelse(speed < deltadepth,deltadepth,sqrt(speed^2 - deltadepth^2))
    print("horizontal speed is being calculated using speed and change in depth")
  }
  
  # 1) Course Steered, by dead reckoning
  
  if (angle == "degree"){ xx <- speed*sin(heading/(180/pi)) ; yy <- speed*cos(heading/(180/pi)) }
  if (angle == "radian"){ xx <- speed*sin(heading) ; yy <- speed*cos(heading)  }    
  xx[1] <- 0
  yy[1] <- 0
  xs <- cumsum(xx)
  ys <- cumsum(yy)
  
  # 2) Course Made Good
  
  # 2.1 )Drift values
  # For a return dive (i.e. animal returns to the start location)
  # For a non-return dive (i.e. animal travels to a location other than the start location)
  # If the location is known (e.g. satellite tag location fix) then the "endcoords" will be used as the end location.
  # If it is a non-return dive, with an unkown end location, course made good cannot be calculated and a course steered is
  # the only option. This will result in CMGx and CMGy values "NA". 
  if (ret == TRUE){endcoords <- c(xs[1],ys[1])}
  
  driftx <- (xs[dive_len] - endcoords[1])/dive_len ; drifty <- (ys[dive_len] - endcoords[2])/dive_len 
  drift <-  sqrt((xs[dive_len] - endcoords[1])^2 + (ys[dive_len] - endcoords[2])^2)/dive_len
  errordistance <- sqrt((xs[dive_len] - endcoords[1])^2 + (ys[dive_len] - endcoords[2])^2)  
  set <- atan2((endcoords[2] - ys[dive_len]),(endcoords[1] - xs[dive_len])) 
  set <- ((2*pi) - (set-(pi/2))) %% (2*pi)


  #Provide initial values for the course made good
  xcmg <- xs[1]
  ycmg <- ys[1]
  
  #Calculate course made good
  tty <- 1:dive_len
  xcmg <- xs - (driftx*tty)
  ycmg <- ys - (drifty*tty)
  xcmg[1] <- 0
  ycmg[1] <- 0

  #Speed made good
  speedmg <- c(0,sqrt(diff(xcmg)^2 + diff(ycmg)^2))
  if (is.null(depth) == FALSE){speedmg <- c(0,sqrt(diff(xcmg)^2 + diff(ycmg)^2 + diff(depth)^2))}
  
  navlist <- list("CSx" = xs,"CSy" = ys, "CMGx" = xcmg, "CMGy" = ycmg, "speedh" = speed, "speedmg" = speedmg,
                  "drift" = drift, "errordistance" = errordistance, "set" = set)
  class(navlist) <- "navigate"
  return(navlist)
  
}


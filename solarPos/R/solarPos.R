
##  Julian day
julianDay <- function(year, month, day, hour=12, min=0, sec=0, tz=0, dut1=0){
  day_decimal = day + (hour - tz + (min + (sec + dut1)/60.0)/60.0)/24.0
  month.adj <- month + (month<3)*12
  year.adj <- year + (month<3)*(-1)
  julian_day = trunc(365.25*(year.adj+4716.0)) + trunc(30.6001*(month.adj+1)) + day_decimal - 1524.5
  B <- (2 - trunc(year.adj/100) + trunc(trunc(year.adj/100)/4))
  julian_day <- julian_day+(julian_day > 2299160.0)*B
  return(julian_day)
}

##  Solar zenith angle and azimuth angle
solarPosition <- function(jd, lon, lat, delta_t=32.184, elev=0, temp=16, pres=1013.25){
  limit_degrees <- function(degrees){
    limited <- degrees/360.0
    return(360.0*(limited-floor(limited)))
  }
  d2r <- function(degrees){
    return(degrees*pi/180)
  }
  r2d <- function(radians){
    return(radians*180/pi)
  }
  julian_ephemeris_day <- function(jd, delta_t=32.184){
    jde <- jd+delta_t/86400.0
    return(jde)
  }
  julian_century <- function(jd){
    jc <- (jd-2451545.0)/36525.0
    return(jc)
  }
  julian_ephemeris_century <- function(jde){
    jce <- (jde - 2451545.0)/36525.0
    return(jce)
  }
  julian_ephemeris_millennium <- function(jce){
    jme <- jce/10.0
    return(jme)
  }
  
  sunGeocentricPosition <- function(jme){
    EPT_L <- list(
      c(175347046.0,0,0,
        3341656.0,4.6692568,6283.07585,
        34894.0,4.6261,12566.1517,
        3497.0,2.7441,5753.3849,
        3418.0,2.8289,3.5231,
        3136.0,3.6277,77713.7715,
        2676.0,4.4181,7860.4194,
        2343.0,6.1352,3930.2097,
        1324.0,0.7425,11506.7698,
        1273.0,2.0371,529.691,
        1199.0,1.1096,1577.3435,
        990,5.233,5884.927,
        902,2.045,26.298,
        857,3.508,398.149,
        780,1.179,5223.694,
        753,2.533,5507.553,
        505,4.583,18849.228,
        492,4.205,775.523,
        357,2.92,0.067,
        317,5.849,11790.629,
        284,1.899,796.298,
        271,0.315,10977.079,
        243,0.345,5486.778,
        206,4.806,2544.314,
        205,1.869,5573.143,
        202,2.458,6069.777,
        156,0.833,213.299,
        132,3.411,2942.463,
        126,1.083,20.775,
        115,0.645,0.98,
        103,0.636,4694.003,
        102,0.976,15720.839,
        102,4.267,7.114,
        99,6.21,2146.17,
        98,0.68,155.42,
        86,5.98,161000.69,
        85,1.3,6275.96,
        85,3.67,71430.7,
        80,1.81,17260.15,
        79,3.04,12036.46,
        75,1.76,5088.63,
        74,3.5,3154.69,
        74,4.68,801.82,
        70,0.83,9437.76,
        62,3.98,8827.39,
        61,1.82,7084.9,
        57,2.78,6286.6,
        56,4.39,14143.5,
        56,3.47,6279.55,
        52,0.19,12139.55,
        52,1.33,1748.02,
        51,0.28,5856.48,
        49,0.49,1194.45,
        41,5.37,8429.24,
        41,2.4,19651.05,
        39,6.17,10447.39,
        37,6.04,10213.29,
        37,2.57,1059.38,
        36,1.71,2352.87,
        36,1.78,6812.77,
        33,0.59,17789.85,
        30,0.44,83996.85,
        30,2.74,1349.87,
        25,3.16,4690.48),
      c(628331966747.0,0,0,
        206059.0,2.678235,6283.07585,
        4303.0,2.6351,12566.1517,
        425.0,1.59,3.523,
        119.0,5.796,26.298,
        109.0,2.966,1577.344,
        93,2.59,18849.23,
        72,1.14,529.69,
        68,1.87,398.15,
        67,4.41,5507.55,
        59,2.89,5223.69,
        56,2.17,155.42,
        45,0.4,796.3,
        36,0.47,775.52,
        29,2.65,7.11,
        21,5.34,0.98,
        19,1.85,5486.78,
        19,4.97,213.3,
        17,2.99,6275.96,
        16,0.03,2544.31,
        16,1.43,2146.17,
        15,1.21,10977.08,
        12,2.83,1748.02,
        12,3.26,5088.63,
        12,5.27,1194.45,
        12,2.08,4694,
        11,0.77,553.57,
        10,1.3,6286.6,
        10,4.24,1349.87,
        9,2.7,242.73,
        9,5.64,951.72,
        8,5.3,2352.87,
        6,2.65,9437.76,
        6,4.67,4690.48),
      c(52919.0,0,0,
        8720.0,1.0721,6283.0758,
        309.0,0.867,12566.152,
        27,0.05,3.52,
        16,5.19,26.3,
        16,3.68,155.42,
        10,0.76,18849.23,
        9,2.06,77713.77,
        7,0.83,775.52,
        5,4.66,1577.34,
        4,1.03,7.11,
        4,3.44,5573.14,
        3,5.14,796.3,
        3,6.05,5507.55,
        3,1.19,242.73,
        3,6.12,529.69,
        3,0.31,398.15,
        3,2.28,553.57,
        2,4.38,5223.69,
        2,3.75,0.98),
      c(289.0,5.844,6283.076,
        35,0,0,
        17,5.49,12566.15,
        3,5.2,155.42,
        1,4.72,3.52,
        1,5.3,18849.23,
        1,5.97,242.73),
      c(114.0,3.142,0,
        8,4.13,6283.08,
        1,3.84,12566.15),
      c(1,3.14,0)
    )
    EPT_L <- lapply(EPT_L, function(x) matrix(data=x, ncol=3, byrow=TRUE))
    
    EPT_B <- list(
      c(280.0,3.199,84334.662,
        102.0,5.422,5507.553,
        80,3.88,5223.69,
        44,3.7,2352.87,
        32,4,1577.34),
      c(9,3.9,5507.55,
        6,1.73,5223.69)
    )
    EPT_B <- lapply(EPT_B, function(x) matrix(data=x, ncol=3, byrow=TRUE))
    
    EPT_R <- list(
      c(100013989.0,0,0,
        1670700.0,3.0984635,6283.07585,
        13956.0,3.05525,12566.1517,
        3084.0,5.1985,77713.7715,
        1628.0,1.1739,5753.3849,
        1576.0,2.8469,7860.4194,
        925.0,5.453,11506.77,
        542.0,4.564,3930.21,
        472.0,3.661,5884.927,
        346.0,0.964,5507.553,
        329.0,5.9,5223.694,
        307.0,0.299,5573.143,
        243.0,4.273,11790.629,
        212.0,5.847,1577.344,
        186.0,5.022,10977.079,
        175.0,3.012,18849.228,
        110.0,5.055,5486.778,
        98,0.89,6069.78,
        86,5.69,15720.84,
        86,1.27,161000.69,
        65,0.27,17260.15,
        63,0.92,529.69,
        57,2.01,83996.85,
        56,5.24,71430.7,
        49,3.25,2544.31,
        47,2.58,775.52,
        45,5.54,9437.76,
        43,6.01,6275.96,
        39,5.36,4694,
        38,2.39,8827.39,
        37,0.83,19651.05,
        37,4.9,12139.55,
        36,1.67,12036.46,
        35,1.84,2942.46,
        33,0.24,7084.9,
        32,0.18,5088.63,
        32,1.78,398.15,
        28,1.21,6286.6,
        28,1.9,6279.55,
        26,4.59,10447.39),
      c(103019.0,1.10749,6283.07585,
        1721.0,1.0644,12566.1517,
        702.0,3.142,0,
        32,1.02,18849.23,
        31,2.84,5507.55,
        25,1.32,5223.69,
        18,1.42,1577.34,
        10,5.91,10977.08,
        9,1.42,6275.96,
        9,0.27,5486.78),
      c(4359.0,5.7846,6283.0758,
        124.0,5.579,12566.152,
        12,3.14,0,
        9,3.63,77713.77,
        6,1.87,5573.14,
        3,5.47,18849.23),
      c(145.0,4.273,6283.076,
        7,3.92,12566.15),
      c(4,2.56,6283.08)
    )
    EPT_R <- lapply(EPT_R, function(x) matrix(data=x, ncol=3, byrow=TRUE))
    
    earth_heliocentric_longitude <- function(jme){
      LX <- do.call(rbind,lapply(EPT_L, function(y) rowSums(matrix(data=apply(y,1,function(x) x[1]*cos((x[2]+x[3]*jme))), nrow=length(jme)))))
      pow <- apply(matrix(jme, ncol=1), 1, function(x) x^(0:(nrow(LX)-1)))
      L <- r2d(colSums(LX*pow)/10^8)
      return(limit_degrees(L))
    }
    earth_heliocentric_latitude <- function(jme){
      BX <- do.call(rbind,lapply(EPT_B, function(y) rowSums(matrix(data=apply(y,1,function(x) x[1]*cos((x[2]+x[3]*jme))), nrow=length(jme)))))
      pow <- apply(matrix(jme, ncol=1), 1, function(x) x^(0:(nrow(BX)-1)))
      B <- r2d(colSums(BX*pow)/10^8)
      return(B)
    }
    earth_heliocentric_radius <- function(jme){
      RX <- do.call(rbind,lapply(EPT_R, function(y) rowSums(matrix(data=apply(y,1,function(x) x[1]*cos((x[2]+x[3]*jme))), nrow=length(jme)))))
      pow <- apply(matrix(jme, ncol=1), 1, function(x) x^(0:(nrow(RX)-1)))
      R <- colSums(RX*pow)/10^8
      return(R)
    }
    sun_geocentric_longitude <- function(L){
      Theta <- L+180
      return(limit_degrees(Theta))
    }
    sun_geocentric_latitude <- function(B){
      beta <- -B
      return(limit_degrees(beta))
    }
    
    L <- earth_heliocentric_longitude(jme)
    B <- earth_heliocentric_latitude(jme)
    R <- earth_heliocentric_radius(jme)
    
    Theta <- sun_geocentric_longitude(L)
    beta <- sun_geocentric_latitude(B)
    
    return(matrix(data=c(Theta, beta, R), ncol=3))
  }
  nutation <- function(jce){
    Y_TERMS <- c(
      0,0,0,0,1,
      -2,0,0,2,2,
      0,0,0,2,2,
      0,0,0,0,2,
      0,1,0,0,0,
      0,0,1,0,0,
      -2,1,0,2,2,
      0,0,0,2,1,
      0,0,1,2,2,
      -2,-1,0,2,2,
      -2,0,1,0,0,
      -2,0,0,2,1,
      0,0,-1,2,2,
      2,0,0,0,0,
      0,0,1,0,1,
      2,0,-1,2,2,
      0,0,-1,0,1,
      0,0,1,2,1,
      -2,0,2,0,0,
      0,0,-2,2,1,
      2,0,0,2,2,
      0,0,2,2,2,
      0,0,2,0,0,
      -2,0,1,2,2,
      0,0,0,2,0,
      -2,0,0,2,0,
      0,0,-1,2,1,
      0,2,0,0,0,
      2,0,-1,0,1,
      -2,2,0,2,2,
      0,1,0,0,1,
      -2,0,1,0,1,
      0,-1,0,0,1,
      0,0,2,-2,0,
      2,0,-1,2,1,
      2,0,1,2,2,
      0,1,0,2,2,
      -2,1,1,0,0,
      0,-1,0,2,2,
      2,0,0,2,1,
      2,0,1,0,0,
      -2,0,2,2,2,
      -2,0,1,2,1,
      2,0,-2,0,1,
      2,0,0,0,1,
      0,-1,1,0,0,
      -2,-1,0,2,1,
      -2,0,0,0,1,
      0,0,2,2,1,
      -2,0,2,0,1,
      -2,1,0,2,1,
      0,0,1,-2,0,
      -1,0,1,0,0,
      -2,1,0,0,0,
      1,0,0,0,0,
      0,0,1,2,0,
      0,0,-2,2,2,
      -1,-1,1,0,0,
      0,1,1,0,0,
      0,-1,1,2,2,
      2,-1,-1,2,2,
      0,0,3,2,2,
      2,-1,0,2,2)
    Y_TERMS <- matrix(data=Y_TERMS, ncol=5, byrow=TRUE)
    
    PE_TERMS <- c(
      -171996,-174.2,92025,8.9,
      -13187,-1.6,5736,-3.1,
      -2274,-0.2,977,-0.5,
      2062,0.2,-895,0.5,
      1426,-3.4,54,-0.1,
      712,0.1,-7,0,
      -517,1.2,224,-0.6,
      -386,-0.4,200,0,
      -301,0,129,-0.1,
      217,-0.5,-95,0.3,
      -158,0,0,0,
      129,0.1,-70,0,
      123,0,-53,0,
      63,0,0,0,
      63,0.1,-33,0,
      -59,0,26,0,
      -58,-0.1,32,0,
      -51,0,27,0,
      48,0,0,0,
      46,0,-24,0,
      -38,0,16,0,
      -31,0,13,0,
      29,0,0,0,
      29,0,-12,0,
      26,0,0,0,
      -22,0,0,0,
      21,0,-10,0,
      17,-0.1,0,0,
      16,0,-8,0,
      -16,0.1,7,0,
      -15,0,9,0,
      -13,0,7,0,
      -12,0,6,0,
      11,0,0,0,
      -10,0,5,0,
      -8,0,3,0,
      7,0,-3,0,
      -7,0,0,0,
      -7,0,3,0,
      -7,0,3,0,
      6,0,0,0,
      6,0,-3,0,
      6,0,-3,0,
      -6,0,3,0,
      -6,0,3,0,
      5,0,0,0,
      -5,0,3,0,
      -5,0,3,0,
      -5,0,3,0,
      4,0,0,0,
      4,0,0,0,
      4,0,0,0,
      -4,0,0,0,
      -4,0,0,0,
      -4,0,0,0,
      3,0,0,0,
      -3,0,0,0,
      -3,0,0,0,
      -3,0,0,0,
      -3,0,0,0,
      -3,0,0,0,
      -3,0,0,0,
      -3,0,0,0)
    PE_TERMS <- matrix(data=PE_TERMS, ncol=4, byrow=TRUE) 
    
    mean_elongation_moon_sun <- function(jce){
      X0 <- 297.85036 + 445267.11148*jce - 0.0019142*jce^2 + (1.0/189474.0)*jce^3
      return(X0)
    }
    mean_anomaly_sun <- function(jce){
      X1 <- 357.52772 + 35999.05034*jce - 0.0001603*jce^2 - (1.0/300000.0)*jce^3
      return(X1)
    }
    mean_anomaly_moon <- function(jce){
      X2 <- 134.96298 + 477198.867398*jce + 0.0086972*jce^2 + (1.0/56250.0)*jce^3
      return(X2)
    }
    argument_latitude_moon <- function(jce){
      X3 <- 93.27191 + 483202.017538*jce - 0.0036825*jce^2 +(1.0/327270.0)*jce^3
      return(X3)
    }
    ascending_longitude_moon <- function(jce){
      X4 <- 125.04452 - 1934.136261*jce + 0.0020708*jce^2 + (1.0/450000.0)*jce^3
      return(X4)
    }
    nutation_longitude <- function(jce, X0, X1, X2, X3, X4){
      X <- c(X0, X1, X2, X3, X4)
      f1 <- apply(PE_TERMS[,1:2], 1, function(x) x[1]+x[2]*jce)
      f2 <- apply(Y_TERMS, 1, function(Y) sin(sum(X*Y)*pi/180))
      delta_psi <- matrix(data=f1, nrow=length(jce)) %*% matrix(data=f2,ncol=1)/36000000
      return(as.numeric(delta_psi))
    }
    nutation_obliquity <- function(jce, X0, X1, X2, X3, X4){
      X <- c(X0, X1, X2, X3, X4)
      f1 <- apply(PE_TERMS[,3:4], 1, function(x) x[1]+x[2]*jce)
      f2 <- apply(Y_TERMS, 1, function(Y) cos(sum(X*Y)*pi/180))
      delta_epsilon <- matrix(data=f1, nrow=length(jce)) %*% matrix(data=f2,ncol=1)/36000000
      return(as.numeric(delta_epsilon))
    }
    
    X0 <- mean_elongation_moon_sun(jce)
    X1 <- mean_anomaly_sun(jce)
    X2 <- mean_anomaly_moon(jce)
    X3 <- argument_latitude_moon(jce)
    X4 <- ascending_longitude_moon(jce)
    delta_psi <- nutation_longitude(jce, X0, X1, X2, X3, X4)
    delta_epsilon <- nutation_obliquity(jce, X0, X1, X2, X3, X4)
    
    return(matrix(data=c(delta_psi,delta_epsilon), ncol=2))
  }
  
  ecliptic_mean_obliquity <- function(jme, delta_epsilon){
    u <- jme/10.0
    epsilon0 <- 84381.448 + -4680.93*u + -1.55*u^2 + 1999.25*u^3 -51.38*u^4 -249.67*u^5 -
      39.05*u^6 + 7.12*u^7 + 27.87*u^8 + 5.79*u^9 + 2.45*u^10
    epsilon <- epsilon0/3600 + delta_epsilon
    return(epsilon)
  }
  aberration_correction <- function(R){
    delta_tau <- -20.4898 / (3600.0*R)
    return (delta_tau)
  }
  apparent_sun_longitude <- function(Theta, delta_psi, delta_tau){
    lambda <- Theta + delta_psi + delta_tau
    return(lambda)
  }
  greenwich_sidereal_time <- function(jd, jc, delta_psi, epsilon){
    nu0 <- 280.46061837 + 360.98564736629 * (jd - 2451545.0) + 0.000387933*jc^2 - jc^3/38710000.0
    nu0 <- limit_degrees(nu0)
    nu <- nu0 + delta_psi*cos(d2r(epsilon))
    return(nu)
  }
  geocentric_right_ascension <- function(lamda, epsilon, beta){
    lamda_rad <- d2r(lamda)
    epsilon_rad <- d2r(epsilon)
    beta_rad <- d2r(beta)
    alpha <- r2d(atan2(sin(lamda_rad)*cos(epsilon_rad) -tan(beta_rad)*sin(epsilon_rad),cos(lamda_rad)))
    return(limit_degrees(alpha))
  }
  geocentric_declination <- function(beta, epsilon, lamda){
    beta_rad <- d2r(beta)
    epsilon_rad <- d2r(epsilon)
    lambda_rad <- d2r(lambda)
    delta <- r2d(asin(sin(beta_rad)*cos(epsilon_rad) + cos(beta_rad)*sin(epsilon_rad)*sin(lambda_rad)))
    return(delta)
  }
  observer_hour_angle <- function(nu, longitude, alpha){
    H <- nu + longitude - alpha
    H <- limit_degrees(H)
    return(H)
  }
  sun_equatorial_horizontal_parallax <- function(R){
    xi <- 8.794/(3600*R)
    return(xi)
  }
  parallax_right_ascension <- function(latitude, elevation, xi, H, delta){
    lat_rad <- d2r(latitude)
    xi_rad <- d2r(xi)
    H_rad <- d2r(H)
    delta_rad <- d2r(delta)
    u <- atan(0.99664719 * tan(lat_rad))
    x <- cos(u) + elevation*cos(lat_rad)/6378140.0
    y <- 0.99664719 * sin(u) + elevation*sin(lat_rad)/6378140.0
    delta_alpha <- atan2(-x*sin(xi_rad)*sin(H_rad),cos(delta_rad)-x*sin(xi_rad)*cos(H_rad))
    return(list(r2d(delta_alpha),x,y))
  }
  topocentric_right_ascension <- function(alpha, delta_alpha){
    alpha_prime <- alpha + delta_alpha
    return(alpha_prime)
  }
  topocentric_declination <- function(xi, H, delta, delta_alpha, x, y){
    xi_rad <- d2r(xi)
    H_rad <- d2r(H)
    delta_rad <- d2r(delta)
    delta_alpha_rad <- d2r(delta_alpha)
    delta_prime <- (atan2((sin(delta_rad) - y*sin(xi_rad))*cos(delta_alpha_rad),
                          cos(delta_rad) - x*sin(xi_rad) *cos(H_rad)))
    return(r2d(delta_prime))
  }
  topocentric_local_hour_angle <- function(H, delta_alpha){
    H_prime <- H - delta_alpha
    return(H_prime)
  }
  topocentric_elevation_angle <- function(latitude, delta_prime, H_prime){
    lat_rad <- d2r(latitude)
    delta_prime_rad <- d2r(delta_prime)
    H_prime_rad <- d2r(H_prime)
    e0 <- asin(sin(lat_rad)*sin(delta_prime_rad) + cos(lat_rad)*cos(delta_prime_rad)*cos(H_prime_rad))
    return(r2d(e0))
  }
  atmospheric_refraction_correction <- function(pressure, temperature, e0){
    delta_e <- (pressure/1010.0)*(283.0/(273.0 + temperature))*
      1.02/(60.0*tan(d2r(e0 + 10.3/(e0 + 5.11))))  
    #if (e0 >= -1*(SUN_RADIUS + atmos_refract))
    return(delta_e)
  }
  topocentric_elevation_angle_corrected <- function(e0, delta_e){
    e <- e0 + delta_e
    return(e)
  }
  topocentric_zenith_angle <- function(e){
    theta <- 90-e
    return(theta)
  }
  topocentric_azimuth_angle <- function(H_prime, latitude, delta_prime){
    H_prime_rad <- d2r(H_prime)
    lat_rad <- d2r(latitude)
    delta_prime_rad <- d2r(delta_prime)
    Phi <- 180+r2d(atan2(sin(H_prime_rad),cos(H_prime_rad)*sin(lat_rad) - tan(delta_prime_rad)*cos(lat_rad)))
    return(limit_degrees(Phi))
  }
  
  jde <- julian_ephemeris_day(jd,delta_t=delta_t)
  jc <- julian_century(jd)
  jce <- julian_ephemeris_century(jde)
  jme <- julian_ephemeris_millennium(jce)
  sgp <- sunGeocentricPosition(jme)
  nut <- nutation(jce)
  epsilon <- ecliptic_mean_obliquity(jme, nut[,2])
  delta_tau <- aberration_correction(sgp[,3])
  lambda <- apparent_sun_longitude(sgp[,1], nut[,1], delta_tau)
  nu <- greenwich_sidereal_time(jd, jc, nut[,1], epsilon)
  alpha <- geocentric_right_ascension(lambda, epsilon, sgp[,2])
  delta <- geocentric_declination(sgp[,2], epsilon, lambda)
  H <- observer_hour_angle(nu, lon, alpha)
  xi <- sun_equatorial_horizontal_parallax(sgp[,3])
  PRA <- parallax_right_ascension(lat, elev, xi, H, delta)
  delta_alpha <- PRA[[1]]
  x <- PRA[[2]]
  y <- PRA[[3]]
  alhpa_prime <- topocentric_right_ascension(alpha, delta_alpha)
  delta_prime <- topocentric_declination(xi, H, delta, delta_alpha,x,y)
  H_prime <- topocentric_local_hour_angle(H, delta_alpha)
  e0 <- topocentric_elevation_angle(lat, delta_prime, H_prime)
  delta_e <- atmospheric_refraction_correction(pres, temp, e0)
  e <- topocentric_elevation_angle_corrected(e0, delta_e)
  theta <- topocentric_zenith_angle(e)
  Phi <- topocentric_azimuth_angle(H_prime, lat, delta_prime)
  
  return(matrix(data=c(theta,Phi),ncol=2, dimnames=list(NULL, c("zenith","azimuth"))))
}



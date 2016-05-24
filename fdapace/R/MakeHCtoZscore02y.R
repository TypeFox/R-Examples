#' Z-score head-circumference for age 0 to 24 months based on WHO standards
#'
#' Make vector of age and height measurement to z-scores based on WHO standards using mu and sigma (not LMS)
#' 
#' @param sex A character 'M' or 'F' indicating the sex of the child. 
#' @param age A vector of time points of size Q.
#' @param hc A vector of head circumference readings of size Q (in cm).
#' 
#' @return A vector of Z-scores of size Q.
#' @export
MakeHCtoZscore02y <- function(sex, age, hc){
  time = 0:24
  
  if(length(age) != length(hc)){
    stop('Number of readings for age and head circ. are not the same.')
  }
  if(!all(  (age <= time[24]) && (time[1] <= age))){
    stop('Age requested is outside the [0,24] months.')
  }
  if(sex == 'F'){  
    # http://www.who.int/childgrowth/standards/second_set/tab_hcfa_girls_p_0_5.txt
    muGhc = c(33.8787, 36.5463, 38.2521, 39.5328, 40.5817, 41.4590, 42.1995, 42.8290, 43.3671,
              43.8300, 44.2319, 44.5844, 44.8965, 45.1752, 45.4265, 45.6551, 45.8650, 46.0598, 
              46.2424, 46.4152, 46.5801, 46.7384, 46.8913, 47.0391, 47.1822)
    sdGhc = c(1.18440, 1.17314, 1.21183, 1.24133, 1.26574, 1.28606, 1.30270, 1.31699, 1.32833, 
              1.33813, 1.34642, 1.35314, 1.35902, 1.36384, 1.36825, 1.37239, 1.37549, 1.37857, 
              1.38126, 1.38410, 1.38669, 1.38907, 1.39126, 1.39330, 1.39518);
    return(  (hc - spline(x=time, y = muGhc, xout = age)$y) / spline(x=time, y = sdGhc, xout = age)$y)
  } else if(sex == 'M'){  
    # http://www.who.int/childgrowth/standards/second_set/tab_hcfa_boys_p_0_5.txt
    muBhc = c(34.4618, 37.2759, 39.1285, 40.5135, 41.6317, 42.5576, 43.3306, 43.9803, 44.5300, 
              44.9998, 45.4051, 45.7573, 46.0661, 46.3395, 46.5844, 46.8060, 47.0088, 47.1962, 
              47.3711, 47.5357, 47.6919, 47.8408, 47.9833, 48.1201, 48.2515)
    sdBhc = c(1.27026, 1.16785, 1.17268, 1.18218, 1.19400, 1.20736, 1.22062, 1.23321, 1.24506, 
              1.25639, 1.26680, 1.27617, 1.28478, 1.29241, 1.30017, 1.30682, 1.31390, 1.32008, 
              1.32639, 1.33243, 1.33823, 1.34433, 1.34977, 1.35554, 1.36667)
    return(  (hc - spline(x=time, y = muBhc, xout = age)$y) / spline(x=time, y = sdBhc, xout = age)$y)
  } else{
    stop("Sex type undefined.")
  }
}
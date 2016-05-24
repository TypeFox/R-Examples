#' Z-score height for age 0 to 24 months based on WHO standards
#'
#' Make vector of age and height measurement to z-scores based on WHO standards using mu and sigma (not LMS)
#' 
#' @param sex A character 'M' or 'F' indicating the sex of the child. 
#' @param age A vector of time points of size Q.
#' @param ln A vector of body-length readings of size Q (in cm).
#' 
#' @return A vector of Z-scores of size Q.
#' @export
MakeLNtoZscore02y <- function(sex, age, ln){
  time = 0:24
  
  if(length(age) != length(ln)){
    stop('Number of readings for age and length are not the same.')
  }
  if(!all(  (age <= time[24]) && (time[1] <= age))){
    stop('Age requested is outside the [0,24] months.')
  }
  if(sex == 'F'){
    #http://www.who.int/childgrowth/standards/tab_lhfa_girls_p_0_2.txt
    muGln = c(49.1477, 53.6872, 57.0673, 59.8029, 62.0899, 64.0301, 
              65.7311, 67.2873, 68.7498, 70.1435, 71.4818, 72.771, 
              74.015, 75.2176, 76.3817, 77.5099, 78.6055, 79.671, 
              80.7079, 81.7182, 82.7036, 83.6654, 84.604, 85.5202, 86.4153)
    sdGln = c(1.8627, 1.9542, 2.0362, 2.1051, 2.1645, 2.2174, 
              2.2664, 2.3154, 2.365, 2.4157, 2.4676, 2.5208, 
              2.575,  2.6296, 2.6841, 2.7392, 2.7944, 2.849, 
              2.9039, 2.9582, 3.0129, 3.0672, 3.1202, 3.1737, 3.2267);
    return(  (ln - spline(x=time, y = muGln, xout = age)$y) / spline(x=time, y = sdGln, xout = age)$y)
  } else if(sex == 'M'){
    #http://www.who.int/childgrowth/standards/tab_lhfa_boys_p_0_2.txt
    muBln = c(49.8842, 54.7244, 58.4249, 61.4292, 63.886, 65.9026, 
              67.6236, 69.1645, 70.5994, 71.9687, 73.281, 74.5388, 
              75.7488, 76.9186, 78.0497, 79.1458, 80.211, 81.2487, 
              82.2587, 83.2418, 84.1996, 85.1348, 6.0477, 86.941, 87.8161)
    sdBln = c(1.8931, 1.9465, 2.0005, 2.0444, 2.0808, 2.1115, 
              2.1403, 2.1711, 2.2055, 2.2433, 2.2849, 2.3293, 
              2.3762, 2.426, 2.4773, 2.5303, 2.5844, 2.6406, 
              2.6973, 2.7553, 2.814, 2.8742, 2.9342, 2.9951, 3.0551)
    return(  (ln - spline(x=time, y = muBln, xout = age)$y) / spline(x=time, y = sdBln, xout = age)$y)
  } else{
    stop("Sex type undefined.")
  }
}
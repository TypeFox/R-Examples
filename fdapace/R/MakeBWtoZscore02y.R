#' Z-score body-weight for age 0 to 24 months based on WHO standards
#'
#' Make vector of age and body-weight to z-scores based on WHO standards using LMS
#' 
#' @param sex A character 'M' or 'F' indicating the sex of the child. 
#' @param age A vector of time points of size Q.
#' @param bw A vector of body-weight readings of size Q.
#' 
#' @return A vector of Z-scores of size Q.
#' @export
MakeBWtoZscore02y <- function(sex, age, bw){
  time = 0:24
  if(length(age) != length(bw)){
    stop('Number of readings for age and body-weight are not the same.')
  }
  if(!all(  (age <= time[24]) && (time[1] <= age))){
    stop('Age requested is outside the [0,24] months.')
  }
  
  if(sex == 'F'){
    # http://www.who.int/childgrowth/standards/tab_wfa_girls_p_0_5.txt
    sGbw = c( 0.14171, 0.13724, 0.13000, 0.12619, 0.12402, 0.12274, 0.12204, 0.12178, 0.12181, 0.12199, 0.12223,
              0.12247, 0.12268, 0.12283, 0.12294, 0.12299, 0.12303, 0.12306, 0.12309, 0.12315, 0.12323, 0.12335,
              0.12350, 0.12369, 0.12390)
    mGbw = c( 3.2322,  4.1873, 5.1282, 5.8458, 6.4237, 6.8985, 7.2970, 7.6422, 7.9487, 8.2254, 8.4800,
              8.7192,  8.9481, 9.1699, 9.3870, 9.6008, 9.8124, 10.0226, 10.2315, 10.4393, 10.6464, 10.8534,
              11.0608, 11.2688, 11.4775 )
    lGbw = c( 0.3809, 0.1714, 0.0962, 0.0402, -0.0050, -0.0430, -0.0756, -0.1039, -0.1288, -0.1507, -0.1700,
              -0.1872, -0.2024, -0.2158, -0.2278, -0.2384, -0.2478, -0.2562, -0.2637, -0.2703, -0.2762, -0.2815,
              -0.2862, -0.2903, -0.2941)
    return( ( ((bw/spline(x=time, y = mGbw, xout = age)$y)^(spline(x=time, y = lGbw, xout = age)$y)) -1 ) / 
              (spline(x=time, y = sGbw, xout = age)$y * spline(x=time, y = lGbw, xout = age)$y) )
  } else if(sex == 'M'){
    # http://www.who.int/childgrowth/standards/tab_wfa_boys_p_0_5.txt
    sBbw = c(0.14602, 0.13395, 0.12385, 0.11727, 0.11316, 0.11080, 0.10958, 0.10902, 0.10882, 0.10881, 0.10891,
             0.10906, 0.10925, 0.10949, 0.10976, 0.11007, 0.11041, 0.11079, 0.11119, 0.11164, 0.11211, 0.11261,
             0.11314, 0.11369,0.11426)
    mBbw = c(3.3464, 4.4709, 5.5675, 6.3762, 7.0023, 7.5105, 7.9340, 8.2970, 8.6151, 8.9014, 9.1649,
             9.4122, 9.6479, 9.8749, 10.0953, 10.3108, 10.5228, 10.7319, 10.9385, 11.1430, 11.3462, 11.5486,
             11.7504, 11.9514, 12.1515)
    lBbw = c(0.3487, 0.2297, 0.1970, 0.1738, 0.1553, 0.1395, 0.1257, 0.1134, 0.1021, 0.0917, 0.0820,
             0.0730, 0.0644, 0.0563, 0.0487, 0.0413, 0.0343, 0.0275, 0.0211, 0.0148, 0.0087, 0.0029,
             -0.0028, -0.0083,-0.0137)
    return( ( ((bw/spline(x=time, y = mBbw, xout = age)$y)^(spline(x=time, y = lBbw, xout = age)$y)) -1 ) / 
              (spline(x=time, y = sBbw, xout = age)$y * spline(x=time, y = lBbw, xout = age)$y) )
  } else{
    stop("Sex type undefined.")
  }
}
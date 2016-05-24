ACC <- function(var_exp, var_obs, lon = NULL, lat = NULL,
                lonlatbox = NULL, conf = TRUE, conftype = "parametric") {
  #library(abind)

  # Security checks and getting dimensions
  dimsvar <- dim(var_exp)
  if (length(dimsvar) == 5) {
    checkfirst <- 2
  }else if (length(dimsvar) == 6) {
    checkfirst <- 3
    nmembexp <- dimsvar[2]
    nmembobs <- dim(var_obs)[2]
  }else{
    stop("var_exp & var_obs should have dimensions (nexp/nsobs, nsdates, nltimes, nlat, nlon) 
                          or  dimensions (nexp/nsobs, nmembers, nsdates, nltimes, nlat, nlon) ")
  }
  for (iind in checkfirst:length(dimsvar)) {
    if (dim(var_obs)[iind] != dimsvar[iind]) {
      stop("var_exp & var_obs must have same dimensions except the first one (number of experiments or number of observational datasets) ")
    }
  }
  nexp <- dimsvar[1]
  nobs <- dim(var_obs)[1]
  nsdates <- dimsvar[checkfirst]
  nltimes <- dimsvar[checkfirst+1]
  nlat <- dimsvar[checkfirst+2]
  nlon <- dimsvar[checkfirst+3]

  # Selecting the domain
  if (is.null(lon) == FALSE & is.null(lat) == FALSE & 
      is.null(lonlatbox) == FALSE) {
    for (jind in 1:2) {
      while (lonlatbox[jind] < 0) {
        lonlatbox[jind] <- lonlatbox[jind] + 360
      }
      while (lonlatbox[jind] > 360) {
        lonlatbox[jind] <- lonlatbox[jind] - 360
      }
    }
    indlon <- which((lon >= lonlatbox[1] & lon <= lonlatbox[2]) | 
                    (lonlatbox[1] > lonlatbox[2] & (lon > lonlatbox[1] | 
                                                    lon < lonlatbox[2])))
    indlat <- which(lat >= lonlatbox[3] & lat <= lonlatbox[4])
  } else {
    indlon <- 1:nlon
    indlat <- 1:nlat
  }

  # Defining the outputs
  if(conf == TRUE) {
    ACC <- array(NA, dim = c(nexp, nobs, nsdates, nltimes, 4))
  } else {
    ACC <- array(NA, dim = c(nexp, nobs, nsdates, nltimes))
  }         
  MACCaux <- array(0, dim = c(nexp, nobs, nsdates, nltimes, 3))

  # Selecting the domain and preparing the ensemble-mean
  if (length(dimsvar) == 6) {
    var_exp <- array(var_exp[,,,,indlat, indlon], 
                     dim = c(nexp, nmembexp, nsdates, nltimes, 
                                length(indlat), length(indlon)))
    var_obs <- array(var_obs[,,,,indlat, indlon], 
                     dim = c(nobs, nmembobs, nsdates, nltimes, 
                                length(indlat), length(indlon)))
    tmp01 <- Mean1Dim(var_exp,2)
    tmp02 <- Mean1Dim(var_obs,2)
  }else{
    var_exp <- array(var_exp[,,,indlat, indlon], 
                     dim = c(nexp, nsdates, nltimes, 
                                length(indlat), length(indlon)))
    var_obs <- array(var_obs[,,,indlat, indlon], 
                     dim = c(nobs, nsdates, nltimes, 
                                length(indlat), length(indlon)))
    tmp01 <- var_exp
    tmp02 <- var_obs
  } 

  for( iobs in 1:nobs) {
    for( iexp in 1:nexp) {

      # tmp1 and tmp2 are splitted to handle NA before building tmp
      tmp1 <- array(tmp01[iexp, , , , ], dim = c(1, nsdates, nltimes,
                                    length(indlon) * length(indlat))) 
      tmp2 <- array(tmp02[iobs, , , , ], dim = c(1, nsdates, nltimes,
                                    length(indlon) * length(indlat)))

      # Variance(tmp1)should not take into account any point 
      # that is not available in tmp2 and therefore not accounted for 
      # in covariance(tmp1,tmp2) and vice-versa 
      tmp1[ is.na(tmp2) ] <- NA
      tmp2[ is.na(tmp1) ] <- NA

      tmp <- abind(tmp1, tmp2, along = 1)

      top <- apply(tmp, c(2, 3), function(x)
                      sum(x[1, ]*x[2, ], na.rm = TRUE) )

      bottom1 <- apply(tmp, c(2, 3), function(x)
                      sum(x[1, ]*x[1, ], na.rm = TRUE) )

      bottom2 <- apply(tmp, c(2, 3), function(x)
                      sum(x[2, ]*x[2, ], na.rm = TRUE) )
    
      bottom <- sqrt(bottom1 * bottom2 )
      
      ACCaux <- top / bottom

      #handle NA
      tmpallNA <- which(is.na(bottom) | bottom == 0)
      ACCaux[tmpallNA] <- NA
      
      top[tmpallNA] = NA
      bottom1[tmpallNA] = NA
      bottom2[tmpallNA] = NA

      #store the value to calculate the MACC
      MACCaux[iexp, iobs, , , 1] <- top
      MACCaux[iexp, iobs, , , 2] <- bottom1
      MACCaux[iexp, iobs, , , 3] <- bottom2
      
      if (conf == TRUE) {
        ACC[iexp, iobs, , , 2] <- ACCaux
        
        #calculate parametric confidence interval
        if (conftype == "parametric") {
        
          eno <- Mean1Dim( Eno(tmp2, 4), 1)
          t <- apply(eno, c(1, 2), 
                   function(x) qt(0.95, x - 2))
          enot <- abind(eno, t, along = 3)

          ACC[iexp, iobs, , , 4] <- apply(enot, c(1, 2), function(x)
                 sqrt((x[2] * x[2]) / ((x[2] * x[2]) + x[1] - 2)) )
        
          correno <- abind(ACCaux, eno, along = 3)
        
          ACC[iexp, iobs, , , 1] <- apply(correno, c(1, 2), function(x)
                  tanh(atanh(x[1]) + qnorm(0.975) / sqrt(x[2] - 3)) )
          ACC[iexp, iobs, , , 3] <- apply(correno, c(1, 2), function(x)
                  tanh(atanh(x[1]) + qnorm(0.025) / sqrt(x[2] - 3)) )
        }
      } else {
        ACC[iexp, iobs, , ] <- ACCaux
      }

    }
  }

  #   #na.rm should be TRUE to obtain a MACC even if a few
  #   #start dates are missing 
  topfinal <- apply(MACCaux, c(1, 2, 4), function(x) 
                                 sum(x[, 1], na.rm = TRUE) )
               
  bottomfinal <- apply(MACCaux, c(1, 2, 4), function(x) 
               sqrt(sum(x[, 2], na.rm = TRUE) * sum(x[, 3], na.rm = TRUE)))

  #to avoid that some NA are called NaN or Inf
  tmpNA <- which(is.na(bottomfinal) | bottomfinal == 0)

  MACC <- topfinal / bottomfinal 
  MACC[tmpNA] <- NA  
  
  if (conf == TRUE & conftype == "bootstrap") {
    if (length(dimsvar) != 6) {
      stop("Var_exp and var_obs must have a member dimension")
    }
    ndraw <- 100
    #create the matrix to store the random values
    ACC_draw  = array(dim=c(nexp,nobs,nsdates,nltimes,ndraw))
    MACC_draw = array(dim=c(nexp,nobs,nltimes,ndraw))

    #put the member dimension first
    var_exp <- aperm(var_exp, c(2, 1, 3, 4, 5, 6))
    var_obs <- aperm(var_obs, c(2, 1, 3, 4, 5, 6))

    for (jdraw in 1:ndraw) {

      #choose a randomly member index for each point of the matrix 
      indexp <- array(sample(nmembexp, size = (nexp*nmembexp*nsdates*nltimes), 
                             replace = TRUE), dim = c(nmembexp, nexp, nsdates, nltimes, 
                                                    length(indlat), length(indlon)) )
      indobs <- array(sample(nmembobs, size = (nobs*nmembobs*nsdates*nltimes), 
                             replace = TRUE), dim = c(nmembobs, nobs, nsdates, nltimes,
                                                    length(indlat), length(indlon)) )

      #combine maxtrix of data and random index
      varindexp <- abind(var_exp, indexp, along = 7 )
      varindobs <- abind(var_obs, indobs, along = 7 )

      #select randomly the members for each point of the matrix
      varexpdraw <- aperm( array( 
                    apply( varindexp, c(2, 3, 4, 5, 6), function(x) x[,1][x[,2]] ),
                                 dim = c(nmembexp, nexp, nsdates, nltimes, 
                                       length(indlat), length(indlon))),
                           c(2, 1, 3, 4, 5, 6)) 
      varobsdraw <- aperm( array(
                    apply( varindobs, c(2, 3, 4, 5, 6), function(x) x[,1][x[,2]] ),
                                 dim = c(nmembobs, nobs, nsdates, nltimes, 
                                       length(indlat), length(indlon))),
                           c(2, 1, 3, 4, 5, 6)) 

      #calculate the ACC of the randomized field
      tmpACC <- ACC(varexpdraw, varobsdraw, conf = FALSE)
      ACC_draw[,,,,jdraw] <- tmpACC$ACC
      MACC_draw[,,,jdraw] <- tmpACC$MACC
    }

    #calculate the confidence interval
    ACC[ , , , , 3] <- apply(ACC_draw, c(1, 2, 3, 4), function(x)
                                       quantile(x, 0.975, na.rm = TRUE))  
    ACC[ , , , , 1] <- apply(ACC_draw, c(1, 2, 3, 4), function(x)
                                       quantile(x, 0.025, na.rm = TRUE))  

    MACC <- InsertDim(MACC, 4, 3)
    MACC[ , , , 3] <- apply(MACC_draw, c(1, 2, 3), function(x)
                                       quantile(x, 0.975, na.rm = TRUE))  
    MACC[ , , , 1] <- apply(MACC_draw, c(1, 2, 3), function(x)
                                       quantile(x, 0.025, na.rm = TRUE))  
  }

  invisible(list(ACC = ACC, MACC = MACC))
}

"calcCirf" <-
function (k, x, irfpar, mirf = FALSE, measured_irf = vector(), 
    convalg = 1, shiftmea = vector(), lamb = 1, reftau = 0, 
    doublegaus = FALSE, multiplegaus = FALSE,
    streak = FALSE, streakT = 0, irffun = "gaus") 
{
   
  if (!mirf) {
    if(irffun == "gaus" || irffun == "doublegaus") {
      mu <- irfpar[1]
      tau <- irfpar[2]
      m <- convGausExp(k, x, mu, tau)
      #m <- rep(0, length(x) * length(k))
      #m <- as.matrix(.C("calcCirf", m = as.double(m), as.double(k), 
      #                  as.double(x), as.double(tau), as.double(mu),
      #                  as.integer(length(k)), 
      #                  as.integer(length(x)), PACKAGE="TIMP")$m)
      #dim(m) <- c(length(x), length(k))
      ## for streak synchronscan data get backsweep term 
      if(streak) {
        backsweep <- getStreakBacksweep(streakT, k, mu,x)
        m <- m + backsweep
      }
            if(doublegaus || irffun == "doublegaus") {
              scal <- irfpar[4] #note in multiplegaus, scal comes first
              tau <- irfpar[3]
              m2 <- convGausExp(k, x, mu, tau)
              #m2 <- rep(0, length(x) * length(k))
              #m2 <- as.matrix(.C("calcCirf", m = as.double(m2),
              #                   as.double(k), as.double(x), as.double(tau2),
              #                   as.double(mu), as.integer(length(k)),
              #                   as.integer(length(x)), PACKAGE="TIMP")$m)
              #dim(m2) <- c(length(x), length(k))
              if(streak) m2 <- m2 + backsweep
              m <- m + (scal * m2)
            }
    }
    
    if(multiplegaus || irffun == "multiplegaus") {
      mu <- irfpar[1]
      tau <- irfpar[2]
      m <- convGausExp(k, x, mu, tau)
      ## for streak synchronscan data get backsweep term 
      if(streak) {
        backsweep <- getStreakBacksweep(streakT, k, mu,x)
        m <- m + backsweep
      }   
          if(length(irfpar)>4) {
              for( i in 1:((length(irfpar)-2)/3)) {
              # scal shift tau
              scal <- irfpar[1+(i-1)*3+2]
              shift <- irfpar[2+(i-1)*3+2]
              tau <- irfpar[3+(i-1)*3+2]
              m2 <- convGausExp(k, x, mu+shift, tau)              
              if(streak) m2 <- m2 + backsweep
              m <- m + (scal * m2)
              }
           }
            
    }
    
    
    
    if(irffun == "step") {
      m <- convolveUnitStep(k, x, irfpar) 
      
    }
  }
  else {
    m <- matrix(0, nrow = length(x), ncol = length(k))
    xspace <- x[2] - x[1]
    if (length(shiftmea) != 0) {
      if (length(shiftmea) == 1) 
        lamb <- 1
      measured_irf <- .C("ShiftCurve", source = as.double(measured_irf), 
                         as.double(measured_irf),
                         as.double(shiftmea[lamb]/xspace), 
                         as.integer(length(x)), PACKAGE="TIMP")$source
    }
    for (i in 1:length(k)) {
      m[, i] <- switch(convalg, .C("Conv1", result = as.double(m[, i]),
                                   as.double(measured_irf),
                                   as.integer(length(x)), 
                                   as.double(k[i]), as.double(xspace), 
                                   PACKAGE="TIMP")$result, 
                       .C("Conv2", result = as.double(m[, i]),
                          as.double(measured_irf), 
                          as.integer(length(x)), as.double(k[i]),
                          as.double(xspace), 
                          PACKAGE="TIMP")$result, 
                       .C("Conv3", result = as.double(m[, 
                                     i]), as.double(measured_irf),
                          as.integer(length(x)), 
                          as.double(k[i]), as.double(xspace),
                          as.double(reftau),  
                          PACKAGE="TIMP")$result)
      
        }
    }
    m
}

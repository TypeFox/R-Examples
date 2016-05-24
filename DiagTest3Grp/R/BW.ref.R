###################################################
#This function calcualtes the bandwidth either from the reference rule or sheather-Jones plug in method
###################################################
BW.ref <-
function(x,method="KS-SJ")
  {
    sd0 <- sd(x,na.rm=TRUE)
    iqr <- IQR(x,na.rm=TRUE)##inter-qunatile range
    
    if(method=="KS-SJ")
      {
        ###Sheather-Jones bandwidth calling KernSmooth package
        scale0 <- ifelse(min(sd0,iqr)==0,"stdev","minim")##dpik uses by default "minim" , i.e., min(sd0,iqr), will have problem if it's 0.        
        bw <- dpik(x,scalest=scale0)
        
      }else if(method=="KS")
      {##reference rule

        nobs <- length(x)
        bw <- 1.06*min(sd0,iqr/1.34)/(nobs^0.2)
      }else stop ("wrong bandwidth method in BW.ref()!")

    return(bw)
  }


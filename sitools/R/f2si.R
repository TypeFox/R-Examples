f2si <-
function(number, unit=""){

  sifactor <-  c(1e-24, 1e-21, 1e-18, 1e-15,
                 1e-12, 1e-09, 1e-06, 1e-03,
                 1e+00,
                 1e+03, 1e+06, 1e+09, 1e+12,
                 1e+15, 1e+18, 1e+21, 1e+24)

  pre <- c(" y", " z", " a", " f",
           " p", " n", " u", " m",
           "",
           " k", " M", " G", " T",
           " P", " E", " Z", " Y")

  
  
  absolutenumber <- number * sign(number)      # findInterval should get positive numbers
                                        # we will fix the sign later again
  ix <- findInterval(absolutenumber, sifactor) # ix is a list like (1 0 0 1 1) 
  
  if (length(ix) > 0 ) {
    sistring <- paste(number/sifactor[ix], pre[ix], sep="", unit=unit)
  } else {
    sistring <- as.character(number)
  }
  return(sistring)
}

# test vector
#a <- 2.34*10^(-22:22)
#b <-  2.34*10^((-8:8)*3)
#c <-  -2.34*10^((-8:8)*3)


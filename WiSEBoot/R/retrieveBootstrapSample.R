retrieveBootstrapSample <-
function(WiSEObj){
  ##Check the WiSEObj is from WiSE.bootstrap
  if(inherits(WiSEObj, "WiSEBoot")==FALSE){
    stop("WiSEObj should be of class wise_bootstrap (i.e. generated from the WiSE.bootstrap function)")
  }
  if(WiSEObj$wavBC!="periodic"){
    stop("WiSEObj needs the periodic boundary condition to retrieve sample.")
  }

  BootSample <- array(dim=dim(WiSEObj$BootWavelet))
  waveletMatrix <- GenW(n=dim(WiSEObj$BootWavelet)[2], family=WiSEObj$wavFam, filter.number=WiSEObj$wavFil,
                        bc=WiSEObj$wavBC)
  J <- log(dim(WiSEObj$BootWavelet)[2], base=2)
  for(ss in 1:dim(WiSEObj$BootWavelet)[3]){
    for(r in 1:dim(WiSEObj$BootWavelet)[1]){
      holdCoef <- rev(WiSEObj$BootWavelet[r, , ss])
      bootWave <- holdCoef[2^J]
      holdCoef <- holdCoef[-(2^J)]
      for(j in (J-1):0){
        bootWave <- c(bootWave, rev(holdCoef[seq(1, 2^j)]))
        holdCoef <- holdCoef[-seq(1, 2^j)]
      }
      bootWave <- matrix(bootWave, ncol=1)
      if(dim(WiSEObj$BootWavelet)[3]==1){
        BootSample[r, , ss] <- (waveletMatrix %*% bootWave + WiSEObj$BootIntercept[r] + 
                                WiSEObj$BootSlope[r]*seq(1, 2^J) )
      }else{ 
        BootSample[r, , ss] <- (waveletMatrix %*% bootWave + WiSEObj$BootIntercept[r, ss] + 
                                WiSEObj$BootSlope[r, ss]*seq(1, 2^J) )
      }
    }
  }
  return(invisible(list(BootSample=BootSample)))
}

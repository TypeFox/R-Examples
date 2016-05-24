`longpstart` <-
function(NPP=6, asta="",acomp="", theday=1, hr=0 )
{
  if(missing(NPP)) { NPP=6  }
  if(missing(asta)) {  asta="TEST"  }
  if(missing(acomp)) {  acomp="T"  }

  
  if(missing(theday)) {  theday=1 }
  if(missing(hr)) {  hr=0  }
  
  print("Doing postscript")
  plname = paste(sep=".", "longfft", asta, acomp, formatC(theday, format="d", width=3,flag="0") , formatC(hr, format="d", width=2,flag="0"))
  RPMG::jpostscript(plname , P=c(14, 14) )
  par(mfrow=c(NPP,1))
  par(mai=c(0.3,0.35,0.3,0.1))
  kplot = 0
  return(kplot)
}


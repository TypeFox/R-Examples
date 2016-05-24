`longreset` <-
function(NPP, PS)
{

  if(PS==FALSE)
    {
      locator(1)
      kplot = 0
      par(mfrow=c(NPP,1))
    }
  else
    {
     
      print("Ending postscript")
      dev.off()
      kplot = 0
    }

  return(kplot)
}


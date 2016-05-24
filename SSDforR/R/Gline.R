Gline <-
  function() {
    yo<-readline("Y ordinate for your goal line  " )
abline(h=yo,col="gray",lwd=3)
u<-readline("accept line? (y/n) ")
if (u=="n")
{replayPlot(ab)}

ab<-NULL
ab<<-recordPlot()
{assign('ab', recordPlot())}

}
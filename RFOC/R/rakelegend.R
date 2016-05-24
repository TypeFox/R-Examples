rakelegend<-function(corn="topright", pal=1)
{
  #######   add a legend to a plot with the
  ######   focal colors identified  based on rake
  fleg = 1:7
  flegc = foc.color(fleg, pal=pal)
  flab = focleg(fleg)

  legend(corn, legend=flab, col=flegc, pch=19, bg="white" )

}

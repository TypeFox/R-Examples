SinMod <-
function(x, myEx, dC2)
{
  ###  function to optimize
  ###    model:  x = c(Phs, Pos, Amp, Prd)
  ####      Phase, Position, Amplitude, Period
  ###   myEx = externally defined X-values
  ###   dC2 = externally defined  observations at X

  Phs = x[1]
  Wmid  = myEx
  Pos    = x[2]
  Amp = x[3]
  Prd = x[4]

  wi = (Amp/2) * sin( (Wmid -Phs)*2*pi/Prd ) + Pos

  sum( (dC2-wi)^2 )

}


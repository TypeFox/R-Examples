`Mmorlet` <-
function(UB=-4, LB=4, N=256, plot=FALSE)
{

  ###  create a morlet function based on the matlab style routines
  if(missing(plot) ) { plot=FALSE }
  if(missing(UB) ) { UB = -4 }
  if(missing(LB) ) { LB = 4 }
  if(missing(N) ) { N = 256 }
  
###  this is the morlet function as set up by MATLAB
  out2 = seq(from=UB, to=LB, length=N)
  out1 = exp(-(out2^2)/2) * cos(5*out2)
  if(plot==TRUE)
    {
      plot(out2, out1, type='l')
    }
  invisible(list(xval=out2, morl=out1))
}


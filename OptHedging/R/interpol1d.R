interpol1d = function(x,F0,minS,maxS)
{
 
  m = length(F0)
  
  out0 = .C("interpolation1d",
            interpol = double(1),
            as.double(x),
            as.double(F0),
            as.integer(m),
            as.double(maxS),
            as.double(minS)
           
  )
  
  out = out0$interpol
  
  
}

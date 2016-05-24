bartl <-
function(m)
{
  w = 2*(0:((m-1)/2) )/(m-1)


  if(pracma::rem(m,2) == 1)
    {
      w = c(w,  w[seq(from=(m-1)/2, by=-1, to=1) ]) 
    }
  else
    {
      w = c(w, w[seq(from=m/2, by=-1, to=1)]   )
    }


  return(w)

}

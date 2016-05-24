`Gcols` <-
function(plow=10, phi=10,  N=100, pal="rainbow", mingray=0.5)
{
  ###   get a palette with the upper or lower parts replaced
  if(missing(plow)) { plow = 10 }
  if(missing(phi)) { phi = 10 }
  if(missing(N)) { N = 100 }
  if(missing(pal)) { pal = "rainbow" }
  if(missing(mingray)) { mingray=0.5 }

  
  nlow = floor(plow*N/100)
  nhi = floor(phi*N/100)
  LOW = grey(seq(from=mingray, to =1, length=nlow))
  HI  = grey(seq(from=mingray, to =1, length=nhi))
  K = N-nlow-nhi

     FUN = match.fun(pal)
     Z = FUN(K)
  #####  Z = rainbow(K)
  return(c(LOW  , Z, HI) )
}


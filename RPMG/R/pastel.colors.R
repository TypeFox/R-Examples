`pastel.colors` <-
function(num, seed=0)
{
  if(missing(seed)) { seed=0 }

  if(seed>0)
    {
      set.seed(seed)
    }

  
  v =1

  s =  (fmod( runif(num, min=0, max=262144)  , 9) + 12) / 100.0
  h = ( runif(num, min=0, max=360) )/60

  i = floor(h)
  f = h - i
  p1 =v*(1-s)
  p2 = {v*(1-(s*f))}
  p3 = {v*(1-(s*(1-f)))}
  pastel = vector()
  for(j in 1:length(i))
    {
      I = i[j]
      if(I==0) {  r =v  ;  g =p3[j] ;  b =p1[j] }
      if(I== 1) {  r =p2[j] ;  g =v  ;  b =p1[j] }
      if(I== 2) {  r =p1[j] ;  g =v  ;  b =p3[j] }
      if(I== 3) {  r =p1[j] ;  g =p2[j] ;  b =v  }
      if(I== 4) {  r =p3[j] ;  g =p1[j] ;  b =v  }
      if(I== 5) {  r =v  ;  g =p1[j] ;  b =p2[j] }
      
      ## print(c(r,b,g))
      pastel[j] = rgb(r,g,b)	
    }
  return(pastel)

}


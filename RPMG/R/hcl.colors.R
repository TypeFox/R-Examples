`hcl.colors` <-
function(n, k=260)
  {
    if(missing(k)) k = 260
    
  S1 =  hcl(h=seq(from=0, to=k, length=n) )
   return(S1)
  }


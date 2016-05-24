Likelihoodfun <-
function(p1m, p2f, p1f, n1m, n0m, n2f, n1f, n0f)
  
{

  if (n1m != 0 & n0m != 0){

  LH <- log(p1m) * n1m + log(1 - p1m) * n0m + log(p2f) * n2f + log(p1f) * n1f + log(1 - p2f - p1f) * n0f
  
  return(LH)

  }

  else{

  LH <- log(p2f) * n2f + log(p1f) * n1f + log(1 - p2f - p1f) * n0f
  
  return(LH)

  }

}

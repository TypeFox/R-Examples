betarand <-
function(a, L0, L, M){
  if(a[1] < 0 | a[2] < 0){
    r <- rbeta(M, 1, 1)
    tr <- round(L0 + (L - L0) * r)    
  } else {
    r <- rbeta(M, a[1], a[2])
    tr <- round(L0 + (L - L0) * r)
  }
}

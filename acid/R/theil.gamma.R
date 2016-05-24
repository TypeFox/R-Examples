theil.gamma <-
function(p){ # exact from Kleiber p.164
  Theil<-1/p+digamma(p)-log(p)
  return(Theil)
}

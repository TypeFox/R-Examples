ecospat.binary.model <-function(Pred, Sp.occ.xy, Percentage)
{
  # Calculate threshold
  #####################
  Pres <- extract (Pred, Sp.occ.xy)
  Pres.sort <- data.frame(sort(Pres))
  Num.pres <- nrow(Sp.occ.xy)
  Num <- as.integer(((Num.pres * Percentage)/100) + 0.5)
  #ifelse (Num == 1, Num <- 2, Num)
  Threshold <- sort(Pres) [Num]


  # Generation of binary model
  ############################
  Threshold.Com <- c(0, Threshold, 0, Threshold, 1000, 1)
  Threshold.Com.b <- matrix (Threshold.Com, ncol=3, byrow=T)
  Pred.binary <- reclassify (Pred,Threshold.Com.b)
}

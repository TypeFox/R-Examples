promat <-
function(pc, nage)
{
  #"pc" is a element comtaining the projection matrix components estimated from the
  #function "procomp"
  lcJ <- diag(pc[, 3], nrow = nage, ncol = nage)
  if(nage == 1)
  {
    l <- matrix(rbind(pc[, 2], lcJ), ncol = nage)
  }
  else
  {
    l <- matrix(rbind(pc[, 2], lcJ[-nage, ]), ncol = nage)
    l[nage, nage] <- lcJ[nage, nage]
  }
  #Output
  l
}

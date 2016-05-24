twoplpart <-
function(Km=Km, abpar=abpar)
{
  SO <- exp(Km %*% abpar) / (1 + exp(Km %*% abpar))  
  SO
}

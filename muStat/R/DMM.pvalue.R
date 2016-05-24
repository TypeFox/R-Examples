`DMM.pvalue` <-
function(pP,qP, xX, condvar = TRUE, exact = FALSE)
{
  if(condvar) 
    MCN.pvalue(pP,      qP,      exact)                    # Wittkowski
  else
    MCN.pvalue(pP+xX/2, qP+xX/2, FALSE)                    #     (1989)
}

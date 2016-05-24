unifres.fun <-
function(posEH)
{
expres <-  c(posEH[1], diff(posEH))
expres<-expres-0.99*min(expres)
unires <- exp( -expres)
dist0<-sum(expres==0)
if (dist0>0) cat("Number of   exponential residuals =0: ", dist0, fill = TRUE)
return(list(expres=expres, unires=unires, posEH=posEH))
}

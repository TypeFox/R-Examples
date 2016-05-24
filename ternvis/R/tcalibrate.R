tcalibrate <-
function(tv, # an object of class tverify as produced by tgetcal()
         p)  # a ternary forecast to be calibrated
{
   out <- tv$f(tv$pars,p)
   out
}

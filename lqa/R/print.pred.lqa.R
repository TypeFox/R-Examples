print.pred.lqa <-
function (x, ...)
{
   cat ("\nPredicted response(s): \n")
   print (x$mu.new)
   cat ("\n")
   invisible (x)
}


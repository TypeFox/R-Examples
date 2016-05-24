gcv.loss <-
function (pred.obj)
{
   dev <- pred.obj$deviance
   tr.H <- pred.obj$tr.H
   nobs <- pred.obj$n.newobs

   nobs * dev / ((nobs - tr.H)^2)
}


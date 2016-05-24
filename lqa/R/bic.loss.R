bic.loss <-
function (pred.obj)
{
   dev <- pred.obj$deviance
   tr.H <- pred.obj$tr.H
   nobs <- pred.obj$n.newobs

   dev + log (nobs) * tr.H
}


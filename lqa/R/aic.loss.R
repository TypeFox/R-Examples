aic.loss <-
function (pred.obj)
{
## Bemerkung: Das müssen wir so umständlich machen, da die hat-matrix bei Boosting-Verfahren anders berechnet wird,
## daher können wir nicht family$aic verwenden...

   dev <- pred.obj$deviance
   tr.H <- pred.obj$tr.H

   dev + 2 * tr.H
}


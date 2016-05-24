`PRESS.lm` <-
function(object, ...)sum((object$residuals/(1-hatvalues(object)))^2)


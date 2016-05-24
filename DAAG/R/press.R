"press" <-
function(obj){sum((resid(obj)/(1-hatvalues(obj)))^2)}


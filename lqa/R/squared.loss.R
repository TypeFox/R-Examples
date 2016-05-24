squared.loss <-
function (pred.obj)
   sum ((pred.obj$new.y - pred.obj$mu.new)^2)


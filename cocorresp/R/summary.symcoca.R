"summary.symcoca" <-
function(object, axes = c(1:min(6, object$n.axes)),
         display = c("species", "sites"),
         scaling = 1, ...)
  {
    cocaScores <- scores(object, choices = axes,
                         scaling = scaling, display = display)
    retval <- list(cocaScores = cocaScores, inertia = object$inertia,
                   lambda = object$lambda,
                   call = object$call, namY = object$nam.dat$namY,
                   namX = object$nam.dat$namX, scaling = scaling)
    class(retval) <- "summary.symcoca"
    retval
  }


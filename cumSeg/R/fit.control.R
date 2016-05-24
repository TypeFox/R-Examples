fit.control <-
function (toll = 1e-03, it.max = 5, display = FALSE, last = TRUE,
    maxit.glm = 25, h = 1, stop.if.error=FALSE){
    list(toll = toll, it.max = it.max, visual = display, last = last,
        maxit.glm = maxit.glm, h = h, stop.if.error = stop.if.error)
    }


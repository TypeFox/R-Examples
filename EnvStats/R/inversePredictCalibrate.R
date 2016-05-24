inversePredictCalibrate <-
function (object, obs.y = NULL, n.points = ifelse(is.null(obs.y), 
    100, length(obs.y)), intervals = FALSE, coverage = 0.99, 
    simultaneous = FALSE, individual = FALSE, trace = FALSE) 
{
    if (!inherits(object, "calibrate")) 
        stop("The argument 'object' must inherit from the class 'calibrate'.")
    if (is.null(object$x)) 
        stop(paste("The argument 'object' must contain", "a component called 'x' that contains the model matrix.", 
            " That is, it must be the result of calling the lm", 
            "funciton with x=TRUE"))
    fcn.to.min.for.x <- function(x.weird, y.weird, object.weird, 
        x.name.weird) {
        dum.df <- data.frame(x.weird)
        names(dum.df) <- x.name.weird
        sum((y.weird - predict(object.weird, newdata = dum.df))^2)
    }
    if (is.null(obs.y)) {
        range.y <- range(object$fitted.values)
        obs.y <- seq(range.y[1], range.y[2], by = diff(range.y)/n.points)[-(n.points + 
            1)]
    }
    x.name <- attr(object$terms, "term.labels")[1]
    mean.x <- mean(object$x[, 2])
    pred.x <- numeric(n.points)
    for (i in 1:n.points) {
        if (trace) 
            cat("\t\tInverse prediction for point", i, "out of", 
                n.points, "\n")
        pred.x[i] <- nlminb(start = mean.x, objective = fcn.to.min.for.x, 
            lower = 0, control = list(x.tol = .Machine$double.eps * 
                1e+08), y.weird = obs.y[i], object.weird = object, 
            x.name.weird = x.name)$par
    }
    ret.mat <- cbind(obs.y = obs.y, pred.x = pred.x)
    if (intervals) {
        ll.x <- ul.x <- numeric(n.points)
        fcn.to.min.for.ll <- function(x.weird, y.weird, object.weird, 
            x.name.weird, coverage.weird, simultaneous.weird, 
            individual.weird) {
            dum.df <- data.frame(x.weird)
            names(dum.df) <- x.name.weird
            pred.list <- predict(object.weird, newdata = dum.df, 
                se.fit = TRUE)
            sum((y.weird - pointwise(pred.list, coverage = coverage.weird, 
                simultaneous = simultaneous.weird, individual = individual.weird)$upper)^2)
        }
        for (i in 1:n.points) {
            if (trace) 
                cat("\t\tComputing lower limit for point", i, 
                  "out of", n.points, "\n")
            ll.x[i] <- nlminb(start = pred.x[i], objective = fcn.to.min.for.ll, 
                lower = 0, upper = pred.x[i], control = list(x.tol = .Machine$double.eps * 
                  1e+08), y.weird = obs.y[i], object.weird = object, 
                x.name.weird = x.name, coverage.weird = coverage, 
                simultaneous.weird = simultaneous, individual.weird = individual)$par
        }
        fcn.to.min.for.ul <- function(x.weird, y.weird, object.weird, 
            x.name.weird, coverage.weird, simultaneous.weird, 
            individual.weird) {
            dum.df <- data.frame(x.weird)
            names(dum.df) <- x.name.weird
            pred.list <- predict(object.weird, newdata = dum.df, 
                se.fit = TRUE)
            sum((y.weird - pointwise(pred.list, coverage = coverage.weird, 
                simultaneous = simultaneous.weird, individual = individual.weird)$lower)^2)
        }
        for (i in 1:n.points) {
            if (trace) 
                cat("\t\tComputing upper limit for point", i, 
                  "out of", n.points, "\n")
            ul.x[i] <- nlminb(start = pred.x[i], objective = fcn.to.min.for.ul, 
                lower = pred.x[i], control = list(x.tol = .Machine$double.eps * 
                  1e+08), y.weird = obs.y[i], object.weird = object, 
                x.name.weird = x.name, coverage.weird = coverage, 
                simultaneous.weird = simultaneous, individual.weird = individual)$par
        }
        int.mat <- cbind(ll.x, ul.x)
        if (individual) 
            dimnames(int.mat) <- list(NULL, c("lpl.x", "upl.x"))
        else dimnames(int.mat) <- list(NULL, c("lcl.x", "ucl.x"))
        ret.mat <- cbind(ret.mat, int.mat)
        attr(ret.mat, "coverage") <- coverage
        attr(ret.mat, "simultaneous") <- simultaneous
    }
    ret.mat
}

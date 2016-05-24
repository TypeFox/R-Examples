`pairs.profile.brglm` <-
function (x, colours = 2:3, ...) 
{
    if (is.null(x$profilesBR)) {
        pairs(x$profilesML, colours = colours,
               title = "Ordinary deviance", ...)
    }
    else {
        pairs(x$profilesML, colours = colours,
               title = "Ordinary deviance", ...)
        getOption("device")()
        fit <- x$profilesBR$fit
        tt <- if (fit$pl | all(fit$family$link == "logit"))
              "Penalized deviance"
              else 
              "Modified score statistic"
        pairs(x$profilesBR, colours = colours,
               title = tt, ...)
    }
}

`plot.profile.brglm` <-
function (x, signed = FALSE, interpolate = TRUE,
          n.interpolations = 100, print.grid.points = FALSE, ...) 
{
    if (is.null(x$profilesBR)) {
        plot(x$profilesML,
             cis = NULL, signed = signed,
             interpolate = interpolate,
             n.interpolations = n.interpolations,
             print.grid.points = print.grid.points, 
             title = "Ordinary deviance", ...)
    }
    else {
        plot(x$profilesML,
             cis = NULL, signed = signed,
             interpolate = interpolate, 
             n.interpolations = n.interpolations,
             print.grid.points = print.grid.points,
             title = "Ordinary deviance", ...)
        getOption("device")()
        fit <- x$profilesBR$fit
        tt <- if (fit$pl | all(fit$family$link == "logit"))
              "Penalized deviance"
              else
              "Modified score statistic"
        plot(x$profilesBR,
             cis = NULL, signed = signed,
             interpolate = interpolate, 
             n.interpolations = n.interpolations,
             print.grid.points = print.grid.points, 
             title = tt, ...)
    }
}

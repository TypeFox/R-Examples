"summary.relimplmbootMI" <-
function (object, ..., lev = max(object@level))
{
       if (!(is(object, "relimplmbootMI")))
        stop("obejct must be the output from function mianalyze.relimp")
       if (!(lev >= 0.5 && lev < 1)) 
            stop("invalid confidence levels requested: ", paste(lev, 
                collapse = " "))
       if (is.null(object@MIresult)) stop("Slot MIresult is empty in object, no summary is created.")
       summary(object@MIresult, ..., alpha = 1-lev)
}
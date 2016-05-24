'harmonize.x' <-
function(svy.A, svy.B, form.x, x.tot=NULL, cal.method="linear", ...)
{

# this function harmonizes joint/marginal distribution of the
# common variables in sample survey A and sample survey B according
# to calibration steps introduced in Renssen's 1998 paper.

    data.A <- svy.A$variables
    n.A <- nrow(data.A)
    w0.A <- weights(svy.A)

    data.B <- svy.B$variables
    n.B <- nrow(data.B)
    w0.B <- weights(svy.B)

#    require(survey)

# if population totals are unknown (x.tot=NULL) a 'pooled' estimate id derived

    if(is.null(x.tot)){
        if(cal.method=="poststratify"){
            ff.xA <- paste("w0.A", paste(as.character(form.x), collapse=""))
            ff.xB <- paste("w0.B", paste(as.character(form.x), collapse=""))
            x.tot.A <- xtabs(as.formula(ff.xA), data=data.A)
            x.tot.B <- xtabs(as.formula(ff.xB), data=data.B)
        }
        else{
            x.tot.A <- colSums( model.matrix(form.x, data=data.A) * w0.A)
            x.tot.B <- colSums( model.matrix(form.x, data=data.B) * w0.B)
        }
        lamda <- n.A/(n.A+n.B)
        x.tot <- lamda * x.tot.A + (1-lamda) * x.tot.B
    }

# calibration/poststratification

    if(cal.method=="linear"){
        cal.A <- calibrate(design=svy.A, formula=form.x,
                                population=x.tot, calfun="linear", ...)
        cal.B <- calibrate(design=svy.B, formula=form.x,
                                population=x.tot, calfun="linear", ...)
    }
    if(cal.method=="raking"){
        cal.A <- calibrate(design=svy.A, formula=form.x,
                                population=x.tot, calfun="raking", ...)
        cal.B <- calibrate(design=svy.B, formula=form.x,
                                population=x.tot, calfun="raking", ...)
    }
    if(cal.method=="poststratify"){
        cal.A <- postStratify(design=svy.A, strata=form.x, population=x.tot, partial = TRUE)
        cal.B <- postStratify(design=svy.B, strata=form.x, population=x.tot, partial = TRUE)
    }

    list(cal.A=cal.A, cal.B=cal.B,
            weights.A=weights(cal.A), weights.B=weights(cal.B), call=match.call())
}

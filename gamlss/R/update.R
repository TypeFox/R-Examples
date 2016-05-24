#####################################################################################
#                         update.gamlss
######################################################################################
# MS Monday, October 25, 2004 
update.gamlss <- function (object, 
                          formula., 
                          ..., 
                          what = c("mu", "sigma", "nu", "tau", "All"),
                          parameter= NULL,
                          evaluate = TRUE) 
{
    call <- object$call
    if (is.null(call)) 
        stop("need an object with call component")
    extras <- match.call(expand.dots = FALSE)$...
    if (!missing(formula.)) 
        {
      what <- if (!is.null(parameter))  {
        match.arg(parameter, choices=c("mu", "sigma", "nu", "tau"))} else  match.arg(what)
        if (what=="mu") 
         { call$formula <- update.formula(formula(object,what), formula.) }
        else if (what=="sigma"||what=="nu"||what=="tau")  
         {
          call[[paste(what,"formula",sep=".")]] <- 
         if (length(update.formula(formula(object,what), formula.))==2)
          update.formula(formula(object,what), formula.)
         else
          update.formula(formula(object,what), formula.)[-2]
         }
        else 
         {
            call$formula <- update.formula(formula(object, "mu"), formula.)
            # if sigma in the model
            if (("sigma"%in%object$parameters))
            {
            call[[paste("sigma", "formula", sep = ".")]] <- if (length(update.formula(formula(object, "sigma"), formula.)) == 2)         
            update.formula(formula(object, "sigma"), formula.)
            else update.formula(formula(object, "sigma"), formula.)[-2]
            }
            # if nu
             if (("nu"%in%object$parameters))
            {
            call[[paste("nu", "formula", sep = ".")]] <- if (length(update.formula(formula(object, "nu"), formula.)) == 2)         
            update.formula(formula(object, "nu"), formula.)
            else update.formula(formula(object, "nu"), formula.)[-2]
            }
            # if tau
             if (("tau"%in%object$parameters))
            {
            call[[paste("tau", "formula", sep = ".")]] <- if (length(update.formula(formula(object, "tau"), formula.)) == 2)         
            update.formula(formula(object, "tau"), formula.)
            else update.formula(formula(object, "tau"), formula.)[-2]
            }
          }
        }
    if (length(extras) > 0) 
        {
        existing <- !is.na(match(names(extras), names(call)))
        for (a in names(extras)[existing]) call[[a]] <- extras[[a]]
        if (any(!existing)) 
           {
            call <- c(as.list(call), extras[!existing])
            call <- as.call(call)
           }
        }
    if (evaluate) 
        eval(call, parent.frame())
    else call
}
######################################################################################

################################################################################
### Part of the surveillance package, http://surveillance.r-forge.r-project.org
### Free software under the terms of the GNU General Public License, version 2,
### a copy of which is available at http://www.r-project.org/Licenses/.
###
### Functions and methods to make step() work for twinstim objects
### (restricted to one component at a time)
###
### Copyright (C) 2013 Sebastian Meyer
### $Revision: 645 $
### $Date: 2013-09-08 15:17:42 +0200 (Son, 08 Sep 2013) $
################################################################################


### To make step() work, we are dealing with modified twinstim objects:
### object$formula is replaced by the result of terms(object), which selects only
### one of the two components! The original full formula specification is
### retained in the new "formulae" component.
### We let this special class inherit from "twinstim" such that, e.g.,
### extractAIC.twinstim is used for its objects. However, this is tricky since
### the classes are actually incompatible in the formula specification. Only
### methods which don't use the $formula part work, but this constraint holds
### for what is needed to run step(), if we define some additional specific
### methods for this class.

twinstim_stependemic <- twinstim_stepepidemic <- function (object)
{
    stepClass <- grep("twinstim_step", sys.call()[[1L]], value=TRUE)
    ##<- since sys.call()[[1L]] may also be surveillance:::...
    if (identical(class(object), "twinstim")) {
        component <- sub("twinstim_step", "", stepClass)
        object$formulae <- object$formula
        object$formula <- object$formulae[[component]]
        class(object) <- c(stepClass, "twinstim")
    } else if (!inherits(object, stepClass)) stop("unintended use")
    object
}


## In the first step() loop, object$call$formula is set to terms(object). Since
## there is no "formula" argument to twinstim(), we must remove it from the call
## before update()ing. We also have to convert object$formula to the complete
## formula specification (a named list) and remove the original one ($formulae).
.step2twinstim <- function (object)
{
    ##if (identical(class(object), "twinstim")) return(object)
    component <- sub("^twinstim_step", "", class(object)[1])
    stopifnot(component %in% c("endemic", "epidemic"))
    object$call$formula <- NULL
    object$formula <- modifyList(
        object$formulae,
        setNames(list(formula(object$formula)), component)
        )
    object$formulae <- NULL
    class(object) <- "twinstim"
    object
}


### special update- and terms-methods for use through stepComponent() below

update.twinstim_stependemic <- function (object, endemic, ..., evaluate = TRUE)
{
    object <- .step2twinstim(object)
    res <- NextMethod("update")         # use update.twinstim()
    
    ## we need to keep the special class such that step() will keep invoking
    ## the special update- and terms-methods on the result
    stepClass <- sub("update.", "", .Method, fixed=TRUE)
    ##<- or: .Class[1L], or: grep("step", class(object), value=TRUE)
    if (evaluate) {
        do.call(stepClass, alist(res))
    } else {
        as.call(list(call(":::", as.name("surveillance"), as.name(stepClass)),
                     res))
        ## the call will only be evaluated within stats:::drop1.default() or
        ## stats:::add1.default, where the "stepClass" constructor function
        ## (twinstim_stependemic or twinstim_stepepidemic) is not visible;
        ## we thus have to use ":::".
    }
}

update.twinstim_stepepidemic <- function (object, epidemic, ..., evaluate = TRUE)
{}
body(update.twinstim_stepepidemic) <- body(update.twinstim_stependemic)


terms.twinstim_stependemic <- terms.twinstim_stepepidemic <-
    function (x, ...) terms(x$formula)



### Function to perform AIC-based model selection (component-specific)
### This is essentially a wrapper around stats::step()

stepComponent <- function (object, component = c("endemic", "epidemic"),
                           scope = list(upper=object$formula[[component]]),
                           direction = "both", trace = 2, verbose = FALSE, ...)
{
    component <- match.arg(component)

    ## Convert to special twinstim class where $formula is the component formula
    object_step <- do.call(paste0("twinstim_step", component), alist(object))
    
    ## silent optimizations
    if (trace <= 2) object_step$call$optim.args$control$trace <-
        object_step$optim.args$control$trace <- 0
    object_step$call$verbose <- verbose

    ## Run the selection procedure
    res <- step(object_step, scope = scope, direction = direction,
                trace = trace, ...)
    
    ## Restore original trace and verbose arguments
    if (trace <= 2) {
        res$call$optim.args$control <- object$call$optim.args$control
        res$optim.args$control <- object$optim.args$control
    }
    res$call$verbose <- object$call$verbose

    ## Convert back to original class
    .step2twinstim(res)
}


### add1.default and drop1.default work without problems through the above
### implementation of stepComponent() using the tricky twinstim classes,
### where object$formula is replaced by the requested component's formula.
### However, for stand-alone use of add1 and drop1, we need specialised methods.

add1.twinstim <- drop1.twinstim <-
    function (object, scope,
              component = c("endemic", "epidemic"),
              trace = 2, ...)
{
    component <- match.arg(component)
    
    ## Convert to special twinstim class where $formula is the component formula
    object <- do.call(paste0("twinstim_step", component), alist(object))

    ## Call the default method (unfortunately not exported from stats)
    ## Note that the next method chosen is "unchanged if the class of the
    ## dispatching argument is changed" (see ?NextMethod)
    ## (the "component" argument will be part of "..." and passed further on to
    ## extractAIC.twinstim() where it is unused)
    NextMethod(trace=trace)
}

add1.twinstim_stependemic <- drop1.twinstim_stependemic <-
    function (object, scope, ...) NextMethod(component="endemic")

add1.twinstim_stepepidemic <- drop1.twinstim_stepepidemic <-
    function (object, scope, ...) NextMethod(component="epidemic")

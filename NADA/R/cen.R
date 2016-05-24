#-->> BEGIN Support routines for survival-analysis based functions 

## Generics

setGeneric("flip", function(x) standardGeneric("flip"))

setGeneric("asSurv", function(x) standardGeneric("asSurv"))

## Classes 

setOldClass("Surv")

setClass("Cen", representation(Surv="Surv", flipFactor="numeric"))

## Methods

# Cen() is analgous to Surv() in survival package.
# However, Cen() can be used in the flip() routines below to
# rescale left-censored data into right-censored data 
# for use in the survival package routines.
Cen =
function(obs, censored, type="left")
{
    ff = max(obs)+(diff(range(obs))/2) # arbitrary flip factor > max(obs)
    new("Cen", Surv=Surv(obs, !censored, type=type), flipFactor = ff)
}

setMethod("print", signature(x="Cen"), function(x, ...) print(x@Surv))

## flip() methods for rescaling input for use in survival routines.

# Note that flip()ing a Cen object results in a Surv object.
setMethod("flip", signature(x="Cen"), function(x)
{
    surv = x@Surv
    # Subtract the flipFactor from all observations
    surv[,1] = x@flipFactor - surv[,1]

    # Update the censoring type
    if (attr(surv, "type") == "right") attr(surv, "type") = "left"
    else if (attr(surv, "type") == "left") attr(surv, "type") = "right"
    else stop("Can only flip \"left\" or \"right\" censored Cen objects")

    return(surv)
})

# flip()ing a formula just symbolically updates the 
# response (which should be a Cen object). 
# Result is like: flip(Cen(obs, cen))~groups

setMethod("flip", signature(x="formula"), function(x) update(x, flip(.)~.))

setMethod("asSurv", signature(x="formula"), function(x) update(x, asSurv(.)~.))

setMethod("asSurv", signature(x="Cen"), function(x) x@Surv)

## Begin cencen.* functions
# These routines are allow cenfit, cendiff, and cenmle methods
# to be used with Cen objects, vectors, or formulas as input.
# Note they all convert input to a formula and call the formula method.

cencen.Cen =
function(obs, censored, groups, ...)
{
    cl = match.call()
    f = as.formula(substitute(a~1, list(a=cl[[2]])))
    environment(f) = parent.frame()
    callGeneric(f, ...)
}

cencen.vectors =
function(obs, censored, groups, ...)
{
    cl = match.call()
    f = as.formula(substitute(Cen(a, b)~1, list(a=cl[[2]], b=cl[[3]])))
    environment(f) = parent.frame()
    callGeneric(f, ...)
}

cencen.vectors.groups =
function(obs, censored, groups, ...)
{
    cl = match.call()
    f = substitute(Cen(a, b)~g, list(a=cl[[2]], b=cl[[3]], g=cl[[4]]))
    f = as.formula(f)
    environment(f) = parent.frame()
    callGeneric(f, ...)
}
## End cencen.* routines

## End utility functions

#-->> END Support routines for survival-analysis based functions 

## area of circle
s = nint_space(
    nint_intvDim(-1, 1),
    nint_funcDim(function(x) nint_intvDim(c(-1, 1) * sin(acos(x[1])) ))
)
nint_integrate(function(x) 1, s) # pi
## see nint_integrate's examples for more sophisticated integrals


## prepare for custom recursive implementation
using = TRUE
nfunc = nint_integrateNFunc_recursive(
    function(f, lowerLimit, upperLimit, ...) {
        if (using) { # this function is called many times
            using <<- FALSE
            cat('using integrateA\n')
        }
        integrateA(f, lowerLimit, upperLimit, ..., subdivisions=1)$value
    }
)
unlockBinding('nint_integrateNFunc', environment(nint_integrate))
assign('nint_integrateNFunc', nfunc, envir=environment(nint_integrate))

## integrate with custom recursive implementation
nint_integrate(function(x) 1, s) # pi


## prepare for custom solution
f = function(f, funcs, x0, i0, ...) {
    # add sophisticated code here
    print(list(f=f, funcs=funcs, x0=x0, i0=i0, ...))
    stop('do something')
}
unlockBinding('nint_integrateNFunc', environment(nint_integrate))
assign('nint_integrateNFunc', f, envir=environment(nint_integrate))

## integrate with custom solution
try(nint_integrate(function(x) 1, s))

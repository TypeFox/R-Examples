## integrate with defaults (stats::integrate)
nint_integrate(sin, nint_space(nint_intvDim(pi/4, 3*pi/4)))


## prepare for integrateA
ncube = function(f, lowerLimit, upperLimit, ...) {
    cat('using integrateA\n')
    integrateA(f, lowerLimit, upperLimit, ..., subdivisions=2)
}
ncube = nint_integrateNCube_integrate(ncube)
unlockBinding('nint_integrateNCube', environment(nint_integrate))
assign('nint_integrateNCube', ncube, envir=environment(nint_integrate))

## integrate with integrateA
nint_integrate(sin, nint_space(nint_intvDim(pi/4, 3*pi/4)))


## prepare for cubature
ncube = function(f, lowerLimit, upperLimit, ...) {
    cat('using cubature\n')
    r = cubature::adaptIntegrate(f, lowerLimit, upperLimit, ..., maxEval=1e3)
    return(r$integral)
}
unlockBinding('nint_integrateNCube', environment(nint_integrate))
assign('nint_integrateNCube', ncube, envir=environment(nint_integrate))

## integrate with cubature
nint_integrate(sin, nint_space(nint_intvDim(pi/4, 3*pi/4)))


## prepare for SparseGrid
ncube = function(dimension) {
    cat('using SparseGrid\n')
    SparseGrid::createIntegrationGrid('GQU', dimension, 7)
}
ncube = nint_integrateNCube_SparseGrid(ncube)
unlockBinding('nint_integrateNCube', environment(nint_integrate))
assign('nint_integrateNCube', ncube, envir=environment(nint_integrate))

## integrate with SparseGrid
nint_integrate(sin, nint_space(nint_intvDim(pi/4, 3*pi/4)))

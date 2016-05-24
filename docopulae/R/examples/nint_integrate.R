## discrete
## a) scatter
s = nint_space(nint_scatDim(1:3),
               nint_scatDim(c(0, 2, 5)))
s
## (1, 0), (2, 2), (3, 5)
nint_integrate(function(x) abs(x[1] - x[2]), s) # 1 + 0 + 2 == 3

## b) grid
s = nint_space(nint_gridDim(1:3),
               nint_gridDim(c(0, 2, 5)))
s
## (1, 0), (1, 2), (1, 5), (2, 0), ..., (3, 2), (3, 5)
nint_integrate(function(x) ifelse(sum(x) < 5, 1, 0), s) # 5


## continous
## c)
s = nint_space(nint_intvDim(1, 3),
               nint_intvDim(1, Inf))
s
nint_integrate(function(x) 1/x[2]**2, s) # 2

## d) infinite, no transform
s = nint_space(nint_intvDim(-Inf, Inf))
nint_integrate(sin, s) # 0

## e) infinite, transform
s = nint_space(nint_intvDim(-Inf, Inf),
               nint_intvDim(-Inf, Inf))
## probability integral transform
tt = nint_transform(function(x) prod(dnorm(x)), s, 1:2, list(
    g=function(x) pnorm(x),
    gij=function(y) { t1 = qnorm(y); cbind(t1, 1/dnorm(t1)) }))
tt$space
nint_integrate(tt$f, tt$space) # 1


## functionally dependent
## f) area of triangle
s = nint_space(nint_intvDim(0, 1),
               nint_funcDim(function(x) nint_intvDim(x[1]/2, 1 - x[1]/2)) )
s
nint_integrate(function(x) 1, s) # 0.5

## g) area of circle
s = nint_space(
    nint_intvDim(-1, 1),
    nint_funcDim(function(x) nint_intvDim( c(-1, 1) * sin(acos(x[1])) ))
)
s
nint_integrate(function(x) 1, s) # pi

## h) volume of sphere
s = nint_space(s[[1]],
               s[[2]],
               nint_funcDim(function(x) {
                   r = sin(acos(x[1]))
                   nint_intvDim(c(-1, 1) * r*cos(asin(x[2] / r)))
               }) )
s
nint_integrate(function(x) 1, s) # 4*pi/3

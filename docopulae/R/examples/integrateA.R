f = function(x) ifelse(x < 0, cos(x), sin(x))
#curve(f(x), -1, 1)
try(integrate(f, -1, 1, subdivisions=1)$value)
integrateA(f, -1, 1, subdivisions=1)$value
integrateA(f, -1, 1, subdivisions=2)$value
integrateA(f, -1, 1, subdivisions=3)$value

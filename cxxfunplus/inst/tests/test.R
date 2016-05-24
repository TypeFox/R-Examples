
require(cxxfunplus)
dso <- cxxfunctionplus(signature(x = "integer", y = "numeric"), 
                       'return ScalarReal(INTEGER(x)[0] * REAL(y)[0]);', 
                       save.dso = TRUE)
show(dso) 
fx <- grab.cxxfun(dso)
fx(3L, 4)

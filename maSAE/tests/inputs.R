#!/usr/bin/Rscript --vanilla
## from system.file("test-tools-1.R",      package = "Matrix"):
assertError <- function(expr) {
  d.expr <- deparse(substitute(expr))
  t.res <- tryCatch(expr, error = function(e) e)
  print(t.res)
  if(!inherits(t.res, "error"))
    stop(d.expr, "\n\t did not give an error", call. = FALSE)
  invisible(t.res)
}
assertWarning <- function(expr) {
    d.expr <- deparse(substitute(expr))
    t.res <- tryCatch(expr, warning = function(w)w)
    if(!inherits(t.res, "warning"))
        stop(d.expr, "\n\t did not give a warning", call. = FALSE)
    invisible(t.res)
}

library('maSAE')
library('methods')

## no errors
new("sadObj")
new("saeObj")
new("sadObj", f = y ~ NULL|g, data = data.frame(y=1, g='a'))
new("saeObj", f = y ~ x+z|g, data = data.frame(y=1, x = 1, z = 1, g='a', stp = FALSE, stp0 = FALSE), s2 = 'stp', s1 = 'stp0')
new("saeObj", f = y ~ x+z|g, data = data.frame(y=1, x = 1, z = 1, g='a', stp = TRUE, stp0 = TRUE), s2 = 'stp', s1 = 'stp0')
new("saeObj", f = y ~ x+z|g, data = data.frame(y=1, x = 1, z = 1, g='a', stp = FALSE), s2 = 'stp')
new("saeObj", f = y ~ x+z|g, data = data.frame(y=1, x = 1, z = 1, g='a', stp = TRUE) , s2= 'stp', smallAreaMeans = data.frame(g = 'a', x = 1))
new("saeObj", f = y ~ x+z|g, data = data.frame(y=1, x = 1, z = 1, g='a', cl = 1, stp = TRUE), cluster = 'cl' ,s2 = 'stp')
new("saeObj", f = y ~ x+z|g, data = data.frame(y=1, x = 1, z = 1, g='a', cl = 1, stp = TRUE, inc = TRUE), cluster = 'cl' , s2 = 'stp', include = 'inc')

## asserted errors
assertError(new("sadObj.virtual"))

assertError(new("sadObj", f = y ~ x:z))
assertError(new("sadObj", f = y ~ x*z))
assertError(new("sadObj", f = y ~ x+z))
assertError(new("sadObj", f = y ~ x+z|g))

assertError(new("sadObj", f = y ~ x+z|g, data = data.frame(y='a')))
assertError(new("sadObj", f = y ~ x+z|g, data = data.frame(y=1)))
assertError(new("sadObj", f = y ~ x+z|g, data = data.frame(y=1, g='a')))

assertError(new("saeObj", f = y ~ x:z))
assertError(new("saeObj", f = y ~ x*z))
assertError(new("saeObj", f = y ~ x+z))
assertError(new("saeObj", f = y ~ x+z|g))

assertError(new("saeObj", f = y ~ x+z|g, data = data.frame(y='a')))
assertError(new("saeObj", f = y ~ x+z|g, data = data.frame(y=1)))

assertError(new("saeObj", f = y ~ NULL|g, data = data.frame(y=1, g='a')))
assertError(new("saeObj", f = y ~ x+z|g, data = data.frame(y=1, g='a')))
assertError(new("saeObj", f = y ~ x+z|g, data = data.frame(y=1, x = 1, g='a')))
assertError(new("saeObj", f = y ~ x+z|g, data = data.frame(y=1, x = 1, z = 1, g='a')))
assertError(new("saeObj", f = y ~ x+z|g, data = data.frame(y=1, x = 1, z = 'a', g='a')))
assertError(new("saeObj", f = y ~ x+z|g, data = data.frame(y=1, x = 1, z = 1, g='a', cl = 1), cluster = 'cl'))
assertError(new("saeObj", f = y ~ x+z|g, data = data.frame(y=1, x = 1, z = 1, g='a') , s2 = 'stp'))
assertError(new("saeObj", f = y ~ x+z|g, data = data.frame(y=1, x = 1, z = 1, g='a', stp = 1)    , s2 = 'stp'))
assertError(new("saeObj", f = y ~ x+z|g, data = data.frame(y=1, x = 1, z = 1, g='a', stp = TRUE, stp0 = FALSE), s2 = 'stp', s1 = 'stp0'))
assertError(new("saeObj", f = y ~ x+z|g, data = data.frame(y=1, x = 1, z = 1, g='a', stp = TRUE), cluster = 'cl' ,s2 = 'stp'))

assertError(new("saeObj", f = y ~ x+z|g, data = data.frame(y=1, x = 1, z = 1, g='a') , smallAreaMeans = data.frame(g = 'b', x = 1, z = 1)))
assertError(new("saeObj", f = y ~ x+z|g, data = data.frame(y=1, x = 1, z = 1, g='a') , smallAreaMeans = data.frame(g = 'a', x = 1)))
assertError(new("saeObj", f = y ~ x+z|g, data = data.frame(y=1, x = 1, z = 1, g='a') , smallAreaMeans = data.frame(g = 'a', x = 1, z = 1, k = 1)))
assertError(new("saeObj", f = y ~ x+z|g, data = data.frame(y=1, x = 1, z = 1, g='a') , smallAreaMeans = data.frame(g = 'a', x = 1, z = NA)))
assertError(new("saeObj", f = y ~ x+z|g, data = data.frame(y=1:2, x = 1:2, z = 1:2, g=c('a', 'b')) , smallAreaMeans = data.frame(g = 'a', x = 1, z = 1)))
assertError(new("saeObj", f = y ~ x+z|g, data = data.frame(y=1, x = 1, z = 1, g='a', cl = 1, stp = TRUE), cluster = 'cl', s2 = 'stp', include = 'inc'))
assertError(new("saeObj", f = y ~ x+z|g, data = data.frame(y=1, x = 1, z = 1, g='a', cl = 1, stp = TRUE, inc = 1), cluster = 'cl', s2 = 'stp', include = 'inc'))
assertError(new("saeObj", f = y ~ x+z|g, data = data.frame(y=1, x = 1, z = 1, g='a', cl = 1, stp = TRUE, inc = TRUE), s2 = 'stp', include = 'inc'))

assertError(new("saeObj", f = y ~ x+z|g
		, data = data.frame(y=1, x = 1, z = 1, g='a', s1=TRUE)
		, smallAreaMeans = data.frame(g = 'a', x = 1, z = 1)
		, s1 = 's1')
)

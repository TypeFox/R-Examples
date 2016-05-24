## we use the example polygon and f (exponential decay) from
example(plotpolyf)

## numerical approximation of 'intrfr'
(intISOnum <- polyCub.iso(letterR, f, center=fcenter))

## analytical 'intrfr' (recall: f_r(r)=dexp(r), we need int_0^R r*f(r) dr)
intrfr <- function (R, rate=1) pgamma(R, 2, rate) / rate
(intISOana <- polyCub.iso(letterR, intrfr=intrfr, center=fcenter))

stopifnot(all.equal(intISOana, intISOnum, check.attributes=FALSE))


### polygon area: f(r) = 1, f(x,y) = 1, center does not really matter

intrfr.const <- function (R) R^2/2
(area.ISO <- polyCub.iso(letterR, intrfr=intrfr.const, center=c(0,0)))

stopifnot(all.equal(spatstat::area.owin(letterR),
                    area.ISO,
                    check.attributes=FALSE))
## the hole is subtracted correctly

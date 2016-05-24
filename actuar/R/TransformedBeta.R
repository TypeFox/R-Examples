### ===== actuar: An R Package for Actuarial Science =====
###
### Definition of the {d,p,q,r,m,lev}trbeta functions to compute
### characteristics of the Transformed Beta distribution. The version
### used in these functions has cumulative distribution function
###
###   Pr[X <= x] = Pr[Y <= (x/scale)^shape2 / (1 + (x/scale)^shape2)],
###
### where Y has a Beta distribution with parameters shape3 and shape1.
###
### See Appendix A of Klugman, Panjer & Willmot, Loss Models, Second
### Edition, Wiley, 2004.
###
### AUTHORS:  Mathieu Pigeon, Vincent Goulet <vincent.goulet@act.ulaval.ca>

dtrbeta <-
    function (x, shape1, shape2, shape3, rate = 1, scale = 1/rate,
              log = FALSE)
    .External("actuar_do_dpq", "dtrbeta", x, shape1, shape2, shape3, scale, log)

ptrbeta <-
    function (q, shape1, shape2, shape3, rate = 1, scale = 1/rate,
              lower.tail = TRUE, log.p = FALSE)
    .External("actuar_do_dpq", "ptrbeta", q, shape1, shape2, shape3, scale,
              lower.tail, log.p)

qtrbeta <-
    function (p, shape1, shape2, shape3, rate = 1, scale = 1/rate,
              lower.tail = TRUE, log.p = FALSE)
    .External("actuar_do_dpq", "qtrbeta", p, shape1, shape2, shape3, scale,
              lower.tail, log.p)

rtrbeta <-
    function (n, shape1, shape2, shape3, rate = 1, scale = 1/rate)
    .External("actuar_do_random", "rtrbeta", n, shape1, shape2, shape3, scale)

mtrbeta <-
    function (order, shape1, shape2, shape3, rate = 1, scale = 1/rate)
    .External("actuar_do_dpq", "mtrbeta", order, shape1, shape2, shape3, scale, FALSE)

levtrbeta <-
    function (limit, shape1, shape2, shape3, rate = 1, scale = 1/rate,
              order = 1)
    .External("actuar_do_dpq", "levtrbeta", limit, shape1, shape2, shape3, scale,
              order, FALSE)

## Aliases
dpearson6 <- dtrbeta
ppearson6 <- ptrbeta
qpearson6 <- qtrbeta
rpearson6 <- rtrbeta
mpearson6 <- mtrbeta
levpearson6 <- levtrbeta

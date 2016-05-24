### ===== actuar: An R Package for Actuarial Science =====
###
### Definition of the {d,p,q,r,m,lev}genbeta functions to compute
### characteristics of the Generalized Beta distribution. The version
### used in these functions has cumulative distribution function
###
###   Pr[X <= x] = Pr[Y <= (x/scale)^shape3], 0 < x < scale,
###
### where Y has a Beta distribution with parameters shape1 and shape2.
###
### See Appendix A of Klugman, Panjer & Willmot, Loss Models, Second
### Edition, Wiley, 2004.
###
### AUTHOR: Vincent Goulet <vincent.goulet@act.ulaval.ca>

dgenbeta <-
    function (x, shape1, shape2, shape3, rate = 1, scale = 1/rate,
              log = FALSE)
    .External("actuar_do_dpq", "dgenbeta", x, shape1, shape2, shape3, scale, log)

pgenbeta <-
    function (q, shape1, shape2, shape3, rate = 1, scale = 1/rate,
              lower.tail = TRUE, log.p = FALSE)
    .External("actuar_do_dpq", "pgenbeta", q, shape1, shape2, shape3, scale,
              lower.tail, log.p)

qgenbeta <-
    function (p, shape1, shape2, shape3, rate = 1, scale = 1/rate,
              lower.tail = TRUE, log.p = FALSE)
    .External("actuar_do_dpq", "qgenbeta", p, shape1, shape2, shape3, scale,
              lower.tail, log.p)

rgenbeta <-
    function (n, shape1, shape2, shape3, rate = 1, scale = 1/rate)
    .External("actuar_do_random", "rgenbeta", n, shape1, shape2, shape3, scale)

mgenbeta <-
    function (order, shape1, shape2, shape3, rate = 1, scale = 1/rate)
    .External("actuar_do_dpq", "mgenbeta", order, shape1, shape2, shape3, scale, FALSE)

levgenbeta <-
    function (limit, shape1, shape2, shape3, rate = 1, scale = 1/rate,
              order = 1)
    .External("actuar_do_dpq", "levgenbeta", limit, shape1, shape2, shape3, scale,
              order, FALSE)

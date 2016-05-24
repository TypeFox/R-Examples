pred.int.npar.k.of.m.on.r.conf.level <-
function (n, n.median, k, m, r, lpl.rank, upl.rank, pi.type, 
    integrate.args.list = NULL) 
{
    if (r == 1 && n.median == 1) {
        conf.level <- predIntNpar(x = 1:n, k = k, m = m, lpl.rank = lpl.rank, 
            n.plus.one.minus.upl.rank = n + 1 - upl.rank, pi.type = pi.type)$interval$conf.level
    }
    else {
        j <- ifelse(pi.type == "upper", upl.rank, n + 1 - lpl.rank)
        fcn.to.integrate <- function(x, n.median.weird, k.weird, 
            m.weird, r.weird, n.weird, j.weird) {
            fcn.to.integrate.vec <- function(x, n.median.weird, 
                k.weird, m.weird, r.weird, n.weird, j.weird) {
                if (n.median == 1) {
                  y <- x
                }
                else {
                  med.rank <- (n.median.weird + 1)/2
                  y <- pnbinom(q = n.median.weird - med.rank, 
                    size = med.rank, prob = x)
                }
                pnbinom(q = m.weird - k.weird, size = k.weird, 
                  prob = y)^r.weird * dbeta(x = x, shape1 = j.weird, 
                  shape2 = n.weird + 1 - j.weird)
            }
            sapply(x, fcn.to.integrate.vec, n.median.weird = n.median.weird, 
                k.weird = k.weird, m.weird = m.weird, r.weird = r.weird, 
                n.weird = n.weird, j.weird = j.weird)
        }
        if (!is.null(integrate.args.list)) {
            args.list <- c(list(f = fcn.to.integrate, lower = 0, 
                upper = 1, n.median.weird = n.median, k.weird = k, 
                m.weird = m, r.weird = r, n.weird = n, j.weird = j), 
                integrate.args.list)
            conf.level <- do.call("integrate", args.list)$value
        }
        else {
            conf.level <- integrate(fcn.to.integrate, lower = 0, 
                upper = 1, n.median.weird = n.median, k.weird = k, 
                m.weird = m, r.weird = r, n.weird = n, j.weird = j)$value
        }
    }
    conf.level
}

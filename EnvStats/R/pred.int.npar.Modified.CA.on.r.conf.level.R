pred.int.npar.Modified.CA.on.r.conf.level <-
function (n, n.median, r, lpl.rank, upl.rank, pi.type, integrate.args.list = NULL) 
{
    j <- ifelse(pi.type == "upper", upl.rank, n + 1 - lpl.rank)
    fcn.to.integrate <- function(x, n.median.weird, r.weird, 
        n.weird, j.weird) {
        fcn.to.integrate.vec <- function(x, n.median.weird, r.weird, 
            n.weird, j.weird) {
            if (n.median == 1) {
                y <- x
            }
            else {
                med.rank <- (n.median.weird + 1)/2
                y <- pnbinom(q = n.median.weird - med.rank, size = med.rank, 
                  prob = x)
            }
            q <- 1 - y
            y^r.weird * (1 + q + q^2 - 2 * q^3)^r.weird * dbeta(x = x, 
                shape1 = j.weird, shape2 = n.weird + 1 - j.weird)
        }
        sapply(x, fcn.to.integrate.vec, n.median.weird = n.median.weird, 
            r.weird = r.weird, n.weird = n.weird, j.weird = j.weird)
    }
    if (!is.null(integrate.args.list)) {
        args.list <- c(list(f = fcn.to.integrate, lower = 0, 
            upper = 1, n.median.weird = n.median, r.weird = r, 
            n.weird = n, j.weird = j), integrate.args.list)
        conf.level <- do.call("integrate", args.list)$value
    }
    else {
        conf.level <- integrate(fcn.to.integrate, lower = 0, 
            upper = 1, n.median.weird = n.median, r.weird = r, 
            n.weird = n, j.weird = j)$value
    }
    conf.level
}


clusOpt3fixedPSU <- function(unit.cost, m, delta1, delta2, unit.rv, k1=1, k2=1, CV0=NULL, tot.cost=NULL, cal.sw){
    if (any(delta1 < 0) | any(delta1 > 1)) stop("delta1 must be in [0,1].\n")
    if (any(delta2 < 0) | any(delta2 > 1)) stop("delta2 must be in [0,1].\n")
    if (!is.null(CV0) & !is.null(tot.cost))
        stop("CV0 and C.prime cannot both be non-null.\n")
    if (is.null(CV0) & is.null(tot.cost))
        stop("CV0 and C.prime cannot both be null.\n")

    if (sum(length(m)>1, length(delta1)>1,
        length(delta2)>1, length(unit.rv)>1, length(CV0)>1, length(tot.cost)>1) > 1)
            stop("Only one argument to function can be vector.\n")

    C1 <- unit.cost[1]
    C2 <- unit.cost[2]
    C3 <- unit.cost[3]
    C.prime <- tot.cost/m - C1

    q.opt <- sqrt((1-delta2)/delta2 * C2 / C3)

    if (cal.sw == 1){
        n <- C.prime/(C2 + C3*q.opt)
        if (n < 0) stop(paste("n is negative. Check inputs. n=",n,"\n"))

        tot.cost <- C1*m + C2*m*n + C3*m*n*q.opt
        CV <- sqrt(unit.rv/m/n/q.opt * (k1*delta1*n*q.opt + k2*(1 + delta2*(q.opt-1))))

        output <-
           structure(list(C1 = C1,
                          C2 = C2,
                          C3 = C3,
                          m = m,
                          delta1 = delta1,
                          delta2 = delta2,
                          "unit relvar" = unit.rv,
                          k1 = k1,
                          k2 = k2,
                          "variable budget" = tot.cost-C1*m,
                          "total cost" = round(tot.cost,0),
                          n = round(n,1),
                          q = round(q.opt,1),
                          CV = round(CV,4)),
                     class = "power.htest")
    }
    if (cal.sw == 2) {
        n <- k2*(1 + delta2*(q.opt-1)) / q.opt / (CV0^2*m/unit.rv - k1*delta1)
        if (n < 0) stop(paste("n is negative. Check inputs. n=",n,"\n"))

        tot.cost <- C1*m + C2*m*n + C3*m*n*q.opt
        CV.chk <- sqrt(unit.rv/m/n/q.opt * (k1*delta1*n*q.opt + k2*(1+ delta2*(q.opt-1))))

        output <-
           structure(list(C1 = C1,
                          C2 = C2,
                          C3 = C3,
                          m = m,
                          delta1 = delta1,
                          delta2 = delta2,
                          "unit relvar" = unit.rv,
                          k1 = k1,
                          k2 = k2,
                          "variable budget" = tot.cost-C1*m,
                          "total cost" = round(tot.cost,0),
                          n = round(n,1),
                          q = round(q.opt,1),
                          CV = CV0,
                          "CV check" = round(CV.chk,4)),
                     class = "power.htest")
    }
    output
}

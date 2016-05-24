
clusOpt3 <- function(unit.cost, delta1, delta2, unit.rv, k1=1, k2=1, CV0=NULL, tot.cost=NULL, cal.sw){
    if (length(unit.cost) != 3) stop("cost must be a vector with 3 components")
    if (!is.null(CV0) & !is.null(tot.cost))
        stop("CV0 and tot.cost cannot both be non-null.\n")
    if (sum(length(delta1)>1, length(delta2)>1, length(unit.rv)>1,
            length(CV0)>1, length(tot.cost)>1) > 1)
            stop("Only one argument to function can be vector.\n")

    if (cal.sw==1 & is.null(tot.cost))
            stop("If cal.sw=1, tot.cost must be non-null.\n")
    if (cal.sw==2 & is.null(CV0))
            stop("If cal.sw=2, CV0 must be non-null.\n")

    C1 <- unit.cost[1]
    C2 <- unit.cost[2]
    C3 <- unit.cost[3]

    q.opt <- sqrt( (1-delta2)/delta2 * C2 / C3 )
    n.opt <- 1/q.opt * sqrt((1-delta2)/delta1 * C1 / C3 * k2 / k1)

    if (cal.sw == 1){
        m.opt <- tot.cost / (C1 + C2*n.opt + C3*n.opt*q.opt)
        CV <- sqrt(unit.rv/m.opt/n.opt/q.opt * (delta1*n.opt*q.opt + 1 + delta2*(q.opt-1)))
        cost.chk <- C1*m.opt + C2*m.opt*n.opt + C3*m.opt*n.opt*q.opt

        output <-
           structure(list(C1 = C1,
                      C2 = C2,
                      C3 = C3,
                      delta1 = delta1,
                      delta2 = delta2,
                      "unit relvar" = unit.rv,
                      k1 = k1,
                      k2 = k2,
                      cost = tot.cost,
                      m.opt = round(m.opt,1),
                      n.opt = round(n.opt,1),
                      q.opt = round(q.opt,1),
                      CV = round(CV,4)),
                 class = "power.htest")
    }
    if (cal.sw == 2) {
        m.opt <- unit.rv/CV0^2/n.opt/q.opt * (k1*delta1*n.opt*q.opt + k2*(1 + delta2*(q.opt-1)))
        CV.chk <- sqrt(unit.rv/m.opt/n.opt/q.opt * (delta1*n.opt*q.opt + 1 + delta2*(q.opt-1)))
        cost.chk <- C1*m.opt + C2*m.opt*n.opt + C3*m.opt*n.opt*q.opt

        output <-
           structure(list(C1 = C1,
                      C2 = C2,
                      C3 = C3,
                      delta1 = delta1,
                      delta2 = delta2,
                      "unit relvar" = unit.rv,
                      k1 = k1,
                      k2 = k2,
                      cost  = cost.chk,
                      m.opt = round(m.opt,1),
                      n.opt = round(n.opt,1),
                      q.opt = round(q.opt,1),
                      CV0 = CV0,
                      CV.chk = round(CV.chk,4)),
                 class = "power.htest")
    }
    output
}

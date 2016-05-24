
clusOpt2 <- function(C1, C2, delta, unit.rv, k=1, CV0=NULL, tot.cost=NULL, cal.sw){
    options(warn = -1)
    if (!is.null(CV0) & !is.null(tot.cost))
        stop("CV0 and tot.cost cannot both be non-null.\n")
    if (sum(length(C1)>1, length(C2)>1, length(delta)>1,
        length(unit.rv)>1, length(CV0)>1, length(tot.cost)>1) > 1)
            stop("Only one argument to function can be vector.\n")

    c.ratio <- C1/C2
    n.opt <- sqrt(c.ratio * (1-delta)/delta)

    if (cal.sw == 1){
        m.opt <- tot.cost / (C1 + C2*n.opt)
        if (m.opt < 0) stop(paste("m.opt is negative. Check inputs. m.opt=",m.opt,"\n"))

        CV <- sqrt(unit.rv/m.opt/n.opt*(1 + delta*(n.opt-1)))
        output <-
           structure(list(C1 = C1,
                          C2 = C2,
                          delta = delta,
                          "unit relvar" = unit.rv,
                          k = k,
                          cost = tot.cost,
                          m.opt = round(m.opt,1),
                          n.opt = round(n.opt,1),
                          CV = round(CV,4)),
                     class = "power.htest")
    }
    if (cal.sw == 2) {
        m.opt <- unit.rv * k * (1 + delta*(n.opt-1)) / n.opt / CV0^2
        if (m.opt < 0) stop(paste("m.opt is negative. Check inputs. m.opt=",m.opt,"\n"))

        cost <- C1*m.opt + C2*m.opt*n.opt
        output <-
           structure(list(C1 = C1,
                          C2 = C2,
                          delta = delta,
                          "unit relvar" = unit.rv,
                          k = k,
                          cost = cost,
                          m.opt = round(m.opt,1),
                          n.opt = round(n.opt,1),
                          CV = CV0),
                     class = "power.htest")
    }
    output
}

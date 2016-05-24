
clusOpt2fixedPSU <- function(C1, C2, m, delta, unit.rv, k=1, CV0=NULL, tot.cost=NULL, cal.sw){
    if (any(delta < 0) | any(delta > 1)) stop("delta must be in [0,1].\n")
    if (!is.null(CV0) & !is.null(tot.cost))
        stop("CV0 and tot.cost cannot both be non-null.\n")
    if (is.null(CV0) & is.null(tot.cost))
        stop("CV0 and tot.cost cannot both be null.\n")

    if (sum(length(C1)>1, length(C2)>1, length(m)>1, length(delta)>1,
        length(unit.rv)>1, length(CV0)>1, length(tot.cost)>1) > 1)
            stop("Only one argument to function can be vector.\n")

    if (cal.sw == 1){
        n <- (tot.cost - C1*m)/C2/m
        if (n < 0) stop(paste("n is negative. Check inputs. n=",n,"\n"))
        CV <- sqrt(unit.rv*k/m/n*(1 + delta*(n-1)))
        output <-
           structure(list(C1 = C1,
                          C2 = C2,
                          m = m,
                          delta = delta,
                          "unit relvar" = unit.rv,
                          k = k,
                          cost = tot.cost,
                          n = round(n,1),
                          CV = round(CV,4)),
                     class = "power.htest")
    }
    if (cal.sw == 2) {
        n <- (1 - delta) / (CV0^2*m/unit.rv/k - delta)
        if (n < 0) stop(paste("n is negative. Check inputs. n=",n,"\n"))
        cost <- C1*m + C2*m*n
        output <-
           structure(list(C1 = C1,
                          C2 = C2,
                          delta = delta,
                          m = m,
                          "unit relvar" = unit.rv,
                          k = k,
                          cost = cost,
                          n = round(n,1),
                          CV = CV0),
                     class = "power.htest")
    }
    output
}

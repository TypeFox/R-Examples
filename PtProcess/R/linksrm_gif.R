linksrm_gif <-
function (data, evalpts, params, TT = NA, tplus = FALSE, eta = 0.75) 
{
    if (tplus) 
        FUN <- "<="
    else FUN <- "<"
    n <- sqrt(length(params) + 1) - 1
    a <- params[1:n]
    b <- params[(n + 1):(2 * n)]
    cc <- matrix(params[(2 * n + 1):length(params)], nrow = n, 
        byrow = TRUE)
#   if (!any(search()=="PtProcess.tmp"))
#       attach(NULL, pos=2L, name="PtProcess.tmp", warn.conflicts=TRUE)
    if (FALSE) print("attach new database to search path")
    if (any(is.na(TT))) {
        evaltimes <- evalpts[, "time"]
        evalregion <- evalpts[, "region"]
        if (!is.null(data)) {
            times <- data[, "time"]
            magnitude <- data[, "magnitude"]
            region <- data[, "region"]
#           if (!exists("St1", mode = "numeric")) {
            if (FALSE) print("this loop calculates St1")
                St1 <- matrix(rep(NA, (nrow(evalpts) * n)), 
                  nrow = n)
                for (i in 1:n) {
                  if (length(magnitude[region == i]) == 0) 
                      St1[i, ] <- rep(0, length(evaltimes))
                  else St1[i, ] <- matrix(10^(eta * magnitude[region == 
                      i]), nrow = 1) %*% outer(times[region == 
                      i], evaltimes, FUN = FUN)
                }
#               assign("St1", St1, pos="PtProcess.tmp")
                if (FALSE) print("assign statement for St1")
#           }
            if (FALSE) print("end loop St1")
            ci <- exp(a[evalregion] + b[evalregion] * (evaltimes - 
                matrix(rep(1, n), nrow = 1) %*% (t(cc[evalregion, 
                  ]) * St1)))
        }
        else ci <- exp(a[evalregion] + b[evalregion] * evaltimes)
    }
    else {
        if (!is.null(data)) {
            times <- data[, "time"]
            magnitude <- data[, "magnitude"]
            region <- data[, "region"]
            within <- (times < TT[2]) & (times >= TT[1])
            ti <- c(times[within], TT[2])
            nt <- length(ti)
#           if (!exists("St2", mode = "numeric")) {
            if (FALSE) print("this loop calculates St2")
                St2 <- matrix(rep(NA, (nt * n)), nrow = n)
                for (i in 1:n) {
                    if (length(magnitude[region == i]) == 0) 
                        St2[i, ] <- rep(0, length(ti))
                    else St2[i, ] <- matrix(10^(eta * magnitude[region == 
                        i]), nrow = 1) %*% outer(times[region == i], 
                        ti, FUN = FUN)
                }
#               assign("St2", St2, pos="PtProcess.tmp")
                if (FALSE) print("assign statement for St2")
#           }
            if (FALSE) print("end loop St2")
            ci <- (diag(1/b, nrow = n) %*% exp(matrix(rep(a, nt),
                nrow = n) - diag(b, nrow = n) %*% cc %*% St2) *
                (exp(matrix(b, ncol = 1) %*% matrix(ti, nrow = 1)) -
                exp(matrix(b, ncol = 1) %*% matrix(c(TT[1], ti[-nt]),
                nrow = 1)))) %*% matrix(1, nrow = nt, ncol = 1)
        }
        else ci <- (diag(1/b, nrow = n) %*% exp(matrix(a, ncol = 1)) *
            (exp(matrix(b * TT[2], ncol = 1)) - exp(matrix(b *
            TT[1])))) %*% matrix(1, nrow = nt, ncol = 1)
    }
    return(as.vector(ci))
}
attr(linksrm_gif, "rate") <- "increasing"
attr(linksrm_gif, "regions") <- expression(sqrt(length(params) + 1) - 1)


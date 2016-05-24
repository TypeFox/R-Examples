psi.hat <- function(x, timehoriz, model, diag = FALSE, ...) {
    aux   <- alpha(model)(timehoriz / x)
    model <- tilt(model, -aux[1L])

    Y <- rpsim(model)(timehoriz = timehoriz, ...)

    ind <- apply(X      = Y,
                 MARGIN = 2L,
                 FUN    = function(arg) {
                     tryCatch(expr    = min(which(arg + x <= 0.0)),
                              warning = function(w) NA_integer_,
                              error   = function(e) NA_integer_)
                 })

    LT <- rbind(-Y[cbind(ind, seq_len(ncol(Y)))],
                time(Y)[ind])

    val             <- exp(drop(crossprod(aux, LT)))
    val[is.na(val)] <- 0.0

    psi <- mean(val * vapply(LT[2L, ] <= timehoriz,
                             isTRUE,
                             logical(1L)))

    if (diag) {
        return(list(psi = psi,
                    LT  = LT,
                    Y   = Y,
                    par = c(x = x,
                            t = timehoriz),
                    val = val,
                    ind = ind,
                    aux = aux))
    } else {
        return(list(psi = psi,
                    LT  = LT))
    }
}

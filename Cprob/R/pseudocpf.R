pseudocpf <- function(formula, data, id, subset, na.action, timepoints, 
                      failcode = 1, ...) {
                      
    call <- match.call()
    m <- match.call(expand.dots = FALSE)
    temp <- c("", "formula", "data", "id", "subset", "na.action")
    m <- m[match(temp, names(m), nomatch = 0)]
    Terms <- if (missing(data)) terms(formula)
    else terms(formula, data = data)
    m$formula <- Terms
    m[[1]] <- as.name("model.frame")
    m <- eval(m, parent.frame())

    response <- model.extract(m, "response")
    if (!inherits(response, "Hist")) stop("Response must be a 'Hist' object")
    if (attr(response, "model") != "competing.risks")
        stop("This is not competing risk data")
    if (attr(response, "cens.type") != "rightCensored")
        stop("Works only for right-censored data")
    
    event <- response[, "event"]
    time <- response[, "time"]
    states <- attr(response, "states")
    id <- model.extract(m, "id")
    event <- factor(event, ordered = TRUE)
    cc <- attr(response, "cens.code")
    levels(event) <- c(states, cc)
    daten <- data.frame(id, event, time)
    tmax <- max(daten$time) + 10^-3
    n <- nrow(daten)
    nt <- length(timepoints)
    psd <- matrix(0, nrow = n * nt, ncol = 3)
    ref <- matrix(predict(cpf(Hist(time, event, cens.code = cc) ~ 1, daten,
                              failcode = failcode), timepoints)$cp,
                  ncol = 1, nrow = nt)
    ref <- apply(ref, 2, rep, n)
    temp <- lapply(seq_len(n), function(i) {
        matrix(predict(cpf(Hist(time, event, cens.code = cc) ~ 1, daten[-i, ],
                           failcode = failcode), timepoints)$cp,
               ncol = 1, nrow = nt)
    })
    temp <- do.call(rbind, temp)
    
    psd[, 2] <- n * ref - (n - 1) * temp
    psd[, 1] <- as.vector(mapply(rep, id, nt))
    psd[, 3] <- rep(timepoints, n)
    psd <- as.data.frame(psd)
    names(psd) <- c("id", "pseudo", "time")
    cov <- model.matrix(Terms, m)[, -1, drop = FALSE]
    dat <- cbind(psd, matrix(mapply(rep, cov, nt), dim(psd)[1], dim(cov)[2]))
    names(dat)[4:dim(dat)[2]] <- colnames(cov)
    
    newf <- update(formula, pseudo ~ factor(time) + .)
    fit <- geepack::geese(newf, id = id, data = dat, family = gaussian,
                          mean.link = "logit", ...)
    
    zzz <- list(fit = fit, pseudo = psd, timepoints = timepoints, call = call)
    class(zzz) <- "pseudocpf"
    zzz
}

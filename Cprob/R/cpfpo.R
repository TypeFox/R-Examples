cpfpo <- function(formula, data, subset, na.action, failcode,
                  tis, w, ...) {
    
    m <- match.call(expand.dots = FALSE)
    temp <- c("", "formula", "data", "subset", "na.action")
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
    
    if (missing(failcode)) failcode <- states[1]
    
    event <- factor(event, ordered = TRUE)
    ## levels(event) <- c(attr(response, "cens.code"), states)
    levels(event) <- c(states, attr(response, "cens.code"))
    daten <- data.frame(event, time)
    
    cov <- model.matrix(Terms, m)
    
    tmax <- max(daten$time) + 10^-3
    ev.comp <- levels(event)[levels(event) != failcode &
                             levels(event) != attr(response, "cens.code")]
    
    y <- lapply(seq_len(nrow(daten)), function(i) {
        if (daten$event[i] == failcode) {
            time <- c(0, daten$time[i], tmax)
            cov <- c(0, 1, 1)
        }
        else {
            time <- c(0, tmax)
            cov <- c(0, 0)
        }
        return(as.lgtdl(data.frame(time, cov)))
    })
    
    s <- lapply(seq_len(nrow(daten)), function(i) {
        if (daten$event[i] %in% ev.comp) {
            time <- c(0, daten$time[i], tmax)
            cov <- c(1, 0, 0)
        }
        else {
            time <- c(0, tmax)
            cov <- c(1, 1)
        }
        return(as.lgtdl(data.frame(time, cov)))
    })
    
    if (missing(w)) {
        w <- rep(1, length(tis))
    }
    
    fit <- tpr(y, s, cov[, 1:ncol(cov)], list(),
               cov[, -(1:ncol(cov))], list(),
               w = w, tis = tis, family = binomial(),
               evstr = list(link = 2, v = 2), ...)
    
    class(fit) <- c("tpr", "cpfpo")
    
    fit
}

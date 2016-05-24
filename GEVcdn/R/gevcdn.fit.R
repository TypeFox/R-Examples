gevcdn.fit <-
function (x, y, iter.max = 1000, n.hidden = 2, Th = gevcdn.logistic,
          fixed = NULL, init.range = c(-0.25, 0.25),
          scale.min = .Machine$double.eps, beta.p = 3.3, beta.q = 2,
          sd.norm = Inf, n.trials = 5,
          method = c("BFGS", "Nelder-Mead"), max.fails = 100, ...)
{
    if (!is.matrix(x)) stop("\"x\" must be a matrix")
    if (!is.matrix(y)) stop("\"y\" must be a matrix")
    method <- match.arg(method)
    if (identical(Th, gevcdn.identity)) n.hidden <- 3
    GML <- Inf
    x.min <- apply(x, 2, min)
    x.max <- apply(x, 2, max)
    x <- sweep(sweep(x, 2, x.min, '-'), 2, x.max - x.min, '/')
    y.min <- min(y)
    y.max <- max(y)
    y <- (y - y.min)/(y.max - y.min)
    for (i in seq_len(n.trials)){
        cat(i,"")
        exception <- TRUE
        exception.count <- 0
        while (exception){
            if (identical(names(init.range), c("W1", "W2"))){
                weights <- unlist(init.range) + gevcdn.initialize(x,
                               n.hidden, c(-0.25, 0.25))
            } else{
                weights <- gevcdn.initialize(x, n.hidden, init.range)
            }
            fit.cur <- try(suppressWarnings(optim(weights, gevcdn.cost,
                               method = method,
                               control = list(maxit = iter.max, ...),
                               x = x, y = y, n.hidden = n.hidden,
                               Th = Th, fixed = fixed,
                               scale.min = scale.min,
                               beta.p = beta.p, beta.q = beta.q,
                               sd.norm = sd.norm)),
                           silent = TRUE)
            if (!class(fit.cur) == "try-error"){
                exception <- fit.cur$value > 1e+308
            }
            if (exception) exception.count <- exception.count + 1
            if (exception.count == max.fails){
                stop("exception... check arguments")
            }
        }
        GML.cur <- fit.cur$value
        cat(GML.cur,"")
        if (GML.cur < GML){
            GML <- GML.cur
            output.cdn <- fit.cur
        }
    }
    cat("\n")
    weights <- output.cdn$par
    cost <- gevcdn.cost(weights, x, y, n.hidden, Th, fixed,
                        scale.min, beta.p, beta.q, sd.norm)
    w <- gevcdn.reshape(x, weights, n.hidden)
    attr(w, "x.min") <- x.min
    attr(w, "x.max") <- x.max
    attr(w, "y.min") <- y.min
    attr(w, "y.max") <- y.max
    attr(w, "Th") <- Th
    attr(w, "fixed") <- fixed
    attr(w, "scale.min") <- scale.min
    NLL <- attr(cost, "NLL")
    penalty <- attr(cost, "penalty")
    attr(w, "GML") <- GML
    attr(w, "NLL") <- NLL
    attr(w, "penalty") <- penalty
    if (sd.norm == Inf){
        if (length(fixed)==3){
            k <- 3
        } else{
            if (identical(Th, gevcdn.identity)){
                k <- (3-length(fixed)) * (ncol(x) + 1) + length(fixed)
            } else{
                k <- length(weights) - n.hidden * length(fixed)
            }
        }
        n <- nrow(y)
        BIC <- 2 * NLL + k * log(n)
        AIC <- 2 * NLL + 2 * k
        AICc <- AIC + (2 * k * (k + 1))/(n - k - 1)
        attr(w, "BIC") <- BIC
        attr(w, "AIC") <- AIC
        attr(w, "AICc") <- AICc
        attr(w, "k") <- k
    }
    w
}


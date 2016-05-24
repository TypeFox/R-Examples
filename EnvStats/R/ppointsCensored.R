ppointsCensored <-
function (x, censored, censoring.side = "left", prob.method = "michael-schucany", 
    plot.pos.con = 0.375) 
{
    if (!is.vector(x, mode = "numeric")) 
        stop("'x' must be a numeric vector")
    if (!is.vector(censored, mode = "numeric") & !is.vector(censored, 
        mode = "logical")) 
        stop("'censored' must be a logical or numeric vector")
    if (length(censored) != length(x)) 
        stop("'censored' must be the same length as 'x'")
    if ((bad.obs <- sum(!(ok <- is.finite(x) & is.finite(as.numeric(censored))))) > 
        0) {
        is.not.finite.warning(x)
        is.not.finite.warning(as.numeric(censored))
        x <- x[ok]
        censored <- censored[ok]
        warning(paste(bad.obs, "observations with NA/NaN/Inf in 'x' and 'censored' removed."))
    }
    if (is.numeric(censored)) {
        if (!all(censored == 0 | censored == 1)) 
            stop(paste("When 'censored' is a numeric vector, all values of", 
                "'censored' must be 0 (not censored) or 1 (censored)."))
        censored <- as.logical(censored)
    }
    N <- length(x)
    n <- sum(!censored)
    n.cen <- N - n
    if (n.cen == 0) 
        stop(paste("No censored values indicated by 'censored'. ", 
            "Use the function 'ppoints'."))
    x.no.cen <- x[!censored]
    if (length(unique(x.no.cen)) < 1) 
        stop("'x' must contain at least one non-missing, uncensored value.")
    prob.method <- match.arg(prob.method, c("michael-schucany", 
        "hirsch-stedinger", "kaplan-meier", "modified kaplan-meier", 
        "nelson"))
    censoring.side <- match.arg(censoring.side, c("left", "right"))
    if (censoring.side == "left" && prob.method == "nelson") 
        stop("Nelson method not available for censoring.side='left'")
    if (censoring.side == "right" && prob.method == "modified kaplan-meier") 
        stop("Modified Kaplan-Meier method not available for censoring.side='right'")
    if (!is.vector(plot.pos.con, mode = "numeric") || length(plot.pos.con) != 
        1 || plot.pos.con < 0 || plot.pos.con > 1) 
        stop("'plot.pos.con' must be a numeric scalar between 0 and 1")
    x.cen <- x[censored]
    c.vec <- table(x.cen)
    cen.levels <- sort(unique(x.cen))
    K <- length(cen.levels)
    ord.stats <- sort(x)
    p <- numeric(N)
    temp <- x
    diffs <- diff(ord.stats)
    eps <- min(diffs[diffs > 0])/2
    if (censoring.side == "left") {
        temp[censored] <- x[censored] - eps
        ord.stats.cen <- censored[order(temp)]
        switch(prob.method, `michael-schucany` = {
            j <- rev((1:N)[!ord.stats.cen])
            cond.probs <- (j - plot.pos.con)/(j - plot.pos.con + 
                1)
            p.no.cen <- rev(((N - plot.pos.con + 1)/(N - 2 * 
                plot.pos.con + 1)) * cumprod(cond.probs))
        }, `hirsch-stedinger` = {
            new.cen.levels <- c(-Inf, cen.levels, Inf)
            A <- numeric(K + 1)
            B <- A
            for (j in 0:K) {
                B[j + 1] <- sum(temp < new.cen.levels[j + 1])
                A[j + 1] <- sum(new.cen.levels[j + 1] <= x.no.cen & 
                  x.no.cen < new.cen.levels[j + 2])
            }
            S <- numeric(K + 2)
            S[1] <- 1
            S[K + 2] <- 0
            for (j in K:1) {
                S[j + 1] <- S[j + 2] + (A[j + 1]/(A[j + 1] + 
                  B[j + 1])) * (1 - S[j + 2])
            }
            p.no.cen <- numeric(n)
            i <- 1
            for (j in 0:K) {
                if ((m <- A[j + 1]) > 0) {
                  p.no.cen[(i - 1) + (1:m)] <- (1 - S[j + 1]) + 
                    (S[j + 1] - S[j + 2]) * ppoints(n = m, a = plot.pos.con)
                  i <- i + m
                }
            }
        }, `modified kaplan-meier` = {
            ord.stats.no.cen <- ord.stats[!ord.stats.cen]
            rle.list <- rle(ord.stats)
            y.j <- rle.list$values
            m.j <- rle.list$lengths
            n.j <- sapply(y.j, function(z) sum(z >= ord.stats))
            d.j <- sapply(y.j, function(z) sum(z == ord.stats.no.cen))
            p <- c(rev(cumprod(rev((n.j - d.j)/n.j)))[-1], (N - 
                0.375)/(N + 0.25))
            p <- rep(p, times = m.j)
        }, `kaplan-meier` = {
            ord.stats.no.cen <- ord.stats[!ord.stats.cen]
            rle.list <- rle(ord.stats)
            y.j <- rle.list$values
            m.j <- rle.list$lengths
            n.j <- sapply(y.j, function(z) sum(z >= ord.stats))
            d.j <- sapply(y.j, function(z) sum(z == ord.stats.no.cen))
            p <- c(rev(cumprod(rev((n.j - d.j)/n.j)))[-1], 1)
            p <- rep(p, times = m.j)
        })
    }
    else {
        temp[censored] <- x[censored] + eps
        ord.stats.cen <- censored[order(temp)]
        j <- (1:N)[!ord.stats.cen]
        switch(prob.method, `michael-schucany` = {
            cond.probs <- (N - j - plot.pos.con + 1)/(N - j - 
                plot.pos.con + 2)
            p.no.cen <- 1 - ((N - plot.pos.con + 1)/(N - 2 * 
                plot.pos.con + 1)) * cumprod(cond.probs)
        }, `kaplan-meier` = {
            cond.probs <- (N - j)/(N - j + 1)
            p.no.cen <- 1 - cumprod(cond.probs)
            ord.stats.not.censored <- ord.stats[!ord.stats.cen]
            rle.list <- rle(ord.stats.not.censored)
            lengths <- rle.list$lengths
            values <- rle.list$values
            index <- (1:length(lengths))[lengths > 1]
            if (n.ties <- length(index)) {
                for (i in 1:n.ties) {
                  index2 <- ord.stats.not.censored == values[index[i]]
                  p.no.cen[index2] <- max(p.no.cen[index2])
                }
            }
        }, nelson = {
            cond.probs <- exp(-1/(N - j + 1))
            p.no.cen <- 1 - cumprod(cond.probs)
        }, `hirsch-stedinger` = {
            new.cen.levels <- c(-Inf, cen.levels, Inf)
            A <- numeric(K + 1)
            B <- A
            for (j in 0:K) {
                B[j + 1] <- sum(temp > new.cen.levels[j + 2])
                A[j + 1] <- sum(new.cen.levels[j + 1] < x.no.cen & 
                  x.no.cen <= new.cen.levels[j + 2])
            }
            S <- numeric(K + 2)
            S[1] <- 1
            S[K + 2] <- 0
            for (j in 1:K) {
                S[j + 1] <- S[j] * (1 - A[j]/(A[j] + B[j]))
            }
            p.no.cen <- numeric(n)
            i <- 1
            for (j in 0:K) {
                if ((m <- A[j + 1]) > 0) {
                  p.no.cen[(i - 1) + (1:m)] <- (1 - S[j + 1]) + 
                    (S[j + 1] - S[j + 2]) * ppoints(n = m, a = plot.pos.con)
                  i <- i + m
                }
            }
        })
    }
    if (!(censoring.side == "left" & prob.method %in% c("kaplan-meier", 
        "modified kaplan-meier"))) {
        p[!ord.stats.cen] <- p.no.cen
        if (prob.method == "hirsch-stedinger") {
            if (censoring.side == "right") {
                for (j in 1:K) {
                  index <- ord.stats.cen & (ord.stats == cen.levels[j])
                  p[index] <- rev(1 - S[j + 1] * ppoints(n = c.vec[j], 
                    a = plot.pos.con))
                }
            }
            else {
                for (j in 1:K) {
                  index <- ord.stats.cen & (ord.stats == cen.levels[j])
                  p[index] <- (1 - S[j + 1]) * ppoints(n = c.vec[j], 
                    a = plot.pos.con)
                }
            }
        }
        else {
            ord.stats.no.cen <- ord.stats[!ord.stats.cen]
            if (censoring.side == "right") {
                for (j in 1:K) {
                  Tj <- cen.levels[j]
                  index <- ord.stats.cen & (ord.stats == Tj)
                  index1 <- ord.stats.no.cen <= Tj
                  if (!any(index1)) 
                    p[index] <- 0
                  else p[index] <- p.no.cen[max((1:n)[index1])]
                }
            }
            else {
                for (j in 1:K) {
                  Tj <- cen.levels[j]
                  index <- ord.stats.cen & (ord.stats == Tj)
                  index1 <- ord.stats.no.cen >= Tj
                  if (!any(index1)) 
                    p[index] <- 1
                  else p[index] <- p.no.cen[min((1:n)[index1])]
                }
            }
        }
    }
    ret.list <- list(Order.Statistics = ord.stats, Cumulative.Probabilities = p, 
        Censored = ord.stats.cen, Censoring.Side = censoring.side, 
        Prob.Method = prob.method)
    if (prob.method %in% c("michael-schucany", "hirsch-stedinger")) 
        ret.list <- c(ret.list, list(Plot.Pos.Con = plot.pos.con))
    ret.list
}

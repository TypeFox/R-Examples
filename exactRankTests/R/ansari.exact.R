ansari.exact <- function(x, ...) UseMethod("ansari.exact")

ansari.exact.default <-
function(x, y, alternative = c("two.sided", "less", "greater"),
         exact = NULL, conf.int = FALSE, conf.level = 0.95, ...) 
{
    alternative <- match.arg(alternative)
    if(conf.int) {
        if(!((length(conf.level) == 1)
             && is.finite(conf.level)
             && (conf.level > 0)
             && (conf.level < 1)))
            stop("conf.level must be a single number between 0 and 1")
    }
    DNAME <- paste(deparse(substitute(x)), "and", deparse(substitute(y)))

    if(!is.numeric(x)) stop("`x' must be numeric")
    if(!is.numeric(y)) stop("`y' must be numeric")

    x <- x[complete.cases(x)]
    y <- y[complete.cases(y)]
    m <- length(x)
    if(m < 1)
        stop("not enough x observations")
    n <- length(y)
    if(n < 1)
        stop("not enough y observations")
    N <- m + n

    r <- cscores(c(x, y), type="Ansari")
    STATISTIC <- sum(r[seq(along = x)])

    if(is.null(exact))
        exact <- ((m < 50) && (n < 50))

    if(exact) {
        PVAL <-
            switch(alternative,
                   two.sided = pperm(STATISTIC, r, m, alternative="two.sided"),
                   less = pperm(STATISTIC, r, m, alternative="less"),
                   greater = pperm(STATISTIC, r, m, alternative="greater"))
        if (conf.int) {
            alpha <- 1 - conf.level
            x <- sort(x)
            y <- sort(y)
            ab <- function(sig, zq = 0) {
                rab <- rank(c(x/sig, y))
                sum(pmin(rab, N - rab + 1)[seq(along = x)])
            }
            ratio <- outer(x,y,"/")
            aratio <- ratio[ratio >= 0]
            sigma <- sort(aratio)

            cci <- function(alpha) {
              u <- absigma - qperm(alpha/2,  r, m) 
              l <- absigma - qperm(1 - alpha/2, r, m) 
              ## Check if the statistic exceeds both quantiles first.
              uci <- NULL
              lci <- NULL                    
              if(length(u[u >= 0]) == 0 || length(l[l > 0]) == 0) {
                  warning(paste("Samples differ in location: Cannot",
                                "compute confidence set, returning NA"))
                  return(c(NA, NA))
              } 
              if (is.null(uci)) {
                  u[u < 0] <- NA
                  uci <- min(sigma[which(u == min(u, na.rm=TRUE))])
              }
              if (is.null(lci)) {
                  l[l <= 0] <- NA
                  lci <- max(sigma[which(l == min(l, na.rm=TRUE))])
              }
              ## The process of the statistics does not need to be
              ## monotone in sigma: check this and interchange quantiles.
              if (uci > lci) {
                  l <- absigma - qperm(alpha/2,  r, m)
                  u <- absigma - qperm(1 - alpha/2, r, m)
                  u[u < 0] <- NA
                  uci <- min(sigma[which(u == min(u, na.rm=TRUE))])
                  l[l <= 0] <- NA
                  lci <- max(sigma[which(l == min(l, na.rm=TRUE))])
               }
               c(uci, lci)
            }

            cint <- if(length(sigma) < 1) {
                warning("Cannot compute confidence set, returning NA")
                c(NA, NA)
            }
            else {
                ## Compute statistics directly: computing the steps is
                ## not faster.
                absigma <-
                    sapply(sigma + c(diff(sigma)/2,
                                     sigma[length(sigma)]*1.01), ab)
                switch(alternative, two.sided = {
                    cci(alpha)
                }, greater= {
                    c(cci(alpha*2)[1], Inf)
                }, less= {
                    c(0, cci(alpha*2)[2])
                })
            }
            attr(cint, "conf.level") <- conf.level
            u <- absigma - qperm(0.5, r, m)
            sgr <- sigma[u <= 0]
            if (length(sgr) == 0) sgr <- NA
            else sgr <- max(sgr)
            sle <- sigma[u > 0]
            if (length(sle) == 0) sle <- NA
            else sle <- min(sle)
            ESTIMATE <- mean(c(sle, sgr))
        }
    }
    else {
        normalize <- function(s, r) {
          ss <- sum(r)
          N <- m + n 
          E <- m/N*ss
          V <- m*(N-m)/(N^2*(N-1))*(N*sum(r^2) - ss^2)
          (s - E)/sqrt(V)
        }
        p <- pnorm(normalize(STATISTIC, r))
        PVAL <- switch(alternative,
                       two.sided = 2 * min(p, 1 - p),
                       less = 1 - p,
                       greater = p)
    
        if(conf.int && !exact) {
            alpha <- 1 - conf.level
            ab <- function(sig, zq = 0) {
                r <- rank(c(x / sig, y))
                s <- sum(pmin(r, N -r + 1)[seq(along = x)])
                normalize(s, r) - zq
            }
            ## Use uniroot here.
            ## Compute the range of sigma first.
            srangepos <- NULL
            srangeneg <- NULL
            if (any(x[x > 0]) && any(y[y > 0]))
                srangepos <-
                    c(min(x[x>0], na.rm=TRUE)/max(y[y>0], na.rm=TRUE), 
                      max(x[x>0], na.rm=TRUE)/min(y[y>0], na.rm=TRUE))
            if (any(x[x <= 0]) && any(y[y < 0]))
                srangeneg <-
                    c(min(x[x<=0], na.rm=TRUE)/max(y[y<0], na.rm=TRUE), 
                      max(x[x<=0], na.rm=TRUE)/min(y[y<0], na.rm=TRUE))
            if (any(is.infinite(c(srangepos, srangeneg)))) {
                warning(paste("Cannot compute asymptotic confidence",
                              "set or estimator"))
                conf.int <- FALSE
            } else {
                ccia <- function(alpha) {
                    ## Check if the statistic exceeds both quantiles
                    ## first.
                    statu <- ab(srange[1], zq=qnorm(alpha/2))
                    statl <- ab(srange[2], zq=qnorm(alpha/2, lower.tail=FALSE))
                    if (statu > 0 || statl < 0) {
                        warning(paste("Samples differ in location:",
                                      "Cannot compute confidence set,",
                                      "returning NA"))
                        return(c(NA, NA))
                    }
                    u <- uniroot(ab, srange, tol=1e-4,
                                 zq=qnorm(alpha/2))$root
                    l <- uniroot(ab, srange, tol=1e-4,
                                 zq=qnorm(alpha/2, lower.tail=FALSE))$root
                    ## The process of the statistics does not need to be
                    ## monotone: sort is ok here.
                    sort(c(u, l))
                }
                srange <- range(c(srangepos, srangeneg), na.rm=FALSE)
                cint <- switch(alternative, two.sided = {
                    ccia(alpha)
                }, greater= {
                    c(ccia(alpha*2)[1], Inf)
                }, less= {
                    c(0, ccia(alpha*2)[2])
                })
                attr(cint, "conf.level") <- conf.level
                ## Check if the statistic exceeds both quantiles first.
                statu <- ab(srange[1], zq=0)
                statl <- ab(srange[2], zq=0)
                if (statu > 0 || statl < 0) {
                    ESTIMATE <- NA
                    warning("Cannot compute estimate, returning NA")
                } else
                    ESTIMATE <- uniroot(ab, srange, tol=1e-4, zq=0)$root
            }
        }
    }
    
    names(STATISTIC) <- "AB"
    RVAL <- list(statistic = STATISTIC,
                 p.value = PVAL,
                 null.value = c("ratio of scales" = 1),
                 alternative = alternative,
                 method = "Ansari-Bradley test",
                 data.name = DNAME)
    if(conf.int)
        RVAL <- c(RVAL,
                  list(conf.int = cint,
                       estimate = c("ratio of scales" = ESTIMATE)))
    class(RVAL) <- "htest"
    return(RVAL)
}

ansari.exact.formula <-
function(formula, data, subset, na.action, ...)
{
    if(missing(formula)
       || (length(formula) != 3)
       || (length(attr(terms(formula[-2]), "term.labels")) != 1)
       || (length(attr(terms(formula[-3]), "term.labels")) != 1))
        stop("formula missing or incorrect")
    if(missing(na.action))
        na.action <- getOption("na.action")
    m <- match.call(expand.dots = FALSE)
    if(is.matrix(eval(m$data, parent.frame())))
        m$data <- as.data.frame(data)
    m[[1]] <- as.name("model.frame")
    m$... <- NULL
    mf <- eval(m, parent.frame())
    DNAME <- paste(names(mf), collapse = " by ")
    names(mf) <- NULL
    response <- attr(attr(mf, "terms"), "response")
    g <- factor(mf[[-response]])
    if(nlevels(g) != 2)
        stop("grouping factor must have exactly 2 levels")
    DATA <- split(mf[[response]], g)
    names(DATA) <- c("x", "y")
    y <- do.call("ansari.exact", c(DATA, list(...)))
    y$data.name <- DNAME
    y
}

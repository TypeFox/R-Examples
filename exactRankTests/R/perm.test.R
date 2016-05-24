# $Id: perm.test.R,v 1.17.4.2 2005/06/02 06:42:59 hothorn Exp $

perm.test <- function(x, ...) UseMethod("perm.test")

perm.test.default <-
function(x,y=NULL, paired = FALSE, alternative = c("two.sided", "less", "greater"),
         mu = 0, exact=NULL, conf.int = FALSE, conf.level = 0.95, tol=NULL, ...) {

    if (is.null(x)) stop("x is missing")
    if (is.null(y)) paired <- TRUE

    alternative <- match.arg(alternative)
    if(!missing(mu) && ((length(mu) > 1) || !is.finite(mu)))
        stop("mu must be a single number")
    if(conf.int) {
        if(!((length(conf.level) == 1)
           && is.finite(conf.level)
           && (conf.level < 1)))
           stop("conf.level must be a single number between 0 and 1")
    }

    if(!is.numeric(x)) stop("`x' must be numeric")

    MIDP <- NULL
    if (!is.null(y)) {
        if(!is.numeric(y)) stop("`y' must be numeric")
        DNAME <- paste(deparse(substitute(x)), "and",
                     deparse(substitute(y)))
    } else DNAME <- deparse(substitute(x))

    if (paired) {
        if (!is.null(y)) {
            OK <- complete.cases(x, y)
            x <- x[OK] - y[OK]
        } else {
            x <- x[complete.cases(x)]
        }
        x <- x - mu
        m <- length(x)
        METHOD <- "1-sample Permutation Test"
        if (is.null(exact)) exact <- (m <= 50)
        if (any(x != floor(x)) && exact) {
            if (is.null(tol)) {
                sc <- (-sum(x <= 0)):sum(x > 0)
                x <- sc[rank(x, ties.method = "first")]
                METHOD <- paste(METHOD, "(scores mapped into 1:m")
            }
            METHOD <- paste(METHOD, "using rounded scores)")
        }

        STATISTIC <- sum(x[x > 0])
        if (exact) {
            PVAL <- switch(alternative,
               "two.sided" = {
                    pperm(STATISTIC, abs(x), m, alternative="two.sided",
                          pprob=TRUE, tol=tol, paired = TRUE)
                }, "greater" = {
                    pperm(STATISTIC, abs(x), m, alternative="greater",
                          pprob=TRUE, tol=tol, paired = TRUE)
                }, "less" = {
                    pperm(STATISTIC, abs(x), m, alternative="less", pprob=TRUE,
                          tol=tol, paired = TRUE)
                })
            MIDP <- PVAL$PPROB
            PVAL <- PVAL$PVALUE
        } else {
            METHOD <- paste("Asymptotic", METHOD)
            wmean <- sum(abs(x))/2
            wvar <- sum(abs(x)^2)/4
            PVAL <- pnorm((STATISTIC - wmean) / sqrt(wvar))
            if (alternative == "greater")
                PVAL <- 1 - PVAL
            if(alternative == "two.sided")
                PVAL <- 2 * min(PVAL, 1-PVAL)
            if (conf.int) { 
                conf.int <- FALSE
                warning("cannot compute asymptotic confidence intervals")
            }
        }
        if (conf.int) {
        warning("Cannot compute confidence interval for paired data!")
        # I do not find any related reference, can R"ohmel 96 be adapted?
        if (FALSE)
            alpha <- 1-conf.level
            xscores <- equiscores(x,m)

            Hx <- rep(0, sum(xscores$scores)*m)

            dummy <- rep(1,m)
            Hx <- .Call("cpermdist2", as.integer(m),
                as.integer(sum(xscores$scores)), as.integer(dummy),
                as.integer(xscores$scores), 
                as.logical(FALSE), PACKAGE="exactRankTests")

            Hx <- matrix(Hx, nrow=m, byrow=TRUE)
            Hx <- rbind(Hx, 0)
            Hx[nrow(Hx), ncol(Hx)] <- 1

            D <- c()
  
            for (i in 1:m) {
                xmean <- which(Hx[i+1,] > 0) - 1
                xm <- length(xmean)
                xmeans <- ((xmean + xscores$add*i)/xscores$fact)/i
                diff <- outer(xmeans, xmeans, "+")
                diff <- sort(diff[!lower.tri(diff)]) / 2
                count <- as.vector(outer(Hx[i+1, xmean+1], Hx[i+1, xmean+1], "*")) 
                od <- order(diff)
                count <- count[od]
                diff <- diff[od]
                du <- duplicated(diff)
                for (k in length(count):1) {
                    if (du[k]) count[k-1] <- count[k-1] + count[k]
                }

                diff <- unique(diff)
                count <- count[!du]

                old <- diff %in% D[,1]
                Dold <- D[,1] %in% diff
 
                if (any(Dold))
                    D[Dold,2] <- D[Dold,2] + count[old]
  
                if (any(!old))
                    D <- rbind(D, cbind(diff[!old], count[!old]))

                D <- D[order(D[,1]),]
            } 

            od <- order(D[,1])
            Dcount <- cumsum(D[od,2])
            Dsort <- D[od,1]

            L <- sum(D[,2])

            cint <- switch(alternative,
                "two.sided" = {
                    qu <- floor((1-alpha/2)*L)             
                    qu <- which(abs(Dcount - qu) == min(abs(Dcount - qu)))
                    ql <- ceiling(alpha/2*L)
                    ql <- which(abs(Dcount - ql) == min(abs(Dcount - ql)))
                    c(Dsort[ql], Dsort[qu])
                }, "greater" = {
                    qu <- floor((1-alpha)*L)
                    qu <- which(abs(Dcount - qu) == min(abs(Dcount - qu)))
                    c(-Inf, Dsort[qu])
                }, "less" = {
                    ql <- ceiling(alpha*L)
                    ql <- which(abs(Dcount - ql) == min(abs(Dcount - ql)))
                    c(Dsort[ql], Inf)
                })
            attr(cint, "conf.level") <- conf.level
        }
#    }
     } else {
        x <- x[complete.cases(x)]
        y <- y[complete.cases(y)]
        m <- length(x)
        n <- length(y)
        if(n < m & conf.int) {
            warning("x has more observations than y, returning perm.test(y, x, ...)")
            return(perm.test(y, x, paired, alternative,
                   mu, exact, conf.int, conf.level, tol=NULL,...))
        }
        x <- x - mu
        if (is.null(exact)) exact <- (m <= 50 && n <= 50)
        cxy <- c(x,y)

        METHOD <- "2-sample Permutation Test"
        if (any(cxy != floor(cxy)) && exact) {
            if (is.null(tol)) {
                cxy <- cxy - min(cxy)
                cxy <- round(cxy*(n+m)/max(cxy))
                METHOD <- paste(METHOD, "(scores mapped into 1:(m+n)")
	        x <- cxy[seq(along=x)]
	        y <- cxy[-seq(along=x)]
            }
            METHOD <- paste(METHOD, "using rounded scores)")
        }

        STATISTIC <- sum(x)
        if (exact) {
            PVAL <- switch(alternative,
                "two.sided" = {
                    pperm(STATISTIC, c(x,y), m, alternative="two.sided",
                          pprob=TRUE, tol=tol)
                }, "greater" = {
                    pperm(STATISTIC, c(x,y), m, alternative="greater",
                          pprob=TRUE, tol=tol)
                }, "less" = {
                    pperm(STATISTIC, c(x,y), m, alternative="less", pprob=TRUE,
                          tol=tol)
                })
            MIDP <- PVAL$PPROB
            PVAL <- PVAL$PVALUE
        } else {
            METHOD <- paste("Asymptotic", METHOD)
            N <- m + n
            wmean <- m/N*sum(c(x,y))
            wvar <- m*n/(N*(N-1))*sum((c(x,y) - wmean/m)^2)
            PVAL <- pnorm((STATISTIC - wmean)/sqrt(wvar))
            if (alternative == "greater")
                PVAL <- 1 - PVAL
            if(alternative == "two.sided")
                PVAL <- 2 * min(PVAL, 1 - PVAL)
            if (conf.int) {
                conf.int <- FALSE
                    warning("cannot compute asymptotic confidence intervals")
            }
        }
        if (conf.int) {
            alpha <- 1-conf.level
            xscores <- equiscores(x,m)
            yscores <- equiscores(y,n)

            tsx <- sum(xscores$scores)
            tsy <- sum(yscores$scores)

            dummy <- rep(1,m)
            Hx <- .Call("cpermdist2", as.integer(m),
                as.integer(tsx), as.integer(dummy),
                as.integer(xscores$scores), 
                as.logical(FALSE), PACKAGE="exactRankTests")

            Hx <- matrix(Hx, nrow=m+1, byrow=TRUE)
           
            dummy <- rep(1,n)
            Hy <- .Call("cpermdist2", as.integer(n),
                as.integer(tsy), as.integer(dummy),
                as.integer(yscores$scores), 
                as.logical(FALSE), PACKAGE="exactRankTests")

            Hy <-  matrix(Hy, nrow=n+1, byrow=TRUE)

            L <- sum(choose(m, 1:m)*choose(n, 1:m))
            D <- c()

            for (i in 1:min(m,n)) {
                xmean <- which(Hx[i+1,] > 0) - 1
                ymean <- which(Hy[i+1,] > 0) - 1
                xm <- length(xmean)
                ym <- length(ymean)

                diff <- as.vector(outer(((xmean + xscores$add*i)/xscores$fact)/i,
                       ((ymean + yscores$add*i)/yscores$fact)/i, "-"))
                count <- as.vector(outer(Hx[i+1, xmean+1], Hy[i+1, ymean+1], "*")) 
                od <- order(diff)
                count <- count[od]
                diff <- diff[od]
                du <- duplicated(diff)
                if (any(is.na(du))) print("NA in du")
                if (length(du) < length(count)) print("du < count")

                for (k in length(count):1) {
                    if (du[k]) count[k-1] <- count[k-1] + count[k]
		}


                diff <- unique(diff)
                count <- count[!du]

                old <- diff %in% D[,1]
                Dold <- D[,1] %in% diff

                if (any(Dold))
                     D[Dold,2] <- D[Dold,2] + count[old]

                if (any(!old))
                    D <- rbind(D, cbind(diff[!old], count[!old]))

                D <- D[order(D[,1]),]
            }
            od <- order(D[,1])
            Dcount <- cumsum(D[od,2])
            Dsort <- D[od,1]
 
            cint <- switch(alternative,
                "two.sided" = {
                    qu <- floor((1-alpha/2)*L)             
                    qu <- which(abs(Dcount - qu) == min(abs(Dcount - qu)))
                    ql <- ceiling(alpha/2*L)
                    ql <- which(abs(Dcount - ql) == min(abs(Dcount - ql)))
                    c(Dsort[ql], Dsort[qu])
                }, "greater" = {
                    qu <- floor((1-alpha)*L)
                    qu <- which(abs(Dcount - qu) == min(abs(Dcount - qu)))
                    c(-Inf, Dsort[qu])
                }, "less" = {
                    ql <- ceiling(alpha*L)
                    ql <- which(abs(Dcount - ql) == min(abs(Dcount - ql)))
                    c(Dsort[ql], Inf)
                })
            attr(cint, "conf.level") <- conf.level
        }
    }
    names(STATISTIC) <- "T"
    if (exact) {
        names(MIDP) <- "point prob"
        RVAL <- list(statistic = STATISTIC,
                     p.value = PVAL,
                     pointprob = MIDP,
                     null.value = c(mu = mu),
                     alternative = alternative,
                     method = METHOD,
                     data.name = DNAME)
        if(conf.int)
             RVAL$conf.int <- cint
    } else {
        RVAL <- list(statistic = STATISTIC,
                     p.value = PVAL,
                     null.value = c(mu = mu),
                     alternative = alternative,
                     method = METHOD,
                     data.name = DNAME)
    }
    class(RVAL) <- "htest"
    return(RVAL)
}

perm.test.formula <-
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
    y <- do.call("perm.test", c(DATA, list(...)))
    y$data.name <- DNAME
    y
}

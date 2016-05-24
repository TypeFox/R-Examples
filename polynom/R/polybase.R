polynomial <-
function(coef = c(0, 1))
{
    a <- as.numeric(coef)
    while((la <- length(a)) > 1 && a[la] == 0) a <- a[-la]
    structure(a, class = "polynomial")
}

as.polynomial <-
function(p)
    if(is.polynomial(p)) p else polynomial(p)
is.polynomial <-
function(p)
    inherits(p, "polynomial")

Ops.polynomial <-
function(e1, e2)
{
    if(missing(e2))
        return(switch(.Generic,
                      "+" = e1,
                      "-" = polynomial(NextMethod(.Generic)),
                      stop("unsupported unary operation")))
    e1 <- unclass(e1)
    e2 <- unclass(e2)
    l1 <- length(e1)
    l2 <- length(e2)
    e1.op.e2 <-
        switch(.Generic,
               "+" = ,
               "-" = {
                   e1 <- c(e1, rep.int(0, max(0, l2 - l1)))
                   e2 <- c(e2, rep.int(0, max(0, l1 - l2)))
                   NextMethod(.Generic)
               },
               "*" = if(l1 == 1 || l2 == 1) e1 * e2 else {
                   m <- outer(e1, e2)
                   as.vector(tapply(m, row(m) + col(m), sum))
               },
               "/" = {
                   if(l2 == 0)
                       stop("unsupported polynomial division")
                   if(l2 == 1)
                       e1 / e2
                   else {
                       p <- rev(e1)
                       q <- rev(e2)
                       r <- rep.int(0, length(p))
                       i <- 0
                       while(length(p) >= l2) {
                           i <- i + 1
                           d <- p[1]/q[1]
                           r[i] <- d
                           p[1:l2] <- p[1:l2] - d * q
                           p <- p[-1]
                       }
                       if(i == 0) 0 else r[i:1]
                   }
               },
               "^" = {
                   if(l2 != 1 || e2 < 0 || e2 %% 1 != 0)
                       stop("unsupported polynomial power")
                   switch(as.character(e2),
                          "0" = 1,
                          "1" = e1,
                      {
                          p <- q <- polynomial(e1)
                          for(i in 2:e2)
                              p <- p * q
                          as.numeric(p)
                      })
               },
               "%%" = {
                   if(l2 == 1)
                       0
                   else {
                       p <- rev(e1)
                       q <- rev(e2)
                       while(length(p) >= l2) {
                           d <- p[1]/q[1]
                           p[1:l2] <- p[1:l2] - d * q
                           p <- p[-1]
                       }
                       if(length(p) == 0) 0 else rev(p)
                   }
               },
               "==" = return(l1 == l2 && all(e1 == e2)),
               "!=" = return(l1 != l2 || any(e1 != e2)),
               stop("unsupported operation on polynomials"))
    polynomial(e1.op.e2)
}

Summary.polynomial <-
function(..., na.rm = FALSE)
{
    ok <- switch(.Generic, sum = , prod = TRUE, FALSE)
    if(!ok)
        stop(gettextf("Generic '%s' not defined for \"%s\" objects.",
                      .Generic, .Class))
    switch(.Generic,
           "sum" = accumulate("+", as.polynomial(0), polylist(...)),
           "prod" = accumulate("*", as.polynomial(1), polylist(...)))
}

Math.polynomial <-
function(x, ...)
{
    switch(.Generic,
           round = ,
           signif = ,
           floor = ,
           ceiling = ,
           trunc = polynomial(NextMethod(.Generic)),
           stop(paste(.Generic, "unsupported for polynomials")))
}

as.character.polynomial <- function(x, decreasing = FALSE, ...)
{
    p <- unclass(x)
    lp <- length(p) - 1
    names(p) <- 0:lp
    p <- p[p != 0]

    if(length(p) == 0) return("0")

    if(decreasing) p <- rev(p)

    signs <- ifelse(p < 0, "- ", "+ ")
    signs[1] <- if(signs[1] == "- ") "-" else ""

    np <- names(p)
    p <- as.character(abs(p))
    p[p == "1" & np != "0"] <- ""

    pow <- paste("x^", np, sep = "")
    pow[np == "0"] <- ""
    pow[np == "1"] <- "x"
    stars <- rep.int("*", length(p))
    stars[p == "" | pow == ""] <- ""
    paste(signs, p, stars, pow, sep = "", collapse = " ")
}

print.polynomial <-
function(x, digits = getOption("digits"), decreasing = FALSE, ...)
{
    p <- as.character.polynomial(signif(x, digits = digits),
                                 decreasing = decreasing)
    pc <- nchar(p)
    ow <- max(35, getOption("width"))
    m2 <- 0
    while(m2 < pc) {
        m1 <- m2 + 1
        m2 <- min(pc, m2 + ow)
        if(m2 < pc)
            while(substring(p, m2, m2) != " " && m2 > m1 + 1)
                m2 <- m2 - 1
        cat(substring(p, m1, m2), "\n")
    }
    invisible(x)
}

as.function.polynomial <-
function(x, ...)
{
    a <- rev(coef(x))
    w <- as.name("w")
    v <- as.name("x")
    ex <- call("{", call("<-", w, 0))
    for(i in seq_along(a)) {
        ex[[i + 2]] <- call("<-", w, call("+", a[1], call("*", v, w)))
        a <- a[-1]
    }
    ex[[length(ex) + 1]] <- w
    f <- function(x) NULL
    body(f) <- ex
    f
}

poly.orth <-
function(x, degree = length(unique(x)) - 1, norm = TRUE)
{
    at <- attr(poly(x, degree), "coefs")
    a <- at$alpha
    N <- at$norm2
    x <- polynomial()
    p <- list(polynomial(0), polynomial(1))
    for(j in 1:degree)
        p[[j + 2]] <-
            (x - a[j]) * p[[j + 1]] - N[j + 1]/N[j] * p[[j]]
    p <- p[-1]
    if(norm) {
        sqrtN <- sqrt(N[-1])
        for(j in 1 + 0:degree) p[[j]] <- p[[j]]/sqrtN[j]
    }
    class(p) <- "polylist"
    p
}

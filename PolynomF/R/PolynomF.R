## revision of the polynomial class with a different representation

polynom <- 
### constructor function
function(a = c(0,1), ..., eps = 0) {
  a <- as.numeric(a)
  if(any(is.na(a))) 
    stop("illegal coefficient vector")
  a[abs(a) < eps] <- 0
  while((k <- length(a)) > 1 && 
    abs(a[k]) <= eps) a <- a[-k]
  structure(function(x) {
    p <- 0
    for(a0 in rev(a)) p <- a0 + x*p
    p
  }, class = "polynom")
}    

as.polynom <- 
### coercion to polynom
function(a)
  if(is.polynom(a)) a else 
    polynom(as.vector(a))

is.polynom <- 
### predicate function
function(a)
  inherits(a, "polynom")

## .poly.mult <- local({
## ### workhorse for convolutions 
##   .pm <- function(e1, e2, l1, l2) {
##     r <- numeric(l1 + l2 -1)
##     ij <- 1:l2
##     for(e in e1) {
##       r[ij] <- r[ij] + e*e2
##       ij <- ij + 1
##     }
##     r
##   }  
##   function(e1, e2, l1, l2) {
##     if(l1 == 1 || l2 == 1) return(e1 * e2) 
##     if(l1 < l2) .pm(e1, e2, l1, l2) else 
##                 .pm(e2, e1, l2, l1)
##     }
## })

.poly.mult <- function(e1, e2, l1, l2) 
   .C("poly_mult",
     as.double(e1), as.integer(l1),
     as.double(e2), as.integer(l2),
     e12 = double(l1+l2-1), PACKAGE = "PolynomF")$e12

.poly.quo <- function(e1, e2, l1, l2) {### quotient
  if(l2 == 0)
    stop("unsupported polynom division")
  if(l2 == 1)
    e1 / e2
  else {
    p <- rev(e1)
    q <- rev(e2)
    r <- numeric(length(p))
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
}

.poly.rem <- function(e1, e2, l1, l2) { ### remainder
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
}
         
Ops.polynom <- function(e1, e2) {                 
  if(missing(e2))   ### unary operations
    return(switch(.Generic,
            "+" = e1,
            "-" = polynom(-coef(e1)),
            stop("unsupported unary operation")))
  e1 <- if(is.polynom(e1)) coef(e1) else as.numeric(e1)
  e2 <- if(is.polynom(e2)) coef(e2) else as.numeric(e2)
  l1 <- length(e1)
  l2 <- length(e2)
  e1.op.e2 <-
    switch(.Generic,
             ### addition and subtraction
         "+" = ,
         "-" = {  
           e1 <- c(e1, rep.int(0, max(0, l2 - l1)))
           e2 <- c(e2, rep.int(0, max(0, l1 - l2)))
           NextMethod(.Generic)
         },
             ### product of two polynomials
         "*" = .poly.mult(e1, e2, l1, l2),
             ### quotient
         "/" =, 
         "%/%" = .poly.quo(e1, e2, l1, l2),
             ### remainder
         "%%" = .poly.rem(e1, e2, l1, l2),
             ### non-negative integer powers
         "^" = { 
           if(l2 != 1 || e2 < 0 || e2 %% 1 != 0)
             stop("unsupported polynom power")
           p <- 1
           while(e2 > 0) {
             if (e2 %% 2 == 1) {
               p <- .poly.mult(p, e1, length(p), l1)
               e2 <- e2 - 1
             }
             e1 <- .poly.mult(e1, e1, l1, l1)
             l1 <- length(e1)
             e2 <- e2 / 2
           }
           p
         },    
             ### equality and inequality
         "==" = return(l1 == l2 && all(e1 == e2)),
         "!=" = return(l1 != l2 || any(e1 != e2)),
         stop("unsupported operation on polynoms"))
  polynom(e1.op.e2)
}

Ops.polylist <- function(e1, e2) {
  if(missing(e2)) 
    return(switch(.Generic,
      "+" = e1,
      "-" = as.polylist(lapply(e1, "-")),
      stop("unknown unary operator!")))
  switch(.Generic,
  "+" =, "-" =, "*" =, "/" =, "%/%" =,
  "%%" = as.polylist(mapply(.Generic, lapply(e1, as.polynom), lapply(e2, as.polynom))),
  "^" = as.polylist(mapply(.Generic, lapply(e1, as.polynom), e2)),
  "==" =,
  "!=" = unlist(mapply(.Generic, lapply(e1, as.polynom), lapply(e2, as.polynom))),
  stop("unsupported operation on polynoms"))
}

.accumulate <- function(f, init, x, right = TRUE) {
## .accumulate a la Abelson and Sussman.
  if(length(x) == 0)
    return(init)
  f <- match.fun(f)
  if(right)
    f(x[[1]], Recall(f, init, x[-1], right = TRUE))
  else
    Recall(f, f(init, x[[1]]), x[-1], right = FALSE)
}

Summary.polynom <- function(..., na.rm = FALSE) {
  ok <- switch(.Generic, sum = , prod = TRUE, FALSE)
  if(!ok)
    stop(gettextf("Generic '%s' not defined for '%s' objects.",
            .Generic, .Class))
  switch(.Generic,
       "sum" = .accumulate("+", polynom(0), polylist(...)),
       "prod" = .accumulate("*", polynom(1), polylist(...)))
}

Summary.polylist <- function(..., na.rm = FALSE) {
  ok <- switch(.Generic, sum = , prod = TRUE, FALSE)
  if(!ok)
    stop(gettextf("Generic '%s' not defined for \"%s\" objects.",
            .Generic, .Class))
  switch(.Generic,
       "sum" = .accumulate("+", polynom(0), c(...)),
       "prod" = .accumulate("*", polynom(1), c(...)))
}

Math.polynom <-
function(x, ...) {
  x <- coef(x)
  switch(.Generic,
       round = ,
       signif = ,
       floor = ,
       ceiling = ,
       trunc = polynom(NextMethod(.Generic)),
       stop(paste(.Generic, "unsupported for polynoms")))
}

Math.polylist <-
  function(x, ...) sapply(x, .Generic, ...)

as.character.polynom <- 
function(x, variable = "x", decreasing = FALSE, ...)
{
  if(is.polynom(x)) p <- coef(x) else p <- unclass(x)
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

  pow <- paste(variable, "^", np, sep = "")
  pow[np == "0"] <- ""
  pow[np == "1"] <- variable
  stars <- rep.int("*", length(p))
  stars[p == "" | pow == ""] <- ""
  paste(signs, p, stars, pow, sep = "", collapse = " ")
}

print.polynom <-
function(x, variable = "x", 
          digits = getOption("digits"), decreasing = FALSE, ...)
{
  x <- coef(x)
  p <- as.character.polynom(signif(x, digits = digits), 
              variable = variable, decreasing = decreasing, ...)
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

as.function.polynom <- 
function (x, variable = "x", ...) {
    a <- rev(coef(x))
    w <- as.name("w")
    v <- as.name(variable)
    ex <- call("{", call("<-", w, 0))
    for (i in seq_along(a)) {
        ex[[i + 2]] <- call("<-", w, 
            call("+", a[1], call("*", v, w)))
        a <- a[-1]
    }
    ex[[length(ex) + 1]] <- w
    f <- quote(function(x) NULL)
    names(f[[2]]) <- variable
    f <- eval(f)
    body(f) <- ex
    f
}
                              
coef.polynom <- function(object,...)
  get("a", envir = environment(object))
  
coef.polylist <- function(object, ...) 
  sapply(object, coef.polynom, ...)

deriv.polynom <- function(expr, ...) {
  stopifnot(require(stats))
  expr <- coef(expr)
  if(length(expr) == 1)
    return(polynom(0))
  expr <- expr[-1]
  polynom(expr * seq(along = expr))
}

integral <- 
function(expr, ...) UseMethod("integral")

integral.default <- function(expr, ...)
stop(gettextf("No 'integral' method for objects of class '%s'",
      class(expr)))

integral.polynom <- function(expr, limits = NULL, ...) {
  expr <- coef(expr)
  p <- polynom(c(0, expr/seq(along = expr)))
  if(is.null(limits))
    p
  else
    diff(p(limits))
}

poly.orth <-
function(x, degree = length(unique(x)) - 1, norm = TRUE)
{
  at <- attr(poly(x, degree), "coefs")
  a <- at$alpha
  N <- at$norm2
  x <- polynom()
  p <- list(polynom(0), polynom(1))
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

.polylist_from_list <-
function(x)
  structure(lapply(x, as.polynom), class = "polylist")

polylist <-
function(...)
  .polylist_from_list(list(...))

is.polylist <-
function(x)
  inherits(x, "polylist")

as.polylist <-
function(x)
{
  if(is.polylist(x)) x
  else if(is.list(x)) .polylist_from_list(x)
  else polylist(x)
}

deriv.polylist <-
function(expr, ...) 
  structure(lapply(expr, deriv), class = class(expr))

integral.polylist <-
function(expr, ...)
{
  result <- lapply(expr, integral, ...)
  if (length(result) > 0 && is.polynom(result[[1]]))
    class(result) <- class(expr)
  result
}

plot.polylist <-
function(x, xlim = 0:1, ylim = range(Px), type = "l", xlab = "x",
          ylab = "P(x)", ..., len = 1000)
{
  p <- x                
  if(missing(xlim)) {
    ## try to cover the "interesting" region
    xlim <- range(Re(unlist(lapply(p, summary.polynom))))
  }
  if(any(is.na(xlim))) {
    warning("summary of polynom fails. Using nominal xlim")
    xlim <- 0:1
  }
  if(diff(xlim) == 0)
    xlim <- xlim + c(-1, 1)/2
  if(length(xlim) > 2)
    x <- xlim
  else {
    eps <- diff(xlim)/100
    xlim <- xlim + c( - eps, eps)
    x <- seq(xlim[1], xlim[2], len = len)
  }
  Px <- unlist(lapply(p, function(Pi, x) Pi(x), x))
  if(!missing(ylim))
    Px[Px < ylim[1]] <- Px[Px > ylim[2]] <- NA
  plot(cbind(x, Px), xlab = xlab, ylab = ylab, type = "n",
     xlim = xlim, ylim = ylim, ...)
  for(i in seq(along = p))
    lines(p[[i]], lty = i, col = i, ...)
  invisible()
}

## print.polylist <- function(x, ...) {
##   cat("List of polynomials:\n")
##   y <- x
##   x <- unclass(x)
##   NextMethod("print", x, ...)
##   invisible(y)
## }

print.polylist <- function(x, ...) {
  y <- x
  x <- unclass(x)
  cat("List of polynomials:\n")
  if(length(x) > 0) {
    nam <- names(x)
    if(is.null(nam)) {
      for(i in 1:length(x)) {
        cat(paste("[[", i, "]]\n", sep=""))
        print(x[[i]], ...)
        cat("\n")
      }
    } else {
      for(n in nam) {
        cat(paste("$\"", n, "\"\n", sep=""))
        print(x[[n]], ...)
        cat("\n")
      }
    }
  } else {
    NextMethod("print", x, ...)
  }  
  invisible(y)
}

c.polylist <-
function(..., recursive = FALSE)
  .polylist_from_list(unlist(lapply(list(...), as.polylist),
                 recursive = FALSE))

"[.polylist" <- function(x, i)
  .polylist_from_list(NextMethod("["))

rep.polylist <- function(x, times, ...)
  .polylist_from_list(NextMethod("rep"))

unique.polylist <-
function(x, incomparables = FALSE, ...)
  .polylist_from_list(NextMethod("unique"))

change.origin <- function(p, o, ...) UseMethod("change.origin")
change.origin.default <- 
  function(p, o, ...) stop("unimplemented method")
change.origin.polynom <- 
  function(p, o, ...) p(polynom() + as.numeric(o)[1])
change.origin.polylist <- function(p, o, ...) 
  structure(lapply(p, change.origin, o = o), class = "polylist")

plot.polynom <- function(x, xlim = 0:1, ylim = range(Px), 
     type = "l", xlab = "x", ylab = "p(x)", ..., len = 1000) {
  p <- x                  
  if(missing(xlim))
    xlim <- range(c(0, Re(unlist(summary(p)))))
  if(any(is.na(xlim))) {
    warning("summary of polynom fails. Using nominal xlim")
    xlim <- 0:1
  }
  if(diff(xlim) == 0)
    xlim <- xlim + c(-1, 1)/2
  if(length(xlim) > 2)
    x <- xlim
  else {
    eps <- diff(xlim)/100
    xlim <- xlim + c(- eps, eps)
    x <- seq(xlim[1], xlim[2], len = len)
  }
  Px <- p(x)
  if(!missing(ylim))
    Px[Px < ylim[1]] <- Px[Px > ylim[2]] <- NA
  plot(x, Px, xlim = xlim, ylim = ylim, type = "n", 
          xlab = xlab, ylab = ylab, ...)
  pu <- par("usr")
  x <- seq(pu[1], pu[2], len = len)
  lines(x, p(x), type = type, ...)
}

lines.polynom <- function(x, ..., len = 1000)  {
  p <- x               
  pu <- par("usr")
  x <- seq(pu[1], pu[2], len = len)
  lines(x, p(x), ...)
}

points.polynom <- function(x, ..., 
    at = seq(pu[1], pu[2], len = len), len = 100)  {
  p <- x               
  pu <- par("usr")
  points(at, p(at), ...)
}

lines.polylist <- function(x, ..., len = 1000)
  for(i in seq(along = x)) lines(x[[i]], col = i, lty = i, len = len)

points.polylist <- function(x, ..., len = 100)
  for(i in seq(along = x)) points(x[[i]], pch = i, col = i, len = len)

poly.calc <- function(x, y, 
      tol = sqrt(.Machine$double.eps), 
      lab = dimnames(y)[[2]]) { 
  if(missing(y)) { ## case 1: polynomial from zeros
    p <- 1
    for(xi in x)
      p <- c(0, p) - c(xi * p, 0)
    return(polynom(p))
  }                ## case 2: Lagrange interpolating polynomial
  if(is.matrix(y)) {
    if(length(x) != nrow(y))
      stop("x and y are inconsistent in size")
    lis <- list()
    if(is.null(lab))
      lab <- paste("p", 1:(dim(y)[2]), sep = "")
    for(i in 1:dim(y)[2])
      lis[[lab[i]]] <- Recall(x, y[, i], tol)
    return(structure(lis, class = "polylist"))
  }
  if(any(toss <- duplicated(x))) {
    crit <- max(tapply(y, x, function(x) diff(range(x))))
    if(crit > tol)
      warning("some duplicated x-points have inconsistent y-values")
    keep <- !toss
    y <- y[keep]
    x <- x[keep]
  }
  if((m <- length(x)) != length(y))
    stop("x and y(x) do not match in length!")
  if(m <= 1)
    return(polynom(y))
  r <- 0
  for(i in 1:m)
    r <- r + (y[i] * coef(Recall(x[ - i])))/prod(x[i] - x[ - i])
  r[abs(r) < tol] <- 0
  polynom(r)
}

poly.from.zeros <- function(...) poly.calc(unlist(list(...)))
poly.from.roots <- poly.from.zeros
poly.from.values <- poly.calc

predict.polynom <- function(object, newdata, ...) object(newdata)
predict.polylist <- function(object, newdata, ...) 
    sapply(object, function(x, .n) x(.n), .n = newdata)

print.summary.polynom <- function(x, ...) {
  cat("\n Summary information for:\n")
  print(attr(x, "originalPolynomial"))
  cat("\n Zeros:\n")
  print(x$zeros)
  cat("\n Stationary points:\n")
  print(x$stationaryPoints)
  cat("\n Points of inflexion:\n")
  print(x$inflexionPoints)
  invisible(x)
}

solve.polynom <- function(a, b, ...) {
  if(!missing(b))
    a <- a - b
  a <- coef(a)
  if(a[1] == 0) {
    z <- rle(a)$lengths[1]
    a <- a[-(1:z)]
    r <- rep(0, z)
  }
  else
    r <- numeric(0)
  switch(as.character(length(a)),
       "0" =,
       "1" = r,
       "2" = sort(c(r,  - a[1]/a[2])),
     {                                       
	   a <- rev(a)
	   a <- (a/a[1])[-1]
	   M <- rbind( - a, cbind(diag(length(a) - 1), 0))     
	   sort(c(r, eigen(M, symmetric = FALSE,
               only.values = TRUE)$values))
     })
}

solve.polylist <- function(a, b, ...) 
  if(!missing(b)) lapply(a, solve.polynom, b) else
    lapply(a, solve.polynom)
    
summary.polynom <-
function(object, ...)
{
  dp <- deriv(object)
  structure(list(zeros = solve(object),
           stationaryPoints = solve(dp),
           inflexionPoints = solve(deriv(dp))),
        class = "summary.polynom",
        originalPolynomial = object)
}

summary.polylist <- 
  function(object, ...) lapply(object, summary.polynom)

.monic <- function(p) {
  a <- coef(p)
  polynom(a/a[length(a)])
}

.degree <-
function(x)
  length(coef(x)) - 1

.effectively_zero <- 
function (p, tolerance = .Machine$double.eps^0.5) 
all(abs(coef(p)) < tolerance)


.GCD2 <-
function(x, y)
{
  if(.effectively_zero(y)) x
  else if(.degree(y) == 0) polynom(1)
  else Recall(y, x %% y)
}

.LCM2 <-
function(x, y)
{
  if(.effectively_zero(x) || .effectively_zero(y))
    return(polynom(0))
  (x / .GCD2(x, y)) * y
}

GCD <- function(...)
  UseMethod("GCD")

GCD.polynom <- function(...) {
  args <- c.polylist(...)
  if(length(args) < 2)
    stop("Need at least two polynoms.")
  .monic(.accumulate(.GCD2, args[[1]], args[-1], FALSE))
}
GCD.polylist <- GCD.polynom
                        
LCM <- function(...)
  UseMethod("LCM")

LCM.polynom <- function(...) {
  args <- c.polylist(...)
  if(length(args) < 2)
    stop("Need at least two polynoms.")
  .monic(.accumulate(.LCM2, args[[1]], args[-1], FALSE))
}

LCM.polylist <- LCM.polynom

as.function.polylist <- function(x, ...) {
  x <- lapply(x, as.function.polynom)
  function(z, ...) sapply(x, function(p) p(z), ...)
}  
  
  


#-----------------------------------------------------------------------------

Polyn <- function(coeffs, base=0) {
  stopifnot(inherits(coeffs, c("polynomial", "numeric", "integer")))
  coeffs <- as.numeric(coeffs)
  stopifnot(inherits(base, c("numeric", "integer")))
  ds <- which(coeffs != 0)
  if(length(ds)>0) {
    p <- coeffs[ds[1]:ds[length(ds)]]
    class(p) <- "Polyn"
    attr(p, "base") <- base+ds[1]-1
    return(p)
  } else {
    p <- 0
    class(p) <- "Polyn"
    attr(p, "base") <- 0
    return(p)  
  }
}

as.Polyn <- function(x,...) UseMethod("as.Polyn")

as.Polyn.Polyn <- function(x, ...) x

as.Polyn.default <- function(x, ...) Polyn(as.numeric(x, ...), ...)

#-----------------------------------------------------------------------------

Pbase <- function(p) UseMethod("Pbase")
Pbase.Polyn <- function(p) attr(p, "base")
Pbase.numeric <- function(p) 0
Pbase.integer <- function(p) 0

#-----------------------------------------------------------------------------

B <- Polyn(1, 1)
A <- Polyn(1, -1)

#-----------------------------------------------------------------------------

as.character.Polyn <- function(x, ..., backward="B", forward="A") {
  p <- x
  d <- attr(p, "base")
  if(identical(p, Polyn(0))) {
    "0"
  } else {
    s <- ""
    for(i in 1:length(p)) {
      di <- d+i-1
      if (p[i]!=0) {
        if(i == 1) {
          if(p[i] < 0) {
            s <- paste(s, "-", sep="")
          }
        } else if(p[i] < 0) {
          s <- paste(s, "-", sep="")
        } else {
          s <- paste(s, "+", sep="")
        }
        if(di == 0) {
          s <- paste(s, as.character(abs(p[i])), sep="")
        } else if(abs(p[i]) != 1) {
          s <- paste(s, as.character(abs(p[i])), "*", sep="")
        }
        if(di < -1) {
          s <- paste(s, forward, "^", as.character(abs(di)), sep="")
        } else if(di == -1) {
          s <- paste(s, forward, sep="")
        } else if(di == 1) {
          s <- paste(s, backward, sep="")
        } else if(di > 1) {
          s <- paste(s, backward, "^", as.character(di), sep="")
        }
      }
    } 
    s
  }
}

print.Polyn <- function(x, ...) {
  print(as.character(x, ...))
  invisible(x)
}

#-----------------------------------------------------------------------------

polynomial_ref <- function(p, r) {
  stopifnot(inherits(p, c("Polyn", "numeric", "integer")))
  b <- Pbase(p)
  stopifnot(b >= r)
  polynomial(c(rep(0, b-r), as.numeric(p)))
}

#-----------------------------------------------------------------------------
# Arithmetic operators

`+.Polyn` <- function(p1, p2) {
  stopifnot(inherits(p1, c("Polyn", "numeric", "integer")))
  stopifnot(inherits(p2, c("Polyn", "numeric", "integer")))
  b = min(Pbase(p1), Pbase(p2))
  Polyn(polynomial_ref(p1, b) + polynomial_ref(p2, b), b)
}

`-.Polyn` <- function(p1, p2) {
  stopifnot(inherits(p1, c("Polyn", "numeric", "integer")))
  stopifnot(inherits(p2, c("Polyn", "numeric", "integer")))
  b = min(Pbase(p1), Pbase(p2))
  Polyn(polynomial_ref(p1, b)-polynomial_ref(p2, b), b)
}

`*.Polyn` <- function(p1, p2) {
  stopifnot(inherits(p1, c("Polyn", "numeric", "integer")))
  stopifnot(inherits(p2, c("Polyn", "numeric", "integer")))
  b <- Pbase(p1) + Pbase(p2)
  Polyn(polynomial(p1)*polynomial(p2), b)
}

`/.Polyn` <- function(p, x) {
  stopifnot(inherits(p, "Polyn"))
  if(inherits(x, "Polyn")) return(Ratio(p, x))
  stopifnot(inherits(x, c("numeric", "integer")))
  stopifnot(length(x)==1)
  Polyn(as.numeric(p)/x, Pbase(p))
}

`^.Polyn` <- function(p, n) {
  stopifnot(inherits(p, "Polyn"))
  stopifnot(inherits(n, c("numeric", "integer")))
  stopifnot(length(n)==1)
  stopifnot(round(n)==n & n>=0)
  if(n==0) Polyn(1, 0)
  else if (n==1) p
  else Polyn(polynomial(p)^n, Pbase(p)*n)
}

`[.Polyn` <- function(p, index, ..., degree) {
  stopifnot(inherits(p, "Polyn"))
  stopifnot(nargs()<=2)
  # stopifnot(missing(index) | missing(degree))
  if(missing(index) & missing(degree)) return(p)
  pnum <- as.numeric(p)
  if(missing(degree)) {
    pnum[index]
  } else {
    b <- attr(p, "base")
    sapply(degree, function(d) { 
      if(d<b | d-b+1>length(p)) 0 else pnum[d-b+1]
    })
  }
}

#-----------------------------------------------------------------------------
# Relational operators

`==.Polyn` <- function(p1, p2) identical(p1, p2)

`!=.Polyn` <- function(p1, p2) !identical(p1, p2)

`<.Polyn` <- function(p1, p2) as.logical(NA)
`<=.Polyn` <- function(p1, p2) as.logical(NA)
`>.Polyn` <- function(p1, p2) as.logical(NA)
`>=.Polyn` <- function(p1, p2) as.logical(NA)

#-----------------------------------------------------------------------------

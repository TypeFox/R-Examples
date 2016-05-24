## due to a change to .C in 2.16.0
## location
.getA.loc <- function(r){
    approx(x = .radius.gitter, y = .A.loc, xout = r, yleft = 1)$y
}
.getb.loc <- function(r){
    approx(x = .radius.gitter, y = .b.loc, xout = r, yleft = Inf)$y
}

## scale
.getA.sc <- function(r){
    approx(x = .radius.gitter, y = .A.sc, xout = r, yleft = 0.5)$y
}
.geta.sc <- function(r){
    approx(x = .radius.gitter, y = .a.sc, xout = r, yleft = 0)$y
}
.getb.sc <- function(r){
    approx(x = .radius.gitter, y = .b.sc, xout = r, yleft = Inf)$y
}

## location and scale
.getA1.locsc <- function(r){
    approx(x = .radius.gitter, y = .A1.locsc, xout = r, yleft = 1)$y
}
.getA2.locsc <- function(r){
    approx(x = .radius.gitter, y = .A2.locsc, xout = r, yleft = 0.5)$y
}
.geta.locsc <- function(r){
    approx(x = .radius.gitter, y = .a.locsc, xout = r, yleft = 0)$y
}
.getb.locsc <- function(r){
    approx(x = .radius.gitter, y = .b.locsc, xout = r, yleft = Inf)$y
}


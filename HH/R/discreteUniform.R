## discrete uniform
pdiscunif <- function(q, size) pmax(0, pmin(1, floor(q)/size))

qdiscunif <- function(p, size) floor(p*size)

ddiscunif <- function(q, size) {
   is.wholenumber <-
      function(x, tol = .Machine$double.eps^0.5)  abs(x - round(x)) < tol
   fq <- floor(q)
   c(0, 1/size)[((fq > 0) & (fq < size+1) & is.wholenumber(q))+1]
}

rdiscunif <- function(n, size) sample.int(size, n, replace=TRUE)

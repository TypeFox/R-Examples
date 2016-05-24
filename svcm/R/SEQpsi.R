##SEQ.psi() computes the matrix product along COLUMNS of A and along
##dimension p of B
##SEQ.psi() is analogon to the rho-function reported in
##Eilers PHC, Currie I and Durban M (2006). Fast and compact smoothing on large multidimensional grids. Computational Statistics & Data Analysis 50, 61-76. doi:10.1016/j.csda.2004.07.008
SEQpsi <- function(A, B, p) {
  
  sa <- dim(A)
  sb <- dim(B)
  n <- length(sb)
  ip <- (1:n)[-p]
  B <- aperm(B, c(p, ip))
  B <- array(B, c(sb[p], prod(sb[ip])))
  C <- A %*% B
  C <- array(C, c(sa[1], sb[ip]))
  C <- aperm(C, order(c(p, ip)))
  invisible(C)
  
}

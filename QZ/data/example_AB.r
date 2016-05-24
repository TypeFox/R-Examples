### Three examples.

exAB1 <- list(
  description = "http://www.nag.com/lapack-ex/node124.html",
  A = matrix(c(-21.10 -22.50i, 53.50 -50.50i, -34.50 +127.50i,   7.50  +0.50i,
                -0.46  -7.78i, -3.50 -37.50i, -15.50  +58.50i, -10.50  -1.50i,
                 4.30  -5.50i, 39.70 -17.10i, -68.50  +12.50i,  -7.50  -3.50i,
                 5.50  +4.40i, 14.40 +43.30i, -32.50  -46.00i, -19.00 -32.50i),
            nrow = 4, byrow = T),
  B = matrix(c(  1.00  -5.00i,  1.60  +1.20i,  -3.00   +0.00i,   0.00  -1.00i,
                 0.80  -0.60i,  3.00  -5.00i,  -4.00   +3.00i,  -2.40  -3.20i,
                 1.00  +0.00i,  2.40  +1.80i,  -4.00   -5.00i,   0.00  -3.00i,
                 0.00  +1.00i, -1.80  +2.40i,   0.00   -4.00i,   4.00  -5.00i),
             nrow = 4, byrow = T)
)

exAB2 <- list(
  description = "http://www.nag.com/lapack-ex/node119.html",
  A = matrix(c(3.9, 12.5, -34.5, -0.5,
               4.3, 21.5, -47.5,  7.5,
               4.3, 21.5, -43.5,  3.5,
               4.4, 26.0, -46.0,  6.0),
             nrow = 4, byrow = T),
  B = matrix(c(1.0, 2.0, -3.0, 1.0,
               1.0, 3.0, -5.0, 4.0,
               1.0, 3.0, -4.0, 3.0,
               1.0, 3.0, -4.0, 4.0),
             nrow = 4, byrow = T)
)

exAB3 <- list(
  description = "http://www.nag.com/numeric/fl/nagdoc_fl23/xhtml/F08/f08yuf.xml",
  S = matrix(c(4.0 +4.0i, 1.0 +1.0i, 1.0 +1.0i, 2.0 -1.0i,
               0.0 +0.0i, 2.0 +1.0i, 1.0 +1.0i, 1.0 +1.0i,
               0.0 +0.0i, 0.0 +0.0i, 2.0 -1.0i, 1.0 +1.0i,
               0.0 +0.0i, 0.0 +0.0i, 0.0 +0.0i, 6.0 -2.0i),
            nrow = 4, byrow = T),
  T = matrix(c(2.0 +0.0i, 1.0 +1.0i, 1.0 +1.0i, 3.0 -1.0i,
               0.0 +0.0i, 1.0 +0.0i, 2.0 +1.0i, 1.0 +1.0i,
               0.0 +0.0i, 0.0 +0.0i, 1.0 +0.0i, 1.0 +1.0i,
               0.0 +0.0i, 0.0 +0.0i, 0.0 +0.0i, 2.0 +0.0i),
            nrow = 4, byrow = T),
  Q = diag(as.complex(1), nrow = 4, ncol = 4),
  Z = diag(as.complex(1), nrow = 4, ncol = 4)
)

exAB4 <- list(
  description = "http://www.nag.com/numeric/fl/nagdoc_fl23/xhtml/F08/f08ygf.xml",
  S = matrix(c(4.0, 1.0, 1.0, 2.0,
               0.0, 3.0, 4.0, 1.0,
               0.0, 1.0, 3.0, 1.0,
               0.0, 0.0, 0.0, 6.0),
             nrow = 4, byrow = T),
  T = matrix(c(2.0, 1.0, 1.0, 3.0,
               0.0, 1.0, 2.0, 1.0,
               0.0, 0.0, 1.0, 1.0,
               0.0, 0.0, 0.0, 2.0),
             nrow = 4, byrow = T),
  Q = diag(1, nrow = 4, ncol = 4),
  Z = diag(1, nrow = 4, ncol = 4)
)

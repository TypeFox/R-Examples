rococo <- function(x, y, similarity=c("linear", "exp", "gauss", "epstol",
                                      "classical"), tnorm="min", r=0)
{
     if (!is.numeric(x) || !is.numeric(y) || length(x) != length(y))
          stop("'x' and 'y' need to be numeric vectors of the same length")

     if (!is.numeric(r) || any(r < 0))
          stop("'r' must be a (vector of) non-negative real number(s)")

     if (missing(similarity))
         similarity <- "linear"

     similarity <- match.arg(similarity, several.ok=TRUE)

     if (length(r) == 1)
         r[2] <- r[1]

     rcorFunc <- ""

     if (is.character(tnorm))
     {
         tnorm <- match.arg(tnorm, c("min", "prod", "lukasiewicz"))

         rcorFunc <- paste0("rcor_", tnorm)
         rcorTestFunc <- paste0("rcor_permtest_", tnorm)
     }
     else if (is.function(tnorm))
     {
          if (length(formals(tnorm)) != 2)
               stop("'tnorm' should be a function of two arguments, ",
                    "e.g. 'tnorm=function(a, b) a * b' for the product t-norm")

          if (abs(tnorm(1, 0.5) - 0.5) > .Machine$double.eps ||
              abs(tnorm(0, 0.25)) > .Machine$double.eps)
               stop("supplied function does not appear to be a valid t-norm")

          if (requireNamespace("compiler", quietly=TRUE))
               tnorm <- compiler::cmpfun(tnorm)
     }
     else
         stop("'tnorm' should be valid string or a function of two arguments, ",
              "e.g. 'tnorm=function(a, b) a * b' for the product t-norm")

     if (length(similarity) > 1 && similarity[1] != similarity[2])
     {
         if (similarity[1] != "classical" && r[1] == 0)
             r[1] <- 0.1 * IQR(x)

         if (similarity[2] != "classical" && r[2] == 0)
             r[2] <- 0.1 * IQR(y)

         xCorFunc <- paste0("rcor_matrix_", similarity[1])
         yCorFunc <- paste0("rcor_matrix_", similarity[2])

         Rx <- .Call(xCorFunc, vx=as.double(x), as.double(r[1]))
         Ry <- .Call(yCorFunc, vx=as.double(y), as.double(r[2]))
     }
     else
     {
         if (similarity[1] != "classical")
         {
             if (r[1] == 0)
                 r[1] <- 0.1 * IQR(x)

             if (r[2] == 0)
                 r[2] <- 0.1 * IQR(y)
         }

         mCorFunc <- paste0("rcor_matrices_", similarity)

         matrices <- .Call(mCorFunc,
                           vx=as.double(x), vy=as.double(y),
                           as.double(r[1]), as.double(r[2]))
         Rx <- matrices$Rx
         Ry <- matrices$Ry
     }

     if (!identical(rcorFunc, ""))
     {
          res <- .Call(rcorFunc, Rx, Ry)
          c <- res$c
          d <- res$d
     }
     else
     {
          c <- sum(mapply(tnorm, Rx, Ry))
          d <- sum(mapply(tnorm, Rx, t(Ry)))
     }

     gamma <- ifelse(c + d == 0, 0, (c - d) / (c + d))

     gamma
}

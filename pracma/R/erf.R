##
##  e r r o r f . R  Error functions (Matlab Style)
##


# Error function
erf     <- function(x) {  # 2*pnorm(sqrt(2)*x)-1
    pchisq(2*x^2,1)*sign(x)
}

# Inverse error function
erfinv  <- function(y) {
        y[abs(y) > 1] <- NA
        sqrt(qchisq(abs(y),1)/2) * sign(y)
}

# Complementary error function
erfc    <- function(x) {  # 1 - erf(x)
    2 * pnorm(-sqrt(2) * x)
}

# Inverse complementary error function
erfcinv <- function(y) {
    y[y < 0 | y > 2] <- NA
    -qnorm(y/2)/sqrt(2)
}

# Scaled complementary error function
erfcx   <- function(x) {
    exp(x^2) * erfc(x)
}

# Complex error function
erfz    <- function(z)
{
    if (is.null(z)) return( NULL )
    else if (!is.numeric(z) && !is.complex(z))
        stop("Argument 'z' must be a numeric or complex scalar or vector.")

    a0 <- abs(z)
    c0 <- exp(-z * z)
    z1 <- ifelse (Re(z) < 0, -z, z) 

	i <- a0 <= 5.8
	work.i <- i
	cer <- rep(NA, length = length(z))
    if ( sum(work.i) > 0) {
        cs <- z1
        cr <- cs
        for (k in 1:120) {
            cr[work.i] <- cr[work.i] * z1[work.i] * z1[work.i]/(k + 0.5)
            cs[work.i] <- cs[work.i] + cr[work.i]
            work.i <- work.i & (abs(cr/cs) >= 1e-15)
	    if (sum(work.i) == 0) break
        }
        cer[i] <- 2 * c0[i] * cs[i]/sqrt(pi)
    }
	work.i <- !i
    if( sum(work.i) > 0) {
        cl <- 1/z1
        cr <- cl
        for (k in 1:13) {
            cr[work.i] <- -cr[work.i] * (k - 0.5)/(z1[work.i] * z1[work.i])
            cl[work.i] <-  cl[work.i] + cr[work.i]
            work.i <- work.i & (abs(cr/cl) >= 1e-15)
	    if (sum(work.i) == 0) break
        }
        cer[!i] <- 1 - c0[!i] * cl[!i]/sqrt(pi)
    }
	cer[ Re(z) < 0] <- -cer[ Re(z) < 0]
    return(cer)
}


# Imaginary error function
erfi <- function(z) {
	-1i * erfz(1i * z)
}


#-- Error function for real values
# erf <- function(x) {
#     eps <- .Machine$double.eps
#     pi <- 3.141592653589793
#     x2 <- x * x
#     if (abs(x) < 3.5) {
#         er <- 1.0
#         r <- 1.0
#         for (k in 1:50) {
#             r <- r * x2 / (k+0.5)
#             er <- er+r
#             if (abs(r) < abs(er)*eps) break
#         }
#         c0 <- 2.0 / sqrt(pi) * x * exp(-x2)
#         err <- c0 * er
#     } else {
#         er <- 1.0
#         r <- 1.0
#         for (k in 1:12) {
#             r<- -r * (k-0.5) / x2
#             er <- er + r
#         }
#         k <- 12+1
#         c0 <- exp(-x2) / (abs(x) * sqrt(pi))
#         err <- 1.0 - c0 * er
#         if (x < 0.0) err <- -err
#     }
#     return(err) 
# }

#-- Error function for complex values
# erfz <- function(z) {
#     if (is.null(z) || length(z) != 1)
#         stop("Argument 'z' must be single complex value.")
# 
#     a0 <- abs(z);
#     c0 <- exp(-z*z)
# 
#     z1 <- if (Re(z) < 0.0) -z else z
# 
#     if(a0 <=  5.8) {
#         cs <- z1
#         cr <- cs
#         for (k in 1:120) {
#             cr <- cr * z1 * z1 / (k+0.5)
#             cs <- cs + cr
#             if (abs(cr/cs) < 1.0e-15) break
#         }
#         cer <- 2.0 * c0 * cs / sqrt(pi)
# 
#     } else {
#         cl <- 1.0 / z1
#         cr <- cl
#         for (k in 1:13) {
#             cr <- -cr * (k-0.5) / (z1 * z1)
#             cl <- cl + cr
#             if (abs(cr/cl) < 1.0e-15) break
#         }
#         cer <- 1.0 - c0 * cl / sqrt(pi)
#     }
# 
#     if(Re(z)< 0.0) cer <- -cer
# 
#     return(cer)
# }


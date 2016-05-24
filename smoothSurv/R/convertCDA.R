###########################################
#### AUTHOR:    Arnost Komarek         ####
####            (2003)                 ####
####                                   ####
#### FILE:      convertCDA.R           ####
####                                   ####
#### FUNCTIONS: c2a                    ####
####            a2c                    ####
####            derivative.expAD       ####
####            find.c                 ####
####            give.c                 ####
####            derivative.cc3         ####
###########################################

### ========================================
### c2a: Compute a coefficients from c's
### ========================================
c2a <- function(ccoef, which.zero = which.max(ccoef), toler = 1e-6){
   ccoef[ccoef < toler] <- toler
   c.zero <- ccoef[which.zero]

   acoef <- log(ccoef/c.zero)
   return(acoef)
}


### ========================================
### a2c: Function to compute c's from a's
### =========================================
a2c <- function(acoef){
   ccoef <- exp(acoef)
   sum.exp.a <- sum(ccoef)
   if (is.nan(sum.exp.a)) return(NULL)
   ccoef <- ccoef/sum.exp.a

   return(ccoef)
}


### =======================================================================================
### derivative.expAD: Function to compute derivatives of non-zero exp(a)'s w.r.t. exp(d)'s
### =======================================================================================
##
##  * there are g - 1 non-zero a's
##  * there are g - 3 d's
##
##  * g - 3 a's are equal to d's
##  * 2 a's are a function of the rest
##  * 1 a is equal to zero
##
## INPUT: knots ....... vector of knots
##        sdspline .... standard deviation of a basis G-spline
##        last.three ... vector with indeces of the three a's which are to be computed from the rest
##             a[last.three[1]] = 0
##             a[last.three[2]] = first function(d's)
##             a[last.three[3]] = second function(d's)
##        all ......... do I want the full  matrix or only two columns w.r.t.
##                      to the two a's

## OUTPUT: Omega.....
##               if all == TRUE, matrix (g - 2) x g (there is one zero column)
##                  all == FALSE, matrix (g - 2) x 2
##                         (the first row is always an intercept)
derivative.expAD <- function(knots, sdspline, last.three, all = TRUE){

   g <- length(knots)
   if (g < 4){
      stop("Too short input 'knots' vector ")
   }

   if (length(last.three) != 3){
      stop("Incorrect 'last.three' parameter ")
   }

   if (sum(last.three %in% (1:g)) != 3){
      stop("Incorrect 'last.three' parameter ")
   }

   if (last.three[1] == last.three[2] || last.three[1] == last.three[3] || last.three[2] == last.three[3]){
      stop("Incorrect 'last.three' parameter ")
   }

   if (sdspline >= 1 || sdspline <= 0){
      stop("Incorrect 'sdspline' parameter ")
   }

   which.zero <- last.three[1]
   l1 <- last.three[2]
   l2 <- last.three[3]
   s02 <- sdspline * sdspline

   if (abs(knots[l1]) < 1e-4){
      stop("Zero reference knot in derivative.expAD ")
   }

   kn2.kn1 <- knots[l2] - knots[l1]
   jsmm <- 1 - s02 + knots[l1]*knots[l2]

## Knots with removed the two ones corresponding to the two special a's
   knotsMin2 <- knots[-c(l1, l2)]

## Index of the zero a in the shorter (by 2) a's sequence
   which.zero2 <- ifelse(which.zero < min(l1, l2),
                         which.zero,
                         ifelse(which.zero < max(l1, l2),
                                which.zero - 1,
                                which.zero - 2))


## Compute the two columns of the resulting matrix corresponding to the two a's
   vec2b <- -(1/kn2.kn1) * (knotsMin2 - knots[l1])
   vec2c <- (1/jsmm) * (1 - s02 + knots[l1] * knotsMin2)
   vec2 <- (vec2b & vec2c)

   vec1a <- knots[l2]/knots[l1]
   vec1d <- (1/knots[l1]) * knotsMin2
   vec1 <- -vec1a * vec2 - vec1d

   int1 <- vec1[which.zero2]
   int2 <- vec2[which.zero2]

   slope1 <- vec1[-which.zero2]
   slope2 <- vec2[-which.zero2]

   slope <- cbind(slope1, slope2)
   intercept <- c(int1, int2)

   Omega <- rbind(intercept, slope)

   if (all){
      leftMat <- rbind(matrix(0, nrow = 1, ncol = g - 3), diag(g - 3))
      zeroCol <- matrix(c(1, rep(0, g - 3)), nrow = g - 2, ncol = 1)
      k <- 1
      kk <- 1
      for (j in 1:g){
         if (j == which.zero){
            Omega <- cbind(Omega, zeroCol)
         }
         else{
            if (j == l1 || j == l2){
               Omega <- cbind(Omega, Omega[,k])
               k <- k + 1
            }
            else{
               Omega <- cbind(Omega, leftMat[,kk])
               kk <- kk + 1
            }
         }
      }
## the first two columns are now obscure -> remove them
      Omega <- Omega[,-c(1,2)]
   }

   return(Omega)
}


### ===================================================================
### find.c: Find mixture proportions that approximate
###         given distribution (dist) by a G-spline mixture with knots
###         and standard deviation sdspline
### ===================================================================
## RETURNS: vector with c coeff.
##          or NULL if problems to find them
find.c <- function(knots, sdspline, dist = "dnorm"){
   nsplines <- length(knots)
   in.knots <- list(x=as.numeric(knots))
   right.side <- do.call(dist, in.knots)
   right.side[abs(right.side) < 1e-10] <- 0
   mus <- matrix(rep(knots, rep(nsplines, nsplines)), nrow = nsplines)
   knotsmat <- matrix(rep(knots, nsplines), nrow = nsplines)
   Cmat <- dnorm(knotsmat, mean=mus, sd=sdspline)
   Cmat <- qr(Cmat, tol = 1e-07)
   if (Cmat$rank == ncol(Cmat$qr)){
        ccoef <- solve(Cmat, right.side)
        tempmin <- min(ccoef[ccoef > 0])
        ccoef[ccoef <= 0] <- tempmin           ## this is not theoretically possible but numerically it is
   }
   else
        ccoef <- NULL

   return(ccoef)
}


### =================================================================
### give.c: Give a vector of all c's satisfying the three constrains
###         if remaining (g-3) c's are given.
### =================================================================
## INPUT: knots ......... knots
##        sdspline ...... standard deviation of a G-spline
##        c.rest ..... remaining g-3 c coefficients
##        last.three... indeces of the three c coefficients which are a function of the rest ones
## RETURN: a vector of all c's
give.c <- function(knots, sdspline, last.three, c.rest)
{
   g <- length(knots)

   if (length(c.rest) != g - 3)
      stop("Incorrect dimension of the input vector ")

   if (length(last.three) != 3) stop("Incorrect 'last.three' parameter ")
   if (sum(last.three %in% 1:g) != 3) stop("Incorrect 'last.three' parameter ")
   if (length(unique(last.three)) != 3) stop("Incorrect 'last.three' parameter ")

   ## Matrix to compute last c's from the first g - 3 ones
   Omega <- derivative.cc3(knots, sdspline, last.three, all = TRUE)

   ## compute all c's
   Omega0 <- Omega[1,]
   Omega1 <- matrix(Omega[2:(g - 2),], nrow = g - 3)
   c.all <- (t(Omega1) %*% c.rest) + Omega0

   return(as.numeric(c.all))
}


### ==========================================================
### derivative.cc3: Derivatives of all c's w.r.t. (g - 3) c's
### ==========================================================
## INPUT: knots ....... knots 
##        sdspline .... standard deviation of a G-spline
##        last.three... indeces of the three c coefficients which are a function of the rest ones
##        all ....... if TRUE, matrix to compute all c's from the first three ones is returned
##                    if FALSE, matrix to compute last three c's from the first three ones is returned
## RETURN: Matrix where the first row is an intercept
##         and remaining (g-3) rows is dc/da
derivative.cc3 <- function(knots, sdspline, last.three, all = TRUE){

   g <- length(knots)
   if (g < 4) stop("Too short input 'knots' vector ")

   if (length(last.three) != 3) stop("Incorrect 'last.three' parameter ")
   if (sum(last.three %in% 1:g) != 3) stop("Incorrect 'last.three' parameter ")
   if (length(unique(last.three)) != 3) stop("Incorrect 'last.three' parameter ")

   last.three <- last.three[order(last.three)]

   s02 = sdspline * sdspline

   ## Indices of the three c coefficients
   l0 <- last.three[1]
   l1 <- last.three[2]
   l2 <- last.three[3]

   ## Compute first the matrix used to compute last three c's from the first g - 3 ones
   Omega = matrix(0, nrow = g - 2, ncol = 3)

   ## 1st row (intercept)
   Omega[1, 1] = (1 - s02 + knots[l2]*knots[l1])/((knots[l2] - knots[l0])*(knots[l1] -knots[l0]))
   Omega[1, 2] = -(1 - s02 + knots[l2]*knots[l0])/((knots[l2] - knots[l1])*(knots[l1] -knots[l0]))
   Omega[1, 3] = 1 - Omega[1, 1] - Omega[1, 2]

   ## the rest  (loop over rows)
   i <- 2
   for (j in 1:g){
       if (j == l0 || j == l1 || j == l2) next;
       Omega[i, 1] = -((knots[l2] - knots[j])*(knots[l1] - knots[j]))/((knots[l2] - knots[l0])*(knots[l1] -knots[l0]));
       Omega[i, 2] = ((knots[l2] - knots[j])*(knots[l0] - knots[j]))/((knots[l2] - knots[l1])*(knots[l1] -knots[l0]));
       Omega[i, 3] = -(1 + Omega[i, 1] + Omega[i, 2])
       i <- i + 1
   }

   ## Add the components to compute all c's from the first g - 3 ones
   if (all){
      leftmat <- rbind(rep(0, g - 3), diag(g - 3))
      OmegaAll <- NULL
      k <- 1; kk <- 1
      for (j in 1:g){
         if (j == l0 || j == l1 || j == l2){
            OmegaAll <- cbind(OmegaAll, Omega[,k])
            k <- k + 1
         }
         else{
            OmegaAll <- cbind(OmegaAll, leftmat[, kk])
            kk <- kk + 1
         }
      }
      Omega <- OmegaAll
   }

   return(Omega)
}


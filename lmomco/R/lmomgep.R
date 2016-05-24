"lmomgep" <-
function(para, byqua=TRUE) {
    if(! are.pargep.valid(para)) return()
    attributes(para$para) <- NULL

    if(byqua) {
       #message("Dispatching to theoLmoms.max.ostat() for the L-moments by way of integration of the quantile function")
       zz <- theoLmoms.max.ostat(para=para, qua=quagep, nmom=6)
       zz$source <- "lmomgep"
       return(zz)
    }
    # It seems important to NOT pull those terms not dependent on "j" out of
    # the summation below for numerical reasons. A lot of sweat has gone into 
    # trying to implement the mathematics here but still might have issues until
    # are.pargep.valid can also support(?) numerical limits of the parameters not
    # just theoretical.
    "GEP.Enn" <- function(para, n=1, eps=1E-12, maxit=500) {
         #message("GEP:Enn ",n)
         attributes(para$para) <- NULL
         B <- 1/para$para[1]; K <- para$para[2]; H <- para$para[3]
         C1 <- log(H) - log(B);  #message("C1: ",C1)
         C <- exp(C1 + lgamma(n+1) - lgamma(n));  #message("C: ",C)
         j <- -1; err <- Inf; PRIOR <- Inf; SJ <- 0
         while(1) {
            j <- j + 1
            F22 <- BEhypergeo(2,2,1,2, H*(j+1))$value;  # message("F22: ", F22)
            # This nKj is a trapped for the warning from lgamma that a 0
            # has been given as the argument, could trap by options()
             nKj <- ifelse((n*K - j) == 0, Inf, lgamma(n*K - j))
            NUM <- lgamma(n*K + 1) - H*(j+1) 
            #message("NUM: ", NUM)
            DEM <- nKj +  lgamma(j+1) + log(n) + (n*K) * log(1-exp(-H))
            #message("DEM: ", DEM)
            TMP <- exp(NUM - DEM + log(F22))
            #message("TMP: ", TMP)
            SJ <- SJ + (-1)^j * TMP; err <- ((SJ - PRIOR)/SJ)^2
            #message("  SJ:", SJ)
            #message(" err:",err)
            if(err < eps) break;
            if(j == maxit) {
               warning("Maximum iterations reached, result might be unreliable")
               break;
            }
            PRIOR <- SJ
         }
         return(list(value=C*SJ, its=(j+1), error=err))
    }
    nmom <- 6
    Enn <- sapply(1:nmom, function(i) { enn <- GEP.Enn(para,n=i);
           ifelse(enn$error > 1E-9, return(NA), return(enn$value)) } )
    L <- rep(NA, nmom)
    for(r in 1:nmom) {
       L[r] <- sum(sapply(1:r,
                 function(k) { ((-1)^(r-k)/k) * choose(r-1,k-1) * choose(r+k-2, r-1)*Enn[k] }))
       if(! is.finite(L[r])) L[r] <- NA
    }
    R <- rep(NA,nmom)
    R[2] <- L[2]/L[1]
    R[3:nmom] <- sapply(3:nmom, function(r) L[r]/L[2])
    z <- list(lambdas = L, ratios = R,
              trim=NULL, leftrim=NULL, rightrim=NULL,
              source = "lmomgep")
    return(z)
}


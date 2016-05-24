"quaaep4" <-
function(f, para, paracheck=TRUE) {

   if(! check.fs(f)) return()
   if(paracheck == TRUE) {
     if(! are.paraep4.valid(para)) return()
   }

   U <- para$para[1]
   A <- para$para[2]
   K <- para$para[3]
   H <- para$para[4]

   H1   <- 1/H
   K2p1 <- (1 + K*K)
   AK   <- A*K
   AoK  <- A/K
   K2p1oKK <- K2p1/(K*K)

   F.of.XI <- cdfaep4(U, para)

   x <- vector(mode="numeric", length=length(f))
   x[f <  F.of.XI] <- U -  AK*(qgamma(K2p1oKK*   f[f <  F.of.XI],  H1, lower.tail=FALSE))^H1
   x[f >= F.of.XI] <- U + AoK*(qgamma(K2p1   *(1-f[f >= F.of.XI]), H1, lower.tail=FALSE))^H1

   names(x) <- NULL
   return(x)
}


Mls  <- function(data, mu = NA, sig = NA, tol = 1e-6, Hessian = FALSE)
{
 # mu is estimate of the mean
 # sig is estimate of the covariance
  if (!is.matrix(data) && class(data) != "orderpattern") {
   cat("Warning: data must have the classes of matrix or orderpattern.\n")
   stop("")
 }
 if (is.matrix(data)) {
   allempty <- which(apply(!is.na(data),1,sum) == 0)
   if (length(allempty) != 0) {
   data <- data[apply(!is.na(data), 1, sum) != 0, ]
   cat("Warning:", length(allempty), "Cases with all variables missing have been removed
         from the data.\n")
   }
  data <- OrderMissing(data)
 }
 if (class(data) == "orderpattern") {
   allempty <- which(apply(!is.na(data$data),1,sum) == 0)
   if (length(allempty) != 0) {
   data <- data$data
   data <- data[apply(!is.na(data), 1, sum) != 0, ]
   cat("Warning:", length(allempty), "Cases with all variables missing have been removed
         from the data.\n")
   data <- OrderMissing(data)
   }
 }
 if (length(data$data)==0)
 {
   cat("Warning: Data is empty")
   stop("")
 }
  if(ncol(data$data)<2)
 {
   cat("Warning: More than 1 variable is required.\n")
   stop("")
 }
 y <- data$data
 patused <- data$patused
 spatcnt <- data$spatcnt
 if (is.na(mu[1])){
   mu <- matrix(0, ncol(y), 1)
   sig <- diag(1, ncol(y))
 }
 itcnt <- 0
 em <- 0
 repeat {
        emtemp <- Sexpect(y, mu, sig, patused, spatcnt)
        ysbar <- emtemp$ysbar
        sstar <- emtemp$sstar
        em <- max(abs(sstar - mu %*% t(mu) - sig), abs(mu - ysbar))
        mu <- ysbar
        sig <- sstar - mu %*% t(mu)
        itcnt <- itcnt + 1
        if(!(em > tol || itcnt < 2)) break()
 }#end repeat
 rownames(mu) <- colnames(y)
 colnames(sig) <- colnames(y)
 if(Hessian)
 {
  templist <- Ddf(y,mu,sig)
  hessian <- templist$dd
  stderror <- templist$se
  return (list(mu = mu, sig = sig, hessian = hessian, stderror = stderror, iteration = itcnt))
 }
 return(list(mu = mu, sig = sig, iteration = itcnt))
}
#------------------------------------------------------------------
Sexpect <- function(y, mu, sig, patused, spatcnt)
{

 n <-  nrow(y)
 pp <- ncol(y)
 sstar <- matrix(0, pp, pp)
 a <- nrow(mu)
 b <- ncol(mu)
 ysbar <- matrix(0, a, b)
 first <- 1
 for (i in 1:length(spatcnt)) {
     ni <- spatcnt[i] - first + 1 
     stemp <- matrix(0, pp, pp)
     indm <- which(is.na(patused[i, ]))
     indo <- which(!is.na(patused[i, ]))
     yo <- matrix(y[first:spatcnt[i], indo], ni, length(indo))
     first <- spatcnt[i] + 1
     muo <- mu[indo]
     mum <- mu[indm]
     sigoo <- sig[indo, indo]
     sigooi <- solve(sigoo)
     soo <- t(yo) %*% yo
     stemp[indo, indo] <- soo
     ystemp <- matrix(0, ni, pp)
     ystemp[, indo] <- yo
     if (length(indm)!= 0) {
        sigmo <- matrix(sig[indm, indo], length(indm), length(indo))
        sigmm <- sig[indm, indm]
        temp1 <- matrix(mum, ni, length(indm), byrow = TRUE)
        temp2 <- yo - matrix(muo, ni, length(indo), byrow = TRUE)
        ym <- temp1 + temp2 %*% sigooi %*% t(sigmo)
        som <- t(yo) %*% ym
        smm <- ni * (sigmm - sigmo %*% sigooi %*% t(sigmo))+ t(ym)%*%ym
        stemp[indo, indm] <- som
        stemp[indm, indo] <- t(som)
        stemp[indm, indm] <- smm
        ystemp[, indm] <- ym
     }# end if 
     sstar <- sstar + stemp;
     if (ni == 1){
        ysbar <- t(ystemp) + ysbar
     }else { 
        ysbar <- apply(ystemp, 2, sum) + ysbar
     }
 }#end for
 ysbar <- (1 / n) * ysbar
 sstar <- (1 / n) * sstar
 sstar <- (sstar + t(sstar))/2
 
 return(list(ysbar = ysbar, sstar = sstar))
}

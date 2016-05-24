Impute <- function(data, mu = NA, sig = NA, imputation.method = "Normal", resid = NA) 
{ # Check if data is not ordered change it to ordered form
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
 if(length(data$data)==0)
 {
   cat("Warning: data is empty")
   return
 }
  if(ncol(data$data)<2)
 {
   cat("More than 1 variable is required.\n")
   stop("")
 }
 y <- data$data
 patused <- data$patused
 spatcnt <- data$spatcnt
 patcnt <- data$patcnt
 g <- data$g
 caseorder <- data$caseorder
 spatcntz <- c(0, spatcnt)
 p <- ncol(y)
 n <- nrow(y)
 yimp <- y
 use.normal <- TRUE

#---------impute the missing data with Servestava method(simple Imputation)--------
 if (imputation.method == "Dist.Free") {
    if (is.na(mu[1])) {
       ybar <- matrix(0, p, 1)
       sbar <- diag(1, p)
       iscomp <- apply(patused, 1, sum, na.rm = TRUE) == p
       
          cind <- which(iscomp)
          ncomp <- patcnt[cind]
          use.normal <- FALSE
          if (ncomp >= 10 && ncomp>=2*p){
             compy <- y[seq(spatcntz[cind] + 1, spatcntz[cind + 1]), ]
             ybar <- matrix(apply(compy, 2, mean))
             sbar <- cov(compy)
             if (is.na(resid[1])){ 
                resid <- (ncomp / (ncomp - 1)) ^ .5 * 
                         (compy - matrix(ybar, ncomp, p, byrow = TRUE))
             }
          } else {
             cat("Warning: There is not sufficient number of complete cases.\n  Dist.free imputation requires a least 10 complete cases\n  or 2*number of variables, whichever is bigger.\n")
             return
          }
       }
      if (!is.na(mu[1])) {
       ybar <- mu
       sbar <- sig
       iscomp <-  apply(patused, 1, sum, na.rm = TRUE) == p
       cind <- which(iscomp)
       ncomp <- patcnt[cind]
       use.normal <- FALSE
       compy <- y[seq(spatcntz[cind] + 1, spatcntz[cind + 1]), ]
       if (is.na(resid[1])){ 
          resid <- (ncomp / (ncomp - 1)) ^ .5 * 
                         (compy - matrix(ybar, ncomp, p, byrow = TRUE))
       }
    }
    indsample <- sample(ncomp, n - ncomp, replace = TRUE)
    resstar <- resid[indsample, ]
    indres1 <- 1
    for (i in 1:g) {
        if (sum(patused[i, ], na.rm = TRUE) != p) # choose a pattern not completely obsered 
        {  test <- y[(spatcntz[i] + 1) : spatcntz[i + 1], ] 
           indres2 <- indres1 + patcnt[i] - 1
           test <- MimputeS(matrix(test,ncol=p), patused[i, ], ybar, sbar,
                            matrix(resstar[indres1:indres2, ],ncol=p))
           indres1 <- indres2 + 1
           
           yimp[(spatcntz[i] + 1) : spatcntz[i + 1], ] <- test
        }
    }#end loop
 }
 if (imputation.method == "Normal" | use.normal) {
    #-----------impute the missing data with normal assumption------------
    if (is.na(mu[1])) {
       emest <- Mls(data, tol = 1e-6)
       mu <- emest$mu
       sig <- emest$sig
    }
    for (i in 1:g) {
        if (sum(patused[i, ], na.rm = TRUE) != p) # choose a pattern not completely obsered 
        {  test <- y[(spatcntz[i] + 1) : spatcntz[i + 1], ] 
           test <- Mimpute(matrix(test,ncol=p), patused[i, ], mu, sig)
           yimp[(spatcntz[i] + 1) : spatcntz[i + 1], ] <- test
        }
    }#end loop
 }
 yimpord <- yimp
 yimp <- yimp[order(caseorder), ]
 imputed <- list(yimp = yimp, yimpOrdered = yimpord, caseorder = caseorder, 
                 patused = patused, patcnt = patcnt)

imputed
}

#--------------------------------------------------------------------------
Mimpute <- function(data, patused, mu, sig) {
# This function imputes the missing data base on multivariate-normal, given the 
# observed and mu and sig for a single pattern it goeas through each patterns and uses the
# conditional distribution of missing given observed and mu and sig to
# impute from the appropriate posterior distribution 
 ni <- nrow(data)
 pp <- ncol(data)
 indm <- which(is.na(patused))
 indo <- which(!is.na(patused))
 pm <- length(indm)
 po <- length(indo)
 muo <- mu[indo]
 mum <- mu[indm]
 sigooi <- solve(sig[indo, indo])
 sigmo <- matrix(sig[indm, indo], pm, po)
 sigmm <- matrix(sig[indm, indm], pm, pm)
 ss1 <- sigmo %*% sigooi
 varymiss <- sigmm - ss1 %*% t(sigmo)
 expymiss <- matrix(mum, ni, pm, byrow = TRUE) + 
             (data[, indo] - matrix(muo, ni, po, byrow = TRUE)) %*% t(ss1)
 if (pm == 1) {
    a <- sqrt(varymiss)
 } else {
    svdvar <- svd(varymiss)
    a <- diag(sqrt(svdvar$d)) %*% t(svdvar$u)
 }
 data[, indm]= matrix(rnorm(ni * pm), ni, pm) %*% a + expymiss
 data
}
#---------------------------------------------------------------------------
MimputeS <- function(data, patused, y1, s1, e)
{
# This function imputes the missing data by Srivastava method,
# given the observed and ybar and s from complete data set
# for a single pattern it goeas through each patterns and uses the
# linear regresion to predict missing given obsrved data and add the
# residual to impute missing data 

 ni <- nrow(data)
 pp <- ncol(data)
 indm <- which(is.na(patused))
 indo <- which(!is.na(patused))
 pm <- length(indm)
 po <- length(indo)
 a <- matrix(s1[indm, indo], pm, po) %*% solve(s1[indo, indo])
 dif <-  data[, indo] - matrix(y1[indo], ni, po, byrow = TRUE)
 z <- matrix(y1[indm], ni, pm, byrow = TRUE) + dif %*% t(a)
 etta <- matrix(e[, indm], ni, pm) - matrix(e[, indo], ni, po) %*% t(a)
 zij <- z + etta
 data[, indm] <- zij 
 data
}



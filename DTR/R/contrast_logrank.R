###################################################
### Referece:
### Guo X: Statistical analysis in two-stage randomization designs in clinical trials. 
### PhD thesis, Department of Statistics, North Carolina State University
###################################################

###################################################
### Referece:
### Feng W, Wahed AS: Supremum weighted log-rank test and sample size for comparing  
### two-stage adaptive treatment strategies. Biometrika 95:695-707, 2008
###################################################

###################################################
### Referece:
### Kidwell KM, Wahed AS: Weighted log-rank statistic to compare shared-path  
### adaptive treatment strategies. Biostatistics, 14(2):299-312, 2013
###################################################

###################################################
### code chunk number 1: chunklibraries
###################################################

require(survival)

###################################################
### code chunk number 2: chunkLogrank
###################################################

contrast_logrank <- function(data # A complete data frame representing the data for two-stage randomization designs
                               # data = data frame {X, TR, R, Z, U, delta}
) {
 
  #Retrieve data
  n <- nrow(data)
  X <- data$X # X=0 for A1, 1 for A2
  TR <- data$TR
  R <- data$R
  Z <- data$Z # Z=0 for B1, 1 for B2
  U <- data$U
  delta <- data$delta
  
  #Chek for errors
  if (is.null(X)) stop("X can not be empty") 
  if (is.null(TR)) stop("TR can not be empty")  
  if (is.null(R)) stop("R can not be empty")
  if (is.null(Z)) stop("Z can not be empty")  
  if (is.null(U)) stop("U can not be empty")  
  if (is.null(delta)) stop("delta can not be empty")
  
  #Estimate probability of being assigned to A2
  pi.x <- sum(X)/n
  #The probability of being assigned to A1 = 1-pi.x
  
  #Known probability of being assigned to B1/B2
  pi.z <- 0.5
    
  #Calculate time-varing weight function for each Ui
  #w11 for A1B1; w12 for A1B2; w21 for A2B1; w22 for A2B2
  #Time as rows, individual as columns
  w11 <- w12 <- w21 <- w22 <- matrix(0, nrow=n, ncol=n)
  
  #Calculate weighted at-risk process Yjk_bar at each Ui
  #Y11_bar for A1B1; Y12_bar for A1B2; Y21_bar for A2B1; Y22_bar for A2B2
  #For example: Y11(Ui)=sum_j{w11j(Ui)*(1-Xj)I(Uj>=Ui)}
  Y11_bar <- Y12_bar <- Y21_bar <- Y22_bar <- rep(0, n)
  
  #Calculate at-risk process for only non-responders
  #Y1_NR for A1; Y2_NR for A2
  Y1_NR <- Y2_NR <- rep(0, n)
  
  #Calculate overall at-risk process for Aj
  #Y1_dot for A1; Y2_dot for A2
  Y1_dot <- Y2_dot <- rep(0, n)
  
  #Calculate overall at-risk process for all patients
  Y <- rep(0, n)
  
  #Calculate sum_p(wjkp(Ui)^2*Yjp(Ui)) for the variance estimator
  s11 <- s12 <- s21 <- s22 <- rep(0, n)
  
  for(i in 1:n) {
    
    #Calculate I(TR <= Ui)
    rind <- as.numeric(TR <= U[i])
    
    #Calculate individual weights for each time Ui
    w11[i,] <- (1-X)*(1 - R*rind + R*rind*(1-Z)/(1-pi.z))/(1-pi.x) # weighting for A1B1
    w12[i,] <- (1-X)*(1 - R*rind + R*rind*Z/pi.z)/(1-pi.x) # weighting for A1B2
    w21[i,] <- X*(1 - R*rind + R*rind*(1-Z)/(1-pi.z))/pi.x # weighting for A2B1
    w22[i,] <- X*(1 - R*rind + R*rind*Z/pi.z)/pi.x # weighting for A2B2
    
    #Calculate Y11_bar(Ui), Y12_bar(Ui), Y21_bar(Ui), Y22_bar(Ui)
    Y11_bar[i] <- sum(w11[i,]*(1-X)*as.numeric(U>=U[i]))
    Y12_bar[i] <- sum(w12[i,]*(1-X)*as.numeric(U>=U[i]))
    Y21_bar[i] <- sum(w21[i,]*X*as.numeric(U>=U[i]))
    Y22_bar[i] <- sum(w22[i,]*X*as.numeric(U>=U[i]))
    
    #Calculate Y1_NR(Ui), Y2_NR(Ui)
    Y1_NR[i] <- sum((1-R*rind)*(1-X)*as.numeric(U>=U[i]))
    Y2_NR[i] <- sum((1-R*rind)*X*as.numeric(U>=U[i]))
    
    #Calculate Y1_dot(Ui), Y2_dot(Ui)
    Y1_dot[i] <- sum((1-X)*as.numeric(U>=U[i]))
    Y2_dot[i] <- sum(X*as.numeric(U>=U[i]))
    
    #Calculate Y(Ui)
    Y[i] <- sum(as.numeric(U>=U[i]))
    
    #Calculate sum_p(wjkp(Ui)^2*Yjp(Ui)) for the variance estimator
    s11[i] <- sum(w11[i,]*w11[i,]*(1-X)*as.numeric(U>=U[i]))
    s12[i] <- sum(w12[i,]*w12[i,]*(1-X)*as.numeric(U>=U[i]))
    s21[i] <- sum(w21[i,]*w21[i,]*X*as.numeric(U>=U[i]))
    s22[i] <- sum(w22[i,]*w22[i,]*X*as.numeric(U>=U[i]))

  }
  
  #Calculate weighted logrank test statistic for H0: A1B1=A1B2
  test_1112 <- sum(Y12_bar*diag(w11)*(1-X)*delta/(Y11_bar+Y12_bar), na.rm=T) -
    sum(Y11_bar*diag(w12)*(1-X)*delta/(Y11_bar+Y12_bar), na.rm=T)
  
  #Calculate weighted logrank test statistic for H0: A1B1=A2B1
  test_1121 <- sum(Y21_bar*diag(w11)*(1-X)*delta/(Y11_bar+Y21_bar), na.rm=T) -
    sum(Y11_bar*diag(w21)*X*delta/(Y11_bar+Y21_bar), na.rm=T)
  
  #Calculate weighted logrank test statistic for H0: A1B1=A2B2
  test_1122 <- sum(Y22_bar*diag(w11)*(1-X)*delta/(Y11_bar+Y22_bar), na.rm=T) -
    sum(Y11_bar*diag(w22)*X*delta/(Y11_bar+Y22_bar), na.rm=T)
  
  #Calculate weighted logrank test statistic for H0: A1B2=A2B1
  test_1221 <- sum(Y21_bar*diag(w12)*(1-X)*delta/(Y12_bar+Y21_bar), na.rm=T) -
    sum(Y12_bar*diag(w21)*X*delta/(Y12_bar+Y21_bar), na.rm=T)
  
  #Calculate weighted logrank test statistic for H0: A1B2=A2B2
  test_1222 <- sum(Y22_bar*diag(w12)*(1-X)*delta/(Y12_bar+Y22_bar), na.rm=T) -
    sum(Y12_bar*diag(w22)*X*delta/(Y12_bar+Y22_bar), na.rm=T)

  #Calculate weighted logrank test statistic for H0: A2B1=A2B2
  test_2122 <- sum(Y22_bar*diag(w21)*X*delta/(Y21_bar+Y22_bar), na.rm=T) -
    sum(Y21_bar*diag(w22)*X*delta/(Y21_bar+Y22_bar), na.rm=T)
  
  #Test for H0: A1B1=A1B2=A2B1=A2B2
  
  est_overall <- c(test_1112, test_1121, test_1122)
  
  sigma_overall <- n * matrix(c(sum((Y12_bar*Y12_bar*s11+Y11_bar*Y11_bar*s12)*delta/((Y11_bar+Y12_bar)*(Y11_bar+Y12_bar)*Y), na.rm=T)/n -
                                  2*sum(Y11_bar*Y12_bar*Y1_NR*delta/((Y11_bar+Y12_bar)*(Y11_bar+Y12_bar)*Y*(1-pi.x)*(1-pi.x)), na.rm=T)/n, 
                                
                                sum(Y12_bar*Y21_bar*s11*delta/((Y11_bar+Y12_bar)*(Y11_bar+Y21_bar)*Y), na.rm=T)/n - 
                                  sum(Y11_bar*Y21_bar*Y1_NR*delta/((Y11_bar+Y12_bar)*(Y11_bar+Y21_bar)*Y*(1-pi.x)*(1-pi.x)), na.rm=T)/n,
                                
                                sum(Y12_bar*Y22_bar*s11*delta/((Y11_bar+Y12_bar)*(Y11_bar+Y22_bar)*Y), na.rm=T)/n - 
                                  sum(Y11_bar*Y22_bar*Y1_NR*delta/((Y11_bar+Y12_bar)*(Y11_bar+Y22_bar)*Y*(1-pi.x)*(1-pi.x)), na.rm=T)/n,
                                
                                sum(Y12_bar*Y21_bar*s11*delta/((Y11_bar+Y12_bar)*(Y11_bar+Y21_bar)*Y), na.rm=T)/n - 
                                  sum(Y11_bar*Y21_bar*Y1_NR*delta/((Y11_bar+Y12_bar)*(Y11_bar+Y21_bar)*Y*(1-pi.x)*(1-pi.x)), na.rm=T)/n,
                                
                                sum((Y21_bar*Y21_bar*s11+Y11_bar*Y11_bar*s21)*delta/((Y11_bar+Y21_bar)*(Y11_bar+Y21_bar)*Y), na.rm=T)/n,
                                
                                sum(Y21_bar*Y22_bar*s11*delta/((Y11_bar+Y21_bar)*(Y11_bar+Y22_bar)*Y), na.rm=T)/n +
                                  sum(Y11_bar*Y11_bar*Y2_NR*delta/((Y11_bar+Y21_bar)*(Y11_bar+Y22_bar)*Y*pi.x*pi.x), na.rm=T)/n,
                                
                                sum(Y12_bar*Y22_bar*s11*delta/((Y11_bar+Y12_bar)*(Y11_bar+Y22_bar)*Y), na.rm=T)/n - 
                                  sum(Y11_bar*Y22_bar*Y1_NR*delta/((Y11_bar+Y12_bar)*(Y11_bar+Y22_bar)*Y*(1-pi.x)*(1-pi.x)), na.rm=T)/n,
                                
                                sum(Y21_bar*Y22_bar*s11*delta/((Y11_bar+Y21_bar)*(Y11_bar+Y22_bar)*Y), na.rm=T)/n +
                                  sum(Y11_bar*Y11_bar*Y2_NR*delta/((Y11_bar+Y21_bar)*(Y11_bar+Y22_bar)*Y*pi.x*pi.x), na.rm=T)/n,
                                
                                sum((Y22_bar*Y22_bar*s11+Y11_bar*Y11_bar*s22)*delta/((Y11_bar+Y22_bar)*(Y11_bar+Y22_bar)*Y), na.rm=T)/n), 
                              nrow=3, ncol=3)
  
  D_overall <- matrix(c(1,-1,0,1,0,-1), nrow=2, ncol=3, byrow=T)
  test_overall <- t(D_overall %*% est_overall) %*% solve(D_overall %*% sigma_overall %*% t(D_overall)) %*% (D_overall %*% est_overall)
  p_overall <- pchisq(test_overall, df=2, lower.tail=F)
    
  #Test for H0: A1B1=A1B2
  #Calculate varince estimate for the test statistic
  sigma_1112 <- sum((Y12_bar*Y12_bar*s11+Y11_bar*Y11_bar*s12)*(1-X)*delta/((Y11_bar+Y12_bar)*(Y11_bar+Y12_bar)*Y1_dot), na.rm=T)/n -
    2*sum(Y11_bar*Y12_bar*Y1_NR*(1-X)*delta/((Y11_bar+Y12_bar)*(Y11_bar+Y12_bar)*Y1_dot*(1-pi.x)*(1-pi.x)), na.rm=T)/n
  #Calculate standardized weighted logrank test statistic
  stest_1112 <- test_1112 / sqrt(sigma_1112*n)  
  #Calculate p-value
  if(stest_1112 >= 0) p_1112 <- 2*pnorm(stest_1112, lower.tail=FALSE) else p_1112 <- 2*pnorm(stest_1112, lower.tail=TRUE)
  
  #Test for H0: A1B1=A2B1
  #Calculate varince estimate for the test statistic
  sigma_1121 <- sum((Y21_bar*Y21_bar*s11+Y11_bar*Y11_bar*s21)*diag(w11)*(1-X)*delta/((Y11_bar+Y21_bar)^3), na.rm=T)/n +
      sum((Y21_bar*Y21_bar*s11+Y11_bar*Y11_bar*s21)*diag(w21)*X*delta/((Y11_bar+Y21_bar)^3), na.rm=T)/n
  #Calculate standardized weighted logrank test statistic
  stest_1121 <- test_1121 / sqrt(sigma_1121*n)  
  #Calculate p-value
  if(stest_1121 >= 0) p_1121 <- 2*pnorm(stest_1121, lower.tail=FALSE) else p_1121 <- 2*pnorm(stest_1121, lower.tail=TRUE)

  #Test for H0: A1B1=A2B2
  #Calculate varince estimate for the test statistic
  sigma_1122 <- sum((Y22_bar*Y22_bar*s11+Y11_bar*Y11_bar*s22)*diag(w11)*(1-X)*delta/((Y11_bar+Y22_bar)^3), na.rm=T)/n +
    sum((Y22_bar*Y22_bar*s11+Y11_bar*Y11_bar*s22)*diag(w22)*X*delta/((Y11_bar+Y22_bar)^3), na.rm=T)/n
  #Calculate standardized weighted logrank test statistic
  stest_1122 <- test_1122 / sqrt(sigma_1122*n)  
  #Calculate p-value
  if(stest_1122 >= 0) p_1122 <- 2*pnorm(stest_1122, lower.tail=FALSE) else p_1122 <- 2*pnorm(stest_1122, lower.tail=TRUE)

  #Test for H0: A1B2=A2B1
  #Calculate varince estimate for the test statistic
  sigma_1221 <- sum((Y21_bar*Y21_bar*s12+Y12_bar*Y12_bar*s21)*diag(w12)*(1-X)*delta/((Y12_bar+Y21_bar)^3), na.rm=T)/n +
    sum((Y21_bar*Y21_bar*s12+Y12_bar*Y12_bar*s21)*diag(w21)*X*delta/((Y12_bar+Y21_bar)^3), na.rm=T)/n
  #Calculate standardized weighted logrank test statistic
  stest_1221 <- test_1221 / sqrt(sigma_1221*n)  
  #Calculate p-value
  if(stest_1221 >= 0) p_1221 <- 2*pnorm(stest_1221, lower.tail=FALSE) else p_1221 <- 2*pnorm(stest_1221, lower.tail=TRUE)

  #Test for H0: A1B2=A2B2
  #Calculate varince estimate for the test statistic
  sigma_1222 <- sum((Y22_bar*Y22_bar*s12+Y12_bar*Y12_bar*s22)*diag(w12)*(1-X)*delta/((Y12_bar+Y22_bar)^3), na.rm=T)/n +
    sum((Y22_bar*Y22_bar*s12+Y12_bar*Y12_bar*s22)*diag(w22)*X*delta/((Y12_bar+Y22_bar)^3), na.rm=T)/n
  #Calculate standardized weighted logrank test statistic
  stest_1222 <- test_1222 / sqrt(sigma_1222*n)  
  #Calculate p-value
  if(stest_1222 >= 0) p_1222 <- 2*pnorm(stest_1222, lower.tail=FALSE) else p_1222 <- 2*pnorm(stest_1222, lower.tail=TRUE)
  
  #Test for H0: A2B1=A2B2
  #Calculate varince estimate for the test statistic
  sigma_2122 <- sum((Y22_bar*Y22_bar*s21+Y21_bar*Y21_bar*s22)*X*delta/((Y21_bar+Y22_bar)*(Y21_bar+Y22_bar)*Y2_dot), na.rm=T)/n -
    2*sum(Y21_bar*Y22_bar*Y2_NR*X*delta/((Y21_bar+Y22_bar)*(Y21_bar+Y22_bar)*Y2_dot*pi.x*pi.x), na.rm=T)/n
  #Calculate standardized weighted logrank test statistic
  stest_2122 <- test_2122 / sqrt(sigma_2122*n)  
  #Calculate p-value
  if(stest_2122 >= 0) p_2122 <- 2*pnorm(stest_2122, lower.tail=FALSE) else p_2122 <- 2*pnorm(stest_2122, lower.tail=TRUE)
  
  #Return results
  TEST <- data.frame(c("A1B1=A1B2=A2B1=A2B2", "A1B1=A1B2", "A1B1=A2B1", "A1B1=A2B2", "A1B2=A2B1", "A1B2=A2B2", "A2B1=A2B2"),
                       round(c(test_overall, stest_1112, stest_1121, stest_1122, stest_1221, stest_1222, stest_2122),4),
                     c(2, rep(1,6)),
                     format.pval(round(c(p_overall, p_1112, p_1121, p_1122, p_1221, p_1222, p_2122),4), eps=0.0001))
  names(TEST) <- c("H0", "(standardized) test statistic", "df", "p")  
  return(TEST)    
    
}  
  

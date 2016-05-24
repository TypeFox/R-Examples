CARE1.print <-
function(z){
z=as.vector(unlist(z))
out1 <- estN.n(z)
out2 <- round(estN.pair(z))
out3 <- rbind(cbind(estN.stat(z),round(estN(z,method="Indep",se=TRUE,nboot=200))),
  cbind(estN.stat(z),round(estN(z,method="HSC",se=TRUE,nboot=200))),
  cbind(estN.stat(z),round(estN(z,method="LSC",se=TRUE,nboot=200))))
rownames(out3) <- c("Nhat-0", "Nhat", "Nhat-1")
out4 <- rbind(round(estN.para(z,estN(z,method="Indep")),2),
  round(estN.para(z,estN(z,method="HSC")),2),
  round(estN.para(z,estN(z,method="LSC")),2))
rownames(out4) <- c("Nhat-0", "Nhat", "Nhat-1")
out=list(out1,out2,out3,out4)


cat("(1) NUMBER OF IDENTIFIED CASES IN EACH LIST: \n")
print(out[[1]])
cat("\n")

cat("(2) ESTIMATES BASED ON ANY PAIR OF SAMPLES: \n")
print(out[[2]]) 
cat("\n")
cat("
Note1: Refer to Seber(1982,pages 59 and 60) for Petersen estimator
       and Chapman estimators as well as s.e formula.
Note2: A log-transformation is used is used to obtain the confidence
       interval so that the lower limit is always greater than the
       number of ascertained. Refer to Chao(1987,Biometrics,43,783-791)
       for the construction of the confidence interval.\n\n")

cat("(3) SAMPLE COVERAGE APPROACH: \n")
print(out[[3]])
cat("\n")
if(Chat(z) <= 0.55)
cat("Warning: The estimated sample coverage(overlapping information)is too 
         low so that Nhat is unstable. Recommend the use of Nhat-0 or Nhat-1.\n")
cat("Parameter estimates: \n")
print(out[[4]]) 
cat("\n")
cat("
Definitions for the sample coverage approach:
M: number of individuals ascertained in at least one list.
D: the average of the number of invididuals listed in the combination 
   of any two lists omitting the other one.
C^: sample coverage estimate, see Equation (14) of Chao and Tsay(1998).
    est: population size estimate.
se: estimated standard error of the population size estimation based on
    bootstrap replications. 
cil: 95% confidence interval lower limit(using a log-transformation).   
ciu: 95% confidence interval upper limit(using a log-transformation).
Nhat: Population size estimate for sufficiently high sample coverage 
      cases, see Equation (20) of Chao and Tsay (1998).
Nhat-1: One-step population size estimate for low sample coverage cases;
        see Equation (2.21) of Chao et al. (1996). This estimator is 
        suggested for use when the estimated se of Nhat is relatively large.
u1,u2,u3: estimated mean probabilities depending on the estimate of N.
r12,r13,r23 etc.: estimated coefficient of covariation(CCV) depending on the 
                  estimate of N.\n ") 
}

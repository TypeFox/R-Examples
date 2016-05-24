ICA.BinCont <- function(Dataset, Surr, True, Treat, Diff.Sigma=TRUE, 
         G_pi_00=seq(from=0, to=1, by=.2), 
         G_rho_01_00=seq(from=0, to=1, by=.2), G_rho_01_01=seq(from=0, to=1, by=.2), 
         G_rho_01_10=seq(from=0, to=1, by=.2), G_rho_01_11=seq(from=0, to=1, by=.2), 
         M=10, Seed=sample(1:100000, size=1)){          
  
  totaal <- M*length(G_pi_00)*length(G_rho_01_00)*length(G_rho_01_10)*length(G_rho_01_01)*length(G_rho_01_11)

  Surr <- (Dataset[,paste(substitute(Surr))])
  True <- (Dataset[,paste(substitute(True))])
  Treat <- Dataset[,paste(substitute(Treat))]
  
if (min(na.exclude(Treat))!=c(-1)) {stop("\nTreatment should be coded as -1=control and 1=experimental treatment.")}
if (na.exclude(max(Treat))!=c(1))  {stop("\nTreatment should be coded as -1=control and 1=experimental treatment.")}
if (length(unique(na.exclude(True)))>2) {stop("\nThe true endpoint should be binary.")}
if (min(na.exclude(True))!=c(0)) {stop("\nThe true endpoint should be coded as 0=no response and 1=response.")}
if (max(na.exclude(True))!=c(1))  {stop("\nThe true endpoint should be coded as 0=no response and 1=response.")}
  
data_no_miss <- data.frame(na.exclude(cbind(Surr, True, Treat)))
data_conttreat <- subset(data_no_miss, Treat=="-1")      
data_exptreat <- subset(data_no_miss, Treat=="1")      
S_0 <- data_conttreat$Surr
S_1 <- data_exptreat$Surr

Surr <- na.exclude(Surr)   #LS
True <- na.exclude(True)   #LS
Treat <- na.exclude(Treat) #LS

pi_punt1 <- mean(data_exptreat$True) # E(T|Z=1)
pi_punt0 <- 1 - pi_punt1 # 1-E(T|Z=1)
pi_1punt <- mean(data_conttreat$True) # E(T|Z=0)
pi_0punt <- 1 - pi_1punt # 1 - E(T|Z=0)


count <- 0
R2_H_all <- PD_OK_all <- pi_00_all <- pi_10_all <- pi_01_all <- pi_11_all <- NULL
G_rho_01_00_all <- G_rho_01_01_all <- G_rho_01_10_all <- G_rho_01_11_all <- NULL
pi_Delta_T_min1_all <- pi_Delta_T_0_all <- pi_Delta_T_1_all <- NULL
pi_0_00_e_all <- pi_0_01_e_all <- pi_0_10_e_all <- pi_0_11_e_all <- NULL 
mu_0_00_all <- mu_0_01_all <- mu_0_10_all <- mu_0_11_all <- NULL
sigma_00_00_all <- sigma_00_01_all <- sigma_00_10_all <- sigma_00_11_all <- NULL
pi_1_00_e_all <- pi_1_01_e_all <- pi_1_10_e_all <- pi_1_11_e_all <- NULL 
mu_1_00_all <- mu_1_01_all <- mu_1_10_all <- mu_1_11_all <- NULL
sigma_11_00_all <- sigma_11_01_all <- sigma_11_10_all <- sigma_11_11_all <- NULL

 for (a in 1: length(G_pi_00)){
  for (b in 1: length(G_rho_01_00)){
   for (c in 1: length(G_rho_01_01)){
    for (d in 1: length(G_rho_01_10)){
     for (e in 1: length(G_rho_01_11)){
       for (i in 1: M){

         count <- count+1
         cat("\n", (count/totaal)*100, "% done. \n", sep="") 
         
G_rho_01_00_hier <- G_rho_01_00[b]
G_rho_01_01_hier <- G_rho_01_01[c]
G_rho_01_10_hier <- G_rho_01_10[d]
G_rho_01_11_hier <- G_rho_01_11[e]

pi_00_hier <- G_pi_00[a]  
pi_10_hier <- pi_punt0 - pi_00_hier #UD
pi_01_hier <- pi_0punt -  pi_00_hier #UD
pi_11_hier <- pi_1punt -  pi_10_hier #UD 


if ((pi_00_hier >= 0 & pi_10_hier >= 0 & pi_01_hier >= 0 & pi_11_hier >= 0) & 
  (pi_00_hier <= 1 & pi_10_hier <= 1 & pi_01_hier <= 1 & pi_11_hier <= 1)){

pi_Delta_T_min1 <- pi_10_hier
pi_Delta_T_0 <- pi_00_hier + pi_11_hier
pi_Delta_T_1 <- pi_01_hier 

Seed <- Seed+1; set.seed(Seed)
mix1 <- invisible(mixtools::normalmixEM(arbvar = Diff.Sigma, x = na.exclude(S_0), 
                              lambda=c(pi_00_hier, pi_01_hier, pi_10_hier, pi_11_hier), k=4, 
                              maxit = 1000000)) 
Seed <- Seed+1; set.seed(Seed)
mix2 <- invisible(mixtools::normalmixEM(arbvar = Diff.Sigma, x = na.exclude(S_1), 
                              lambda=c(pi_00_hier, pi_01_hier, pi_10_hier, pi_11_hier), k=4, 
                              maxit = 1000000))

# mixture components f(S_0)
pi_0_00_e <- mix1$lambda[1]
pi_0_01_e <- mix1$lambda[2] 
pi_0_10_e <- mix1$lambda[3]
pi_0_11_e <- mix1$lambda[4]
mu_0_00 <- mix1$mu[1] 
mu_0_01 <- mix1$mu[2] 
mu_0_10 <- mix1$mu[3] 
mu_0_11 <- mix1$mu[4] 
sigma_00_00 <- mix1$sigma[1]**2  # notatie sigma refereert naar var ipv SD
sigma_00_01 <- mix1$sigma[2]**2 
sigma_00_10 <- mix1$sigma[3]**2 
sigma_00_11 <- mix1$sigma[4]**2 


# mixture components f(S_1)
pi_1_00_e <- mix2$lambda[1]
pi_1_01_e <- mix2$lambda[2] 
pi_1_10_e <- mix2$lambda[3]
pi_1_11_e <- mix2$lambda[4]
mu_1_00 <- mix2$mu[1] 
mu_1_01 <- mix2$mu[2] 
mu_1_10 <- mix2$mu[3] 
mu_1_11 <- mix2$mu[4] 
sigma_11_00 <- mix2$sigma[1]**2 
sigma_11_01 <- mix2$sigma[2]**2 
sigma_11_10 <- mix2$sigma[3]**2 
sigma_11_11 <- mix2$sigma[4]**2


# Check PD covariance matrices 
mat_a <- matrix(c(sigma_00_00, (G_rho_01_00_hier*(sqrt(sigma_00_00*sigma_11_00))), 
                  (G_rho_01_00_hier*(sqrt(sigma_00_00*sigma_11_00))), 
                  sigma_11_00), nrow = 2, byrow = TRUE)
eigen_mat_a <- min(eigen((mat_a))$values)

mat_b <- matrix(c(sigma_00_01, (G_rho_01_01_hier*(sqrt(sigma_00_01*sigma_11_01))), 
                  (G_rho_01_01_hier*(sqrt(sigma_00_01*sigma_11_01))), sigma_11_01), nrow = 2, byrow = TRUE)
eigen_mat_b <- min(eigen((mat_b))$values)

mat_c <- matrix(c(sigma_00_10, (G_rho_01_10_hier*(sqrt(sigma_00_10*sigma_11_10))), 
                  (G_rho_01_10_hier*(sqrt(sigma_00_10*sigma_11_10))), sigma_11_10), nrow = 2, byrow = TRUE)
eigen_mat_c <- min(eigen((mat_c))$values)

mat_d <- matrix(c(sigma_00_11, (G_rho_01_11_hier*(sqrt(sigma_00_11*sigma_11_11))), 
                  (G_rho_01_11_hier*(sqrt(sigma_00_11*sigma_11_11))), sigma_11_11), nrow = 2, byrow = TRUE)
eigen_mat_d <- min(eigen((mat_d))$values)



if ((eigen_mat_a > 0) & (eigen_mat_b > 0) & (eigen_mat_c > 0) & (eigen_mat_d > 0)){

mu_s00 <- mu_1_00 - mu_0_00
mu_s01 <- mu_1_01 - mu_0_01
mu_s10 <- mu_1_10 - mu_0_10
mu_s11 <- mu_1_11 - mu_0_11
omega_1 <- pi_00_hier / (pi_00_hier + pi_11_hier); omega_2 <- 1 - omega_1

sigma_s00 <- sigma_00_00 + sigma_11_00 - (2 * (sqrt(sigma_00_00 * sigma_11_00) * G_rho_01_00_hier)) #var
sigma_s01 <- sigma_00_01 + sigma_11_01 - (2 * (sqrt(sigma_00_01 * sigma_11_01) * G_rho_01_01_hier))
sigma_s10 <- sigma_00_10 + sigma_11_10 - (2 * (sqrt(sigma_00_10 * sigma_11_10) * G_rho_01_10_hier))
sigma_s11 <- sigma_00_11 + sigma_11_11 - (2 * (sqrt(sigma_00_11 * sigma_11_11) * G_rho_01_11_hier))


f_Delta_S <- function(val){
  v <- NA
  v <- (pi_00_hier * dnorm(val, mu_s00, sd = sqrt(sigma_s00))) + 
       (pi_01_hier * dnorm(val, mu_s01, sd = sqrt(sigma_s01))) +
       (pi_10_hier * dnorm(val, mu_s10, sd = sqrt(sigma_s10))) +
       (pi_11_hier * dnorm(val, mu_s11, sd = sqrt(sigma_s11)))
  return(v)
}

# I_10
Seed <- Seed+1; set.seed(Seed)
S1 <- rnorm(n=10000, mean = mu_s10, sd = sqrt(sigma_s10))
I_10 <- mean(log(dnorm(S1, mu_s10, sd = sqrt(sigma_s10)) / f_Delta_S(S1)))

# I_01
Seed <- Seed+1; set.seed(Seed)
S2 <- rnorm(n=10000, mean = mu_s01, sd = sqrt(sigma_s01))
I_01 <- mean(log(dnorm(S2, mu_s01, sd = sqrt(sigma_s01)) / f_Delta_S(S2)))

# I_00
Seed <- Seed+1; set.seed(Seed)
S3 <- rnorm(n=10000, mean = mu_s00, sd = sqrt(sigma_s00))
I_00 <- mean(log(
  ((omega_1 * dnorm(S3, mu_s00, sd = sqrt(sigma_s00))) + 
  ((omega_2 * dnorm(S3, mu_s11, sd = sqrt(sigma_s11))))) / 
    f_Delta_S(S3)))

# I_11
Seed <- Seed+1; set.seed(Seed)
S4 <- rnorm(n=10000, mean = mu_s11, sd = sqrt(sigma_s11))
I_11 <-mean(log(
  ((omega_1 * dnorm(S4, mu_s00, sd = sqrt(sigma_s00))) + 
     ((omega_2 * dnorm(S4, mu_s11, sd = sqrt(sigma_s11))))) / 
    f_Delta_S(S4)))

# R^2_H
I_Delta_T_Delta_S <- (pi_00_hier * I_00) + (pi_01_hier * I_01) + (pi_10_hier * I_10) + (pi_11_hier * I_11)
H_Delta_S <- 
  - ((pi_Delta_T_min1 * log(pi_Delta_T_min1)) + 
     (pi_Delta_T_0 * log(pi_Delta_T_0)) +
     (pi_Delta_T_1 * log(pi_Delta_T_1)))
R2_H <- I_Delta_T_Delta_S / H_Delta_S
R2_H_all <- c(R2_H_all, R2_H)

G_rho_01_00_all <- cbind(G_rho_01_00_all, G_rho_01_00_hier)
G_rho_01_01_all <- cbind(G_rho_01_01_all, G_rho_01_01_hier)
G_rho_01_10_all <- cbind(G_rho_01_10_all, G_rho_01_10_hier)
G_rho_01_11_all <- cbind(G_rho_01_11_all, G_rho_01_11_hier)

pi_00_all <- cbind(pi_00_all, pi_00_hier); pi_10_all <- cbind(pi_10_all, pi_10_hier)
pi_11_all <- cbind(pi_11_all, pi_11_hier); pi_01_all <- cbind(pi_01_all, pi_01_hier)

pi_Delta_T_min1_all <- cbind(pi_Delta_T_min1_all, pi_Delta_T_min1)
pi_Delta_T_0_all <- cbind(pi_Delta_T_0_all, pi_Delta_T_0)
pi_Delta_T_1_all <- cbind(pi_Delta_T_1_all, pi_Delta_T_1)

pi_0_00_e_all <- cbind(pi_0_00_e_all, pi_0_00_e)
pi_0_01_e_all <- cbind(pi_0_01_e_all, pi_0_01_e)
pi_0_10_e_all <- cbind(pi_0_10_e_all, pi_0_10_e)
pi_0_11_e_all <- cbind(pi_0_11_e_all, pi_0_11_e)
mu_0_00_all <- cbind(mu_0_00_all, mu_0_00)
mu_0_01_all <- cbind(mu_0_01_all, mu_0_01)
mu_0_10_all <- cbind(mu_0_10_all, mu_0_10)
mu_0_11_all <- cbind(mu_0_11_all, mu_0_11)
sigma_00_00_all <- cbind(sigma_00_00_all, sigma_00_00)
sigma_00_01_all <- cbind(sigma_00_01_all, sigma_00_01)
sigma_00_10_all <- cbind(sigma_00_10_all, sigma_00_10)
sigma_00_11_all <- cbind(sigma_00_11_all, sigma_00_11)
pi_1_00_e_all <- cbind(pi_1_00_e_all, pi_1_00_e)
pi_1_01_e_all <- cbind(pi_1_01_e_all, pi_1_01_e)
pi_1_10_e_all <- cbind(pi_1_10_e_all, pi_1_10_e)
pi_1_11_e_all <- cbind(pi_1_11_e_all, pi_1_11_e)
mu_1_00_all <- cbind(mu_1_00_all, mu_1_00)
mu_1_01_all <- cbind(mu_1_01_all, mu_1_01)
mu_1_10_all <- cbind(mu_1_10_all, mu_1_10)
mu_1_11_all <- cbind(mu_1_11_all, mu_1_11)
sigma_11_00_all <- cbind(sigma_11_00_all, sigma_11_00)
sigma_11_01_all <- cbind(sigma_11_01_all, sigma_11_01)
sigma_11_10_all <- cbind(sigma_11_10_all, sigma_11_10)
sigma_11_11_all <- cbind(sigma_11_11_all, sigma_11_11)



#flush.console()
#cat("\n \n R2_H =", R2_H, "\n")


} # einde if eigen_mat_a en eigen_mat_b pos def
}  # einde for (a in 1: M) loop
}}}}  # einde for i in 1:length(G_pi_ij) loops
}}

fit <- 
  list(R2_H=R2_H_all, pi_00=as.numeric(pi_00_all), pi_01=as.numeric(pi_01_all), 
       pi_10=as.numeric(pi_10_all), pi_11=as.numeric(pi_11_all), 
       G_rho_01_00=as.numeric(G_rho_01_00_all), 
       G_rho_01_01=as.numeric(G_rho_01_01_all),
       G_rho_01_10=as.numeric(G_rho_01_10_all), 
       G_rho_01_11=as.numeric(G_rho_01_11_all),
       pi_Delta_T_min1=as.numeric(pi_Delta_T_min1_all), 
       pi_Delta_T_0=as.numeric(pi_Delta_T_0_all),
       pi_Delta_T_1=as.numeric(pi_Delta_T_1_all),
       pi_0_00=as.numeric(pi_0_00_e_all),
       pi_0_01=as.numeric(pi_0_01_e_all),
       pi_0_11=as.numeric(pi_0_10_e_all),
       pi_0_11=as.numeric(pi_0_11_e_all), 
       mu_0_00=as.numeric(mu_0_00_all),
       mu_0_01=as.numeric(mu_0_01_all),
       mu_0_10=as.numeric(mu_0_10_all),
       mu_0_11=as.numeric(mu_0_11_all),
       sigma2_00_00=as.numeric(sigma_00_00_all),
       sigma2_00_01=as.numeric(sigma_00_01_all),
       sigma2_00_10=as.numeric(sigma_00_10_all),
       sigma2_00_11=as.numeric(sigma_00_11_all),
       pi_1_00=as.numeric(pi_1_00_e_all),
       pi_1_01=as.numeric(pi_1_01_e_all),
       pi_1_10=as.numeric(pi_1_10_e_all),
       pi_1_11=as.numeric(pi_1_11_e_all), 
       mu_1_00=as.numeric(mu_1_00_all),
       mu_1_01=as.numeric(mu_1_01_all),
       mu_1_10=as.numeric(mu_1_10_all),
       mu_1_11=as.numeric(mu_1_11_all),
       sigma_11_00=as.numeric(sigma_11_00_all),
       sigma2_11_01=as.numeric(sigma_11_01_all),
       sigma2_11_10=as.numeric(sigma_11_10_all),
       sigma2_11_11=as.numeric(sigma_11_11_all),
       Call=match.call())   

class(fit) <- "ICA.BinCont"
fit

}  


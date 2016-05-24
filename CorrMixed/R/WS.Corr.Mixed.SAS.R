WS.Corr.Mixed.SAS <- function(Model, D, Sigma2, Asycov, Rho, Tau2, Alpha=0.05, Time){

  if (missing(Rho)) {Rho <- NULL}
  if (missing(Tau2)) {Tau2 <- NULL}
  if (missing(Asycov)) {Asycov <- NULL}
  Pred.Model3.Loess <- NULL
  
if (Model == "Model 1"){
d <- D
R <- rep(d / (d + Sigma2), times=length(Time))
S <- Sigma2 
R_hier <- R
der_d <- S/((d+S)**2)
der_S <- -(d )/((d+S)**2)
sa <- matrix(data = c(der_d, der_S), ncol=1)
var_rho <- var_p <- t(sa) %*% Asycov %*% (sa)
z <- .5 * log((1+R_hier)/(1 - R_hier))
se_alpha <- 1/(1-R_hier**2) * sqrt(var_rho)
up <- z + (qnorm(c(Alpha/2), mean=0, sd=1, lower.tail=FALSE)*se_alpha)
low <- z - (qnorm(c(Alpha/2), mean=0, sd=1, lower.tail=FALSE)*se_alpha)
CI.Upper <- (exp(2*up)-1)/(exp(2*up)+1)
CI.Lower <- (exp(2*low)-1)/(exp(2*low)+1)
}

if (Model == "Model 2"){
u <- Time; d <- D
tau2 <- T <- Tau2; rho <- Rho; S <- sigma2 <- Sigma2
R <- (d + tau2 * exp((-u**2) / (rho**2)))/
  (d + tau2 + sigma2)

CI.Lower <- CI.Upper <- NULL
for (i in 1: length(R)){
  R_hier <- R[i]
  u <- Time[i]**2
  der_rho <- (T * exp(-u**2 / rho**2))/(d+T+S) * ((2 * (u**2))/rho**3)
  der_S <- -(d + (T * exp(-u**2 / rho**2)))/((d+T+S)**2)
  der_d <- (((1 - exp(-u**2 / rho**2))* T)+S)/((d+T+S)**2)
  der_T <- ((d*(exp(-u**2 / rho**2)-1)) + (S * exp(-u**2 / rho**2))) / ((d+T+S)**2)
  sa <- matrix(data = c(der_d, der_T, der_rho, der_S), ncol=1)
  var_rho <- var_p <- t(sa) %*% Asycov %*% sa
  z <- .5 * log((1+R_hier)/(1 - R_hier))
  se_alpha <- 1/(1-R_hier**2) * sqrt(var_rho)
  up <- z + (qnorm(c(Alpha/2), mean=0, sd=1, lower.tail=FALSE)*se_alpha)
  low <- z - (qnorm(c(Alpha/2), mean=0, sd=1, lower.tail=FALSE)*se_alpha)
  up <- (exp(2*up)-1)/(exp(2*up)+1)
  low <- (exp(2*low)-1)/(exp(2*low)+1)
  CI.Lower <- cbind(CI.Lower, low)
  CI.Upper <- cbind(CI.Upper, up)
 }
CI.Lower <- as.numeric(CI.Lower)
CI.Upper <- as.numeric(CI.Upper)
}


if (Model == "Model 3"){ 
  CI.Lower <- CI.Upper <- NULL
  rho <- Rho
  sigma2 <- Sigma2   
  tau2 <- Tau2
  
  # corrs 
  Time <- Time #LS
  all_cols <- NULL
  for (i in 1: max(Time)){
    t1 <- Time[i]
    R_all <- NULL
    
    for (j in 1: max(Time)){
      t2 <- Time[j]
      u <- t2 - t1
      z_s <- matrix(data = c(1, t1), nrow = 1)
      z_t <- matrix(data = c(1, t2), nrow = 1)
      
      R <- ((z_s %*% D %*% t(z_t)) + (tau2 * exp((-u**2) / (rho**2)))) /
        (sqrt(z_s %*% D %*% t(z_s) + tau2 + sigma2) * sqrt(z_t %*% D %*% t(z_t) + tau2 + sigma2)) 
      
      R_all <- c(R_all, R)
    }
    all_cols <- cbind(all_cols, t(t(R_all)))
  }
  
  R <- Correlatie_matrix <- all_cols
  
  Correlatie_matrix <- 
    matrix(Correlatie_matrix, nrow = dim(Correlatie_matrix)[1])
  alles <- NULL
  for (i in 1: max(Time-1)){
    cor_hier <- Correlatie_matrix[row(Correlatie_matrix) == (col(Correlatie_matrix) - i)]
    erbij <- cbind(i, cor_hier)    
    alles <- rbind(alles, erbij)
    rm(erbij)
  } 
  
}


fit <- 
  list(Model=Model, R=R, Alpha=Alpha, CI.Upper=CI.Upper, CI.Lower=CI.Lower, 
       Time=Time, Call=match.call())   

class(fit) <- "WS.Corr.Mixed.SAS"
fit

}


Fract.Poly <- function(Covariate, Outcome, 
                       S=c(-2,-1,-0.5,0,0.5,1,2,3), Max.M=5, Dataset) {

  if (Max.M > 5){
    cat("The maximum allowed Max.M is 5. The requested Max.M is replaced by Max.M=5\n")
    Max.M <- c(5)
  }
  
  Covariate <- Dataset[,paste(substitute(Covariate))]
  Outcome <- Dataset[,paste(substitute(Outcome))]
  n=length(Outcome)
  Results_m1 <- Results_m2 <- Results_m3 <- Results_m4 <- Results_m5 <- NULL
 
  for (j in 1: Max.M){
  r <- length(S)  #aantal in S
  n <- j    #aantal k
  aantal <- factorial(n+r-1)/
    (factorial(n) * factorial(r-1)) 
  if (aantal > 5000){
  cat("Warning: A total of ", aantal, " models of degree ", i, " can be fitted. Consider using \n", sep= "")
  cat("a more restricted set S or a lower Max.M \n\n", sep= "") }
  }
  
Covariate[Covariate==0] <- 0.00000001
  
#m1  
if (Max.M>=1){
Results_m1 <- NULL
for (i in 1:length(S)){
  if  (S[i] != 0){term1 <- Covariate ** S[i]} 
  if  (S[i] == 0){term1 <- log(Covariate)} 
  fit <- lm(Outcome~term1, data=Dataset)
  AIC <- extractAIC(fit)
  Results_m1 <- rbind(Results_m1, cbind(S[i], AIC[2]))
  colnames(Results_m1) <- c("power1", "AIC") 
 }
}

# m=2
if (Max.M>=2){
Results_m2 <- NULL
for (i in 1:length(S)){
  for (j in 1:length(S)){
    if  (S[i] != 0){term1 <- Covariate ** S[i]} 
    if  (S[i] == 0){term1 <- log(Covariate)} 
    if  (S[j] != 0){term2 <- Covariate ** S[j]}
    if  (S[j] == 0){term2 <- log(Covariate)}
    if  (S[i] == S[j]){term2 <- (Covariate ** S[j]) * log(Covariate)} 
    fit <- 
      lm(Outcome~term1+term2, data=Dataset)
    AIC <- extractAIC(fit)
    Results_m2 <- rbind(Results_m2, cbind(S[i], S[j], AIC[2]))
  }
  colnames(Results_m2) <- c("power1", "power2", "AIC") 
}
}
#m=3
if (Max.M>=3){  
Results_m3 <- NULL
for (i in 1:length(S)){
  for (j in 1:length(S)){
    for (k in 1: length(S)){
      if  (S[i] != 0){term1 <- Covariate ** S[i]} 
      if  (S[i] == 0){term1 <- log(Covariate)} 
      if  (S[j] != 0){term2 <- Covariate ** S[j]}
      if  (S[j] == 0){term2 <- log(Covariate)}
      if  (S[k] != 0){term3 <- Covariate ** S[k]} 
      if  (S[k] == 0){term3 <- log(Covariate)}    
      if  (S[i] == S[j]){term2 <- (Covariate ** S[j]) * log(Covariate)} 
      if  (S[i] == S[k]){term2 <- (Covariate ** S[k]) * log(Covariate)} 
      if  (S[j] == S[k]){term2 <- (Covariate ** S[k]) * log(Covariate)} 
      fit <- lm(Outcome~term1+term2+term3, data=Dataset)
      if (NA %in% fit$coeff == FALSE){
      AIC <- extractAIC(fit)
      Results_m3 <- rbind(Results_m3, cbind(S[i], S[j], S[k], AIC[2]))}
      
    }
  }
  colnames(Results_m3) <- c("power1", "power2", "power3", "AIC")
}
}

# m=4
if (Max.M>=4){
  Results_m4 <- NULL
for (i in 1:length(S)){
  for (j in 1:length(S)){
    for (k in 1: length(S)){
      for (l in 1: length(S)){      
        if  (S[i] != 0){term1 <- Covariate ** S[i]} 
        if  (S[i] == 0){term1 <- log(Covariate)} 
        if  (S[j] != 0){term2 <- Covariate ** S[j]}
        if  (S[j] == 0){term2 <- log(Covariate)}
        if  (S[k] != 0){term3 <- Covariate ** S[k]} 
        if  (S[k] == 0){term3 <- log(Covariate)}    
        if  (S[l] != 0){term4 <- Covariate ** S[l]} 
        if  (S[l] == 0){term4 <- log(Covariate)}    
        if  (S[i] == S[j]){term2 <- (Covariate ** S[j]) * log(Covariate)} 
        if  (S[i] == S[k]){term2 <- (Covariate ** S[k]) * log(Covariate)} 
        if  (S[i] == S[l]){term2 <- (Covariate ** S[l]) * log(Covariate)} 
        if  (S[j] == S[k]){term2 <- (Covariate ** S[k]) * log(Covariate)} 
        if  (S[j] == S[l]){term2 <- (Covariate ** S[l]) * log(Covariate)} 
        if  (S[k] == S[l]){term2 <- (Covariate ** S[l]) * log(Covariate)} 
        fit <- lm(Outcome~term1+term2+term3+term4, data=Dataset)
        
        if (NA %in% fit$coeff == FALSE){
        AIC <- extractAIC(fit)
        Results_m4 <- rbind(Results_m4, cbind(S[i], S[j], S[k], S[l], AIC[2]))}
      }
    }
  }
  colnames(Results_m4) <- c("power1", "power2", "power3", "power4", "AIC")
}
}

if (Max.M>=5){  
Results_m5 <- NULL
for (i in 1:length(S)){
  for (j in 1:length(S)){
    for (k in 1: length(S)){
      for (l in 1: length(S)){      
        for (m in 1: length(S)){      
          if  (S[i] != 0){term1 <- Covariate ** S[i]} 
          if  (S[i] == 0){term1 <- log(Covariate)} 
          if  (S[j] != 0){term2 <- Covariate ** S[j]}
          if  (S[j] == 0){term2 <- log(Covariate)}
          if  (S[k] != 0){term3 <- Covariate ** S[k]} 
          if  (S[k] == 0){term3 <- log(Covariate)}    
          if  (S[l] != 0){term4 <- Covariate ** S[l]} 
          if  (S[l] == 0){term4 <- log(Covariate)}    
          if  (S[m] != 0){term5 <- Covariate ** S[m]} 
          if  (S[m] == 0){term5 <- log(Covariate)}    
          if  (S[i] == S[j]){term2 <- (Covariate ** S[j]) * log(Covariate)} 
          if  (S[i] == S[k]){term2 <- (Covariate ** S[k]) * log(Covariate)} 
          if  (S[i] == S[l]){term2 <- (Covariate ** S[l]) * log(Covariate)} 
          if  (S[i] == S[m]){term2 <- (Covariate ** S[m]) * log(Covariate)} 
          if  (S[j] == S[k]){term2 <- (Covariate ** S[k]) * log(Covariate)} 
          if  (S[j] == S[l]){term2 <- (Covariate ** S[l]) * log(Covariate)} 
          if  (S[j] == S[m]){term2 <- (Covariate ** S[m]) * log(Covariate)} 
          if  (S[k] == S[l]){term2 <- (Covariate ** S[l]) * log(Covariate)} 
          if  (S[k] == S[m]){term2 <- (Covariate ** S[m]) * log(Covariate)} 
          if  (S[l] == S[m]){term2 <- (Covariate ** S[m]) * log(Covariate)} 
          fit <- lm(Outcome~term1+term2+term3+term4+term5, data=Dataset)

          if (NA %in% fit$coeff == FALSE){
            AIC <- extractAIC(fit)
            Results_m5 <- rbind(Results_m5, cbind(S[i], S[j], S[k], S[l], S[m], AIC[2]))}
          
        }
      }
    }
  }
}
colnames(Results_m5) <- c("power1", "power2", "power3", "power4", "power5", "AIC")
}

fit <- list(Results.M1=Results_m1, Results.M2=Results_m2, Results.M3=Results_m3, 
            Results.M4=Results_m4, Results.M5=Results_m5)   

class(fit) <- "Fract.Poly"
fit
}



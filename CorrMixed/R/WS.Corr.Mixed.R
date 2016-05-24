WS.Corr.Mixed <- function(Dataset, Fixed.Part=" ", Random.Part=" ", Correlation=" ", 
                    Id, Time = Time, Model=1, Number.Bootstrap=100, Alpha=.05, Seed=1){

  Time <- Dataset[,paste(substitute(Time))]
  
  Max.Time <- length(unique(Time))
  
  CI.Upper.Loess <- CI.Lower.Loess <- Pred.Model3.Loess.Boot_all <- NULL
  Pred.Model3.Loess=NULL
  CI.Upper.Loess <- CI.Lower.Loess <- NULL

  options(warn=-1)
  
  if (Model == 1){

    Model = "Model 1, Random intercept"
    Tau2=NULL
    Rho=NULL
    model <- lme(fixed = Fixed.Part, random = Random.Part, data = Dataset, na.action = na.omit)
    Coef.Fixed <- model$coefficients$fixed
    Std.Error.Fixed <- sqrt(diag(model$varFix))
    
    D <- d <- nlme::getVarCov(model)[1]
    Sigma2 <- as.numeric(model[6])**2
    AIC <- summary(model)$AIC
    LogLik <- summary(model)$logLik
    
    u <- unique(Time)
    R <- rep(d / (d + Sigma2), times=Max.Time)
    
    if (Number.Bootstrap>0){
    # bootstrap
     all_y <- R
      Dataset <- cbind(Dataset, Dataset[, Id]) 
      colnames(Dataset[dim(Dataset)[2]]) <- c("Id_name")
     all_boot_R <- NULL
      for (i in 1: Number.Bootstrap){
        set.seed(Seed+i)
        boot <- 
          sample(unique(Dataset[,colnames(Dataset)%in%Id]), 
                 size = length(unique(Dataset[,colnames(Dataset)%in%Id])), replace = T)
        sample.boot <- NULL
        for (j in 1: length(boot)){
          Dataset <- data.frame(Dataset)
          samen <- 
            Dataset[Dataset[dim(Dataset)[2]]==boot[j], ]
          unit <- rep(j, dim(samen)[1])
          samen <- cbind(samen, unit)
         sample.boot <- rbind(sample.boot, samen) 
       }
       sample.boot <- 
         sample.boot[, 1: (dim(sample.boot)[2])]   
       
       try(model_hier <- lme(fixed = Fixed.Part, random = (~1 | unit), 
                             data = sample.boot, na.action = na.omit), silent=TRUE)
       
       if (exists("model_hier")==TRUE){
         d_hier <- nlme::getVarCov(model_hier)[1]
         Sigma2_hier <- as.numeric(model_hier[6])**2
         
         u <- unique(Time)
         R_hier <- rep(d_hier / (d_hier + Sigma2_hier), times=Max.Time)
         
         all_boot_R <- rbind(all_boot_R, R_hier)
         
         rm(model_hier)  
       }    
       
       ci <- NULL
       for (i in 1: dim(all_boot_R)[2]){
         ci <- cbind(ci, t(t(quantile(all_boot_R[,i], probs = c(Alpha/2, (1-(Alpha/2))), na.rm = T)))) 
       }
       
       CI.Upper <- ci[2,]
       CI.Lower <- ci[1,]  
       R <- all_y #LS  
   } # end bootstrap CI
  }  
    if (Number.Bootstrap==0){
      CI.Upper <- NULL
      CI.Lower <- NULL
    }
  }
  
  
  if (Model == 2){
  
  Model = "Model 2, Random intercept + serial corr (Gaussian)"
      
  model <- 
    lme(fixed = Fixed.Part, random = Random.Part, correlation = Correlation, data = Dataset, na.action = na.omit)

  Coef.Fixed <- model$coefficients$fixed
  Std.Error.Fixed <- sqrt(diag(model$varFix))
    
  serial <- (coef(model$modelStruct$corStruct, unconstrained=FALSE))
  serial <- matrix(data = serial, nrow=1)
  
  nugget <- serial[2]
  rho <- serial[1]
  d <- nlme::getVarCov(model)[1]
  residual_var <- as.numeric(summary(model)[6]) **2   # sigma**2 + tau**2
  tau2 <- residual_var/ ((nugget/(1-nugget))+1)
  sigma2 <- (nugget/(1-nugget)) * tau2
  D <- d
  Tau2 <- tau2
  Rho <- rho
  Sigma2 <- sigma2
  
  AIC <- summary(model)$AIC
  LogLik <- summary(model)$logLik
  
  # corrs als functie van time lag manueel berekenen
  u <- unique(Time)
  R <- (d + tau2 * exp((-u**2) / (rho**2)))/
    (d + tau2 + sigma2)
  
  # CI 
  # bootstrap
  all_y <- R
  Dataset <- cbind(Dataset, Dataset[, Id]) 
  colnames(Dataset[dim(Dataset)[2]]) <- c("Id_name")
  
  if (Number.Bootstrap>0){
  all_boot_R <- NULL
  for (i in 1: Number.Bootstrap){
    set.seed(Seed+i)
  boot <- 
    sample(unique(Dataset[,colnames(Dataset)%in%Id]), 
           size = length(unique(Dataset[,colnames(Dataset)%in%Id])), replace = T)
  
  sample.boot <- NULL
  for (j in 1: length(boot)){
    Dataset <- data.frame(Dataset)
    samen <- 
      Dataset[Dataset[dim(Dataset)[2]]==boot[j], ]
    unit <- rep(j, dim(samen)[1])
    samen <- cbind(samen, unit)
    sample.boot <- rbind(sample.boot, samen) 
  }
  sample.boot <- 
    sample.boot[, 1: (dim(sample.boot)[2])]  
  
  try(model_hier <- lme(fixed = Fixed.Part, random = (~1 | unit), 
                        correlation = Correlation, 
                        data = sample.boot, na.action = na.omit), silent=TRUE)
  
  if (exists("model_hier")==TRUE){
    serial_hier <- (coef(model_hier$modelStruct$corStruct, unconstrained=FALSE))
    serial_hier <- matrix(data = serial_hier, nrow=1)
    nugget_hier <- serial_hier[2]
    rho_hier <- serial_hier[1]
    d_hier <- nlme::getVarCov(model_hier)[1]
    residual_var_hier <- as.numeric(summary(model_hier)[6]) **2   # sigma**2 + tau**2
    tau2_hier <- residual_var_hier/ ((nugget_hier/(1-nugget_hier))+1)
    sigma2_hier <- (nugget_hier/(1-nugget_hier)) * tau2_hier
    u <- unique(Time)
    R_hier <- (d_hier + tau2_hier * exp((-u**2) / (rho_hier**2)))/
      (d_hier + tau2_hier + sigma2_hier)
    
    all_boot_R <- rbind(all_boot_R, R_hier)
    
    rm(model_hier)  
  }    
  
  ci <- NULL
  for (i in 1: dim(all_boot_R)[2]){
    ci <- cbind(ci, t(t(quantile(all_boot_R[,i], probs = c(Alpha/2, (1-(Alpha/2))), na.rm = T)))) 
  }
  
  CI.Upper <- ci[2,]
  CI.Lower <- ci[1,]  
  
  } # end bootstrap CI
  }
  
  if (Number.Bootstrap==0){
    CI.Upper <- NULL
    CI.Lower <- NULL
  }
  
  pred.loess=NULL
  
  } # einde model = 2
  
  
  if (Model == 3){ 

    Model = "Model 3, Random intercept, slope + serial corr (Gaussian)"
    
    model <- lme(fixed = Fixed.Part, random = Random.Part, correlation = Correlation, 
                 data = Dataset, na.action = na.omit)
    Coef.Fixed <- model$coefficients$fixed
    Std.Error.Fixed <- sqrt(diag(model$varFix))
        
    serial <- (coef(model$modelStruct$corStruct, unconstrained=FALSE))
    serial <- matrix(data = serial, nrow=1)
    
    nugget <- serial[2]
    rho <- serial[1]
    D <- matrix(nlme::getVarCov(model)[1:4], nrow=2)
    residual_var <- as.numeric(summary(model)[6]) **2   
    tau2 <- residual_var/ ((nugget/(1-nugget))+1)
    sigma2 <- (nugget/(1-nugget)) * tau2
    Tau2 <- tau2
    Rho <- rho
    Sigma2 <- sigma2
    
    AIC <- summary(model)$AIC
    LogLik <- summary(model)$logLik
    
    # corrs 
    Time <- unique(Time) #LS
    all_cols <- NULL
    for (i in 1: length(Time)){
      t1 <- Time[i]
      R_all <- NULL
      
      for (j in 1: length(Time)){
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
  for (i in 1: length(Time)-1){
    cor_hier <- Correlatie_matrix[row(Correlatie_matrix) == (col(Correlatie_matrix) - i)]
    erbij <- cbind(i, cor_hier)    
    alles <- rbind(alles, erbij)
    rm(erbij)
  } 
  
  # bootstrap
  if (Number.Bootstrap>0){
  all_y <- R
  Dataset <- cbind(Dataset, Dataset[, Id]) 
  colnames(Dataset[dim(Dataset)[2]]) <- c("Id_name")
  
  all_boot_R <- NULL
  val_sol <- 1 
  for (i in 1: Number.Bootstrap){
    set.seed(Seed+i)
    boot <- 
      sample(unique(Dataset[,colnames(Dataset)%in%Id]), 
             size = length(unique(Dataset[,colnames(Dataset)%in%Id])), replace = T)
    
    sample.boot <- NULL
    for (j in 1: length(boot)){
      Dataset <- data.frame(Dataset)
      samen <- 
        Dataset[Dataset[dim(Dataset)[2]]==boot[j], ]
      unit <- rep(j, dim(samen)[1])
      samen <- cbind(samen, unit)
      sample.boot <- rbind(sample.boot, samen) 
    }
    sample.boot <- 
      sample.boot[, 1: (dim(sample.boot)[2])]  
    
    p_1 <- sub(pattern = Id, replacement = "unit", x = Random.Part)[1]
    p_2 <- sub(pattern = Id, replacement = "unit", x = Random.Part)[2]
    
    try(model_hier <- lme(fixed = Fixed.Part, random = formula(paste(p_1, p_2)), 
                          correlation = Correlation, 
                          data = sample.boot, na.action = na.omit), silent=TRUE)
    
    if (exists("model_hier")==TRUE){
      
      serial_hier <- (coef(model_hier$modelStruct$corStruct, unconstrained=FALSE))
      serial_hier <- matrix(data = serial_hier, nrow=1)
      
      nugget_hier <- serial_hier[2]
      rho_hier <- serial_hier[1]
      D_hier <- 
        matrix(nlme::getVarCov(model_hier)[1:4], nrow=2)
      residual_var_hier <- as.numeric(summary(model_hier)[6]) **2   
      tau2_var_hier <- residual_var_hier/ ((nugget_hier/(1-nugget_hier))+1)
      sigma2_var_hier <- (nugget_hier/(1-nugget_hier)) * tau2_var_hier
      Tau2_var_hier <- tau2_var_hier
      rho_hier <- rho_hier
      Sigma2_var_hier <- sigma2_var_hier
      
      #Time <- 1:(max(Max.Time)+1) #LS
      all_cols <- NULL
      for (im in 1: length(Time)){
        t1 <- Time[im]
        R_all <- NULL
        
        for (j in 1: length(Time)){
          t2 <- Time[j]
          u <- t2 - t1
          z_s <- matrix(data = c(1, t1), nrow = 1)
          z_t <- matrix(data = c(1, t2), nrow = 1)
          
          R_hier <- ((z_s %*% D_hier %*% t(z_t)) + (tau2_var_hier * exp((-u**2) / (rho_hier**2)))) /
            (sqrt(z_s %*% D_hier %*% t(z_s) + tau2_var_hier + sigma2_var_hier) * sqrt(z_t %*% D_hier %*% t(z_t) + tau2_var_hier + sigma2_var_hier)) 
          
          R_all <- c(R_all, R_hier)
        }
        all_cols <- cbind(all_cols, t(t(R_all)))
      }
      all_cols <- cbind(t(t(rep(i, times=dim(all_cols)[1]))), all_cols)
      all_boot_R <- rbind(all_boot_R, all_cols)
            
      rm(model_hier)  
    }    
  }  # einde bootstrap
  
  all_upper <- all_lower <- diag(nrow=dim(Correlatie_matrix)[1]) #LS
  for (k in 1: (dim(all_boot_R)[2]-1)){
    for (m in 1: (dim(Correlatie_matrix)[1])){
      r_rowcol <- NULL
      for (s in 1: Number.Bootstrap){
        try(rm(r_rowcol_val), silent=TRUE)
        dat_hier <- 
          all_boot_R[all_boot_R[,1]==s,]
        try(r_rowcol_val <- dat_hier[m,(k+1)], silent=TRUE)
        if (exists("r_rowcol_val")){
        r_rowcol <- (cbind(r_rowcol, r_rowcol_val))
        }
      }
      upper <- quantile(r_rowcol, probs = c(1-(Alpha/2)), na.rm = T)
      lower <- quantile(r_rowcol, probs = c(Alpha/2), na.rm = T)
      try(all_upper[m,k] <- upper, silent=TRUE)
      try(all_lower[m,k] <- lower, silent=TRUE)
    }
  }
  
  CI.Upper <- all_upper
  CI.Lower <- all_lower  
  }
  
  if (Number.Bootstrap==0){
    CI.Upper <- NULL
    CI.Lower <- NULL
    }
  }

  fit <- 
    list(Model=Model, D=D, Tau2=Tau2, Rho=Rho, Sigma2=Sigma2, AIC=AIC, LogLik=LogLik, 
         R=R, CI.Upper=CI.Upper, CI.Lower=CI.Lower, Alpha=Alpha, 
         Coef.Fixed=Coef.Fixed, Std.Error.Fixed=Std.Error.Fixed,
         Time=unique(Time), Call=match.call())   
  
  class(fit) <- "WS.Corr.Mixed"
  fit
  
}

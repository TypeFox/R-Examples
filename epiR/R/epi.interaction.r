epi.interaction <- function(model, coeff, type = c("RERI", "APAB", "S"), conf.level = 0.95){

   N. <- 1 - ((1 - conf.level)/2)
   z <- qnorm(N., mean = 0, sd = 1)
  
   if(type == "RERI"){
     if (!(class(model)[1] == "glm" & class(model)[2] == "lm") & !(class(model)[1] == "clogit" & class(model)[2] == "coxph"))
       stop("Error: model must be either a glm or coxph object")     

     if(class(model)[1] == "glm" & class(model)[2] == "lm"){
       theta1 <- as.numeric(model$coefficients[coeff[1]])
       theta2 <- as.numeric(model$coefficients[coeff[2]])
       theta3 <- as.numeric(model$coefficients[coeff[3]])
     }
     
     if(class(model)[1] == "clogit" & class(model)[2] == "coxph"){
       theta1 <- as.numeric(model$coefficients[coeff[1]])
       theta2 <- as.numeric(model$coefficients[coeff[2]])
       theta3 <- as.numeric(model$coefficients[coeff[3]])
     }
     
     cov.mat <- vcov(model)
     h1 <- -exp(theta1)
     h2 <- -exp(theta2)
     h3 <- exp(theta3)
     
     reri.var <- (h1^2 * (cov.mat[coeff[1],coeff[1]])) + (h2^2 * (cov.mat[coeff[2],coeff[2]])) + (h3^2 * (cov.mat[coeff[3],coeff[3]])) + (2 * h1 * h2 * cov.mat[coeff[1],coeff[2]]) + (2 * h1 * h3 * cov.mat[coeff[1],coeff[3]]) + (2 * h2 * h3 * cov.mat[coeff[2],coeff[3]])
     reri.se <- sqrt(reri.var)
     
     reri.p <- exp(theta3) - exp(theta1) - exp(theta2) + 1
     reri.l <- reri.p - (z * reri.se)
     reri.u <- reri.p + (z * reri.se)
     
     rval <- data.frame(reri.p, reri.l, reri.u)
     names(rval) <- c("est", "lower", "upper")
   }
   
   if(type == "APAB"){
     if (!(class(model)[1] == "glm" & class(model)[2] == "lm") & !(class(model)[1] == "clogit" & class(model)[2] == "coxph"))
       stop("Error: model must be either a glm or coxph object")      

     if(class(model)[1] == "glm" & class(model)[2] == "lm"){
       theta1 <- as.numeric(model$coefficients[coeff[1]])
       theta2 <- as.numeric(model$coefficients[coeff[2]])
       theta3 <- as.numeric(model$coefficients[coeff[3]])
     }
     
     if(class(model)[1] == "clogit" & class(model)[2] == "coxph"){
       theta1 <- as.numeric(model$coefficients[coeff[1]])
       theta2 <- as.numeric(model$coefficients[coeff[2]])
       theta3 <- as.numeric(model$coefficients[coeff[3]])
     }
     
     cov.mat <- vcov(model)
     h1 <- -exp(theta1 - theta3)
     h2 <- -exp(theta2 - theta3)
     h3 <- (exp(theta1) + exp(theta2) - 1) / exp(theta3)
     
     apab.var <- (h1^2 * (cov.mat[coeff[1],coeff[1]])) + (h2^2 * (cov.mat[coeff[2],coeff[2]])) + (h3^2 * (cov.mat[coeff[3],coeff[3]])) + (2 * h1 * h2 * cov.mat[coeff[1],coeff[2]]) + (2 * h1 * h3 * cov.mat[coeff[1],coeff[3]]) + (2 * h2 * h3 * cov.mat[coeff[2],coeff[3]])
     apab.se <- sqrt(apab.var)
     
     # apab.p <- exp(-theta3) - exp(theta1 - theta3) - exp(theta2 - theta3) + 1
     # Equation 4 (Skrondal 2003):
     apab.p <- (exp(theta3) - exp(theta1) - exp(theta2) + 1) / exp(theta3)
     apab.l <- apab.p - (z * apab.se)
     apab.u <- apab.p + (z * apab.se)
     rval <- data.frame(apab.p, apab.l, apab.u)
     names(rval) <- c("est", "lower", "upper")
   }
   
   
   if(type == "S"){
     if (!(class(model)[1] == "glm" & class(model)[2] == "lm") & !(class(model)[1] == "mle2") & !(class(model)[1] == "clogit" & class(model)[2] == "coxph"))
       stop("Error: model must be either a glm, mle2 or coxph object")     
     
     if(class(model)[1] == "glm" & class(model)[2] == "lm"){
       theta1 <- as.numeric(model$coefficients[coeff[1]])
       theta2 <- as.numeric(model$coefficients[coeff[2]])
       theta3 <- as.numeric(model$coefficients[coeff[3]])
     }
     
     if(class(model)[1] == "clogit" & class(model)[2] == "coxph"){
       theta1 <- as.numeric(model$coefficients[coeff[1]])
       theta2 <- as.numeric(model$coefficients[coeff[2]])
       theta3 <- as.numeric(model$coefficients[coeff[3]])
     }
     
     if(class(model)[1] == "mle2"){
       theta1 <- as.numeric(slot(model, "fullcoef")[coeff[1]])
       theta2 <- as.numeric(slot(model, "fullcoef")[coeff[2]])
       theta3 <- as.numeric(slot(model, "fullcoef")[coeff[3]])
     }
     
     # Calculate S.p:
     S.p <- (exp(theta3) - 1) / (exp(theta1) + exp(theta2) - 2)
     cov.mat <- vcov(model)
     
     # If model type is glm or cph and point estimate of S is negative terminate analysis and advise user to use a linear odds model:
     if(class(model)[1] == "glm" & class(model)[2] == "lm" & S.p < 0){
       message <- paste("Point estimate of synergy index (S) is less than zero (", round(S.p, digits = 2), ").\n  Confidence intervals cannot be calculated using the delta method. Consider re-parameterising as linear odds model.", sep = "")
       stop(message)
     }

     if(class(model)[1] == "clogit" & class(model)[2] == "coxph" & S.p < 0){
       message <- paste("Point estimate of synergy index (S) is less than zero (", round(S.p, digits = 2), ").\n  Confidence intervals cannot be calculated using the delta method. Consider re-parameterising as linear odds model.", sep = "")
       stop(message)
     }
     
     # Use delta method (Hosmer and Lemeshow 1992) if model type it glm or cph:
     if(class(model)[1] == "glm" & class(model)[2] == "lm" & S.p > 0){
       ha <- -exp(theta1) / (exp(theta1) + exp(theta2) - 2)
       hb <- -exp(theta2) / (exp(theta1) + exp(theta2) - 2)
       hab <- exp(theta3) / (exp(theta3) - 1)
       
       lnS.var <- (ha^2 * (cov.mat[coeff[1],coeff[1]])) + (hb^2 * (cov.mat[coeff[2],coeff[2]])) + (hab^2 * (cov.mat[coeff[3],coeff[3]])) + (2 * ha * hb * cov.mat[coeff[1],coeff[2]]) + (2 * ha * hab * cov.mat[coeff[1],coeff[3]]) + (2 * hb * hab * cov.mat[coeff[2],coeff[3]])
       lnS.se <- sqrt(lnS.var)
       
       lnS.p <- log((exp(theta3) - 1)) - log((exp(theta1) + exp(theta2) - 2))
       lnS.l <- lnS.p - (z * lnS.se)
       lnS.u <- lnS.p + (z * lnS.se)
       
       S.p <- exp(lnS.p)
       S.l <- exp(lnS.l)
       S.u <- exp(lnS.u)
       rval <- as.data.frame(cbind(S.p, S.l, S.u))
       names(rval) <- c("est", "lower", "upper")
     }
     
     # Use Skrondal (2003) method if model type is mle2:
     if(class(model)[1] == "mle2"){     
       # Confidence interval for S assuming regression coefficients are from a linear odds model (see appendix of Skrondal, 2003):
       S.p <- (exp(theta3) - 1) / (exp(theta1) + exp(theta2) - 2)
       lnS.p <- log(S.p) 
       
       c <- (1 / (theta1 + theta2 + theta3)) - (1 / (theta1 + theta2))
       d <- 1 / (theta1 + theta2 + theta3)
       
       # Covariance matrix from the model.
       # Diagonals entries are the variances of the regression coefficients.
       # Off-diagonals are the covariance between the corresponding regression coefficients.
       tvcov <- vcov(model)
       
       theta1.var <- tvcov[coeff[1], coeff[1]]
       theta2.var <- tvcov[coeff[2], coeff[2]]
       theta3.var <- tvcov[coeff[3], coeff[3]]
       
       theta12.cov <- tvcov[coeff[1], coeff[2]]
       theta13.cov <- tvcov[coeff[1], coeff[3]]
       theta23.cov <- tvcov[coeff[2], coeff[3]]
       
       lnS.se <- sqrt((c^2 * theta1.var) + (c^2 * theta2.var) + (d^2 * theta3.var) + (2 * c^2 * theta12.cov) + (2 * c * d * theta13.cov) + (2 * c * d * theta23.cov))
       
       lnS.l <- lnS.p - (z * lnS.se)
       lnS.u <- lnS.p + (z * lnS.se)
       
       S.l <- exp(lnS.l)
       S.u <- exp(lnS.u)
       
       rval <- data.frame(S.p, S.l, S.u)
       names(rval) <- c("est", "lower", "upper")
     }
   }
   return(rval)
}

Tboot <-
function(model, significance=0.05, hc=4,
                  double=FALSE, J=NULL, K=NULL, distribution="rademacher"){
  
  if(class(model)!="lm") stop("The argument model must have class lm.")
  if(significance>=1 || significance<=0){
    stop("The significance level should belong to the open interval (0,1).")
  }
  if(hc%in%c(0,2,3,4,5) == FALSE){ 
    warning("The argument hc should be 0, 2, 3, 4 or 5. How did you choose hc = ", hc, " that is not an option,
  the calculation is considering hc = 4.")
    hc=4
  }
  
  # Booth, J.G. and Hall, P. (1994). Monte Carlo approximation and the iterated bootstrap.
  # Biometrika, 81, 331-340.
  if(is.null(J)==TRUE || is.null(K)==TRUE){
    L = length(model$residuals)^3
    gamma2 = ((1/2)*significance^(-2) * (1-significance)*(significance+1/4))^(1/3)
    J = as.integer(gamma2*L^(2/3))
    K = as.integer(gamma2^(-1)*L^(1/3))
  }  

  is.wholenumber <- function(x, tol = .Machine$double.eps^0.5)  abs(x - round(x)) < tol
  
  while(is.wholenumber(K/2)==FALSE && is.wholenumber((J+1)*significance)==FALSE){
    while(is.wholenumber(K/2)==FALSE){
      K = K+1
    }
    while(is.wholenumber((J+1)*significance)==FALSE){
      J = J+1
    }
    while(is.wholenumber((J+1)/K)==FALSE){
      K = K+1
    }  
  }
  number_parameters = length(model$coefficients)
  X = as.matrix(cbind(1,model$model[,-1]))
  n = dim(X)[1]
  beta = as.vector(model$coefficients)
  h = as.vector(hatvalues(model))
  Xbeta = X%*%beta
  error_hat = as.vector(model$residuals)
  root_1_less_h = sqrt(1-h)
  standard_error = as.vector(sqrt(diag(HC(model,method=4))))
  Z = z_star = matrix(,nrow=J,ncol=number_parameters)
  Z_temp = matrix(,nrow=K,ncol=number_parameters)
  
  # Aqui inicia o primeiro bootstrap.
  for(j in 1:J){
    
    if(distribution=="rademacher"){
      t_star = as.vector(sample(c(-1,1),size=length(model$fitted.values),
                                replace=TRUE))
    }else t_star = rnorm(length(model$fitted.values),0,1) 
    
    y_star = as.vector(Xbeta + t_star*error_hat/root_1_less_h)
    model_string = as.character(model$call$formula)
    model_star = lm(formula=as.formula(paste("y_star~",
                    as.character(model_string[3]),sep="")))
    error_hat_star = as.vector(model_star$residuals)
    beta_star = as.vector(model_star$coefficients)
    standard_error_star = sqrt(diag(HC(model_star, method=hc)))
    
    for(m in 1:number_parameters){
      z_star[j,m] = as.numeric((beta_star[m] - beta[m])/standard_error_star[m])
    }
    
    Xbeta_star = X%*%beta_star
    
    if(double==TRUE){
      # Aqui começa o bootstrap duplo
      for(k in 1:K){
        if(distribution=="rademacher"){
          t_star_star = as.vector(sample(c(-1,1),size=length(model$fitted.values),
                                    replace=TRUE))
        }else t_star_star = rnorm(length(model$fitted.values),0,1)
        y_star_star = as.vector(Xbeta_star + t_star_star*error_hat_star/root_1_less_h)
        model_string = as.character(model$call$formula)
        model_star_star = lm(formula=as.formula(paste("y_star_star~",
                                                as.character(model_string[3]),sep="")))
        beta_star_star = as.vector(model_star_star$coefficients)
        standard_error_star_star = sqrt(diag(HC(model_star_star, method=hc)))
      
        for(l in 1:number_parameters){
          z_star_star = as.numeric((beta_star_star[l] - beta_star[l])/standard_error_star_star[l])
          if(z_star_star<=as.numeric(z_star[j,l])) Z_temp[k,l] = 1
          else Z_temp[k,l] = 0
        }
      } # Aqui  termina o bootstrap duplo
      Z[j,] = as.vector(apply(Z_temp,2,mean))
    } # AQUI TERMINA O PRIMEIRO IF DO BOOTSTRAP.
  } # Aqui termina o primeiro bootstrap
  
    ic_inf_simple <- vector()
    ic_sup_simple <- vector()
    ic_inf_double <- vector()
    ic_sup_double <- vector()
  
    for(m in 1:number_parameters){
      ic_inf_simple[m] = beta[m] - as.numeric(quantile(as.vector(z_star[,m]),
                                                       1-significance/2))*standard_error[m]
      ic_sup_simple[m] = beta[m] - as.numeric(quantile(as.vector(z_star[,m]),
                                                       significance/2))*standard_error[m]

      if(double==TRUE){                                               
        ic_inf_double[m] = beta[m] - 
          as.numeric(quantile(as.vector(z_star[,m]),
                   as.numeric(quantile(Z[,m],1-significance/2))))*standard_error[m]
        ic_sup_double[m] = beta[m] - 
          as.numeric(quantile(as.vector(z_star[,m]),
                   as.numeric(quantile(Z[,m],significance/2))))*standard_error[m]
      }
    }
    result = list("beta" = beta,"ci_lower_simple" = ic_inf_simple, "ci_upper_simple" = ic_sup_simple,
                  "ci_lower_double" = ic_inf_double, "ci_upper_double" = ic_sup_double, "J" = J,
                  "K" = K)
    class(result) <- "list"
    return(result)
}

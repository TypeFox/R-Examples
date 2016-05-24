Pboot <- function(model, significance=0.05,
                  double=FALSE, J=NULL, K=NULL, distribution="rademacher"){
  
  if(class(model)!="lm") stop("The argument model must have class lm.")
  if(significance>=1 || significance<=0){
    stop("The significance level should belong to the open interval (0,1).")
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
  U = matrix_beta_star = matrix(,nrow=J,ncol=number_parameters)
  matrix_u = matrix(,nrow=K,ncol=number_parameters)
  
  # Aqui inicia o primeiro bootstrap.
  #browser()
  
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
    matrix_beta_star[j,] = as.vector(beta_star)
    
    Xbeta_star = X%*%beta_star
    
    # Aqui começa o bootstrap duplo.
    if(double==TRUE){
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
        
        for(m in 1:number_parameters){
          if(beta_star_star[m]<=2*beta_star[m]-beta[m]){
            matrix_u[k,m] = 1
          }else matrix_u[k,m] = 0
        }
      } # Aqui termina o segundo bootstrap.
      
      U[j,] = as.vector(apply(matrix_u,2,mean))
    } # Aqui termina o IF 
  } # Aqui termina o primeiro bootstrap.
  
  ic_inf_simple = ic_sup_simple = ic_inf_double = ic_sup_double <- vector()
  
  for(m in 1:number_parameters){  
    ic_inf_simple[m] = quantile(as.vector(matrix_beta_star[,m]), significance/2)
    ic_sup_simple[m] = quantile(as.vector(matrix_beta_star[,m]), 1-significance/2)
    
    if(double==TRUE){
      ic_inf_double[m] = quantile(as.vector(matrix_beta_star[,m]),
                                  quantile(as.vector(U[,m]), significance/2))
      ic_sup_double[m] = quantile(as.vector(matrix_beta_star[,m]),
                                  quantile(as.vector(U[,m]), 1-significance/2))
    }
  }
  
  result = list("beta" = beta,"ci_lower_simple" = ic_inf_simple,
                "ci_upper_simple" = ic_sup_simple, "ci_lower_double" = ic_inf_double,
                "ci_upper_double" = ic_sup_double)
  class(result) <- "list"
  return(result)
}

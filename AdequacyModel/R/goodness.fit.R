goodness.fit <-
function(pdf, cdf, starts, data, method = "L-BFGS-B", domain = c(0,Inf),
                         mle = NULL){
  
  if(missingArg(cdf)==TRUE) stop("Unknown cumulative distribution function. The function needs to be informed.")
  if(missingArg(pdf)==TRUE) stop("Unknown probability density function. The function needs to be informed.")
  if(class(pdf)!="function") stop("The argument pdf must be a function. See the example in the documentation!")
  if(class(cdf)!="function") stop("The argument cdf must be a function. See the example in the documentation!")
  if(missingArg(data)==TRUE) stop("Database missing!")
  if(TRUE%in%is.nan(data)==TRUE) warning("The data have missing information!")
  if(length(domain)!=2) stop("The domain must have two arguments!")
  
  if(is.null(mle)==TRUE){
    if(missingArg(starts)==TRUE) stop("The initial shots were not informed!")
  }else{
    starts = mle 
  }   
  
  # Verifying properties of cumulative distribution function.
  
  if(cdf(par=starts, x = domain[2])!=1) warning("The cdf function informed is not a cumulative distribution function! The function no takes value 1 in Inf.")
  if(cdf(par=starts, x = domain[1])!=0) warning("Check if the cumulative distribution informed is actually a distribution function.")
  
  myintegrate = function(...) tryCatch(integrate(...), error=function(e) NA)
  
  value_int = as.numeric(myintegrate(f=pdf,par=starts,lower=domain[1],
                                     upper=domain[2])[1])
  if(isTRUE(is.na(value_int))==TRUE) warning("Make sure that pdf is a probability density function. The integral in the domain specified is not 					     convergent.")
  
  if(isTRUE(is.na(value_int))!=TRUE){
     #Verifying properties of probability density function.
    if(value_int<0.90){
      warning("The integral from ", domain[1], " to ", domain[2]," of the probability density function has different from 1. Make sure the option 		      domain is correct.")
    }
    
    if(round(value_int)!=1){
      warning("pdf is not a probability density function.")
    }
  }
  
  if(is.null(mle)==TRUE){    
    likelihood = function(par,x){
      -sum(log(pdf(par,x)))
    }
    
    if(method == "Nelder-Mead" || method == "N"){
      result = optim(par = starts, fn = likelihood, x = data,
                        method = "Nelder-Mead", hessian = TRUE)
    }
    
    if(method == "CG" || method == "C"){
      result = optim(par = starts, fn = likelihood, x = data,
                        method = "CG", hessian = TRUE)
    }  
    
    if(method == "L-BFGS-B" || method == "L"){
      result = optim(par=starts, fn = likelihood,method="L-BFGS-B", x = data,
                        lower=c(1e-10,1e-10,1e-10,1e-10,1e-10), upper=c(Inf,Inf,Inf,Inf,Inf), hessian=TRUE)
    }
    
    if(method == "SANN" || method == "S"){
      result = optim(par = starts, fn = likelihood, x = data,
                        method = "SANN", hessian = TRUE)
    }  
    
    if(method == "BFGS" || method == "B"){
      result = optim(par = starts, fn = likelihood, x = data,
                        method = "BFGS", hessian = TRUE)
    }
    
    if((FALSE %in% (method != c("L-BFGS-B", "L", "BFGS", "B",
                                "Nelder-Mead", "N", "SANN", "S", "CG", "C")))==FALSE){
      stop("Valid options are: L-BFGS-B or L, BFGS or B, Nelder-Mead or N, SANN or S, CG or C.")
    }
    
    parameters = result$par
    hessiana = result$hessian # matriz hessiana.
    
    data_orderdenados = sort(data)
    v = cdf(as.vector(parameters), data_orderdenados) # Dados ordenados.
    v[v==1] = 0.99999999999999994
    n = length(data) # Tamanho da amostra.
    y = qnorm(v) # Inversa da acumulada da normal.
    u = pnorm((y-mean(y))/sqrt(var(y)))
    
    W_temp <- vector()
    A_temp <- vector()
    
    for(i in 1:n){
      W_temp[i] = (u[i] - (2*i-1)/(2*n))^2
      A_temp[i] = (2*i-1)*log(u[i]) + (2*n+1-2*i)*log(1-u[i])
    }
    
    A_2 = -n - mean(A_temp)
    W_2 = sum(W_temp) + 1/(12*n)
    W_star = W_2*(1+0.5/n)
    A_star = A_2*(1+0.75/n + 2.25/n^2)
    
    p = length(parameters)
    log.likelihood = -1*likelihood(parameters,data) 
    AICc = -2*log.likelihood + 2*p + 2*(p*(p+1))/(n-p-1)
    AIC  = -2*log.likelihood + 2*p
    BIC  = -2*log.likelihood + p*log(n)
    HQIC = -2*log.likelihood + 2*log(log(n))*p
    ks.testg = function(...) tryCatch(ks.test(...),
                                      warning = function(war) NA)
    KS = ks.test(x = data, y= "cdf", par = as.vector(parameters))
    
    result = (list("W" = W_star,"A" = A_star, "KS" = KS,
                      "mle" = parameters, "AIC" = AIC ,"CAIC " = AICc,
                      "BIC" = BIC, "HQIC" = HQIC, "Erro" = sqrt(diag(solve(hessiana))),
                      "Value" = result$value, "Convergence" = result$convergence))
    class(result) <- "list" 
    return(result)
  }
  
  if(class(mle)=="numeric"){
    
    likelihood = function(par,x){
      -sum(log(pdf(par,x)))
    }
    
    parameters = mle
    data_orderdenados = sort(data)
    v = cdf(as.vector(parameters),data_orderdenados) # Dados ordenados.
    v[v==1] = 0.99999999999999994
    n = length(data) # Tamanho da amostra.
    y = qnorm(v) # Inversa da acumulada da normal.
    u = pnorm((y-mean(y))/sqrt(var(y)))
    
    W_temp <- vector()
    A_temp <- vector()
    
    for(i in 1:n){
      W_temp[i] = (u[i] - (2*i-1)/(2*n))^2
      A_temp[i] = (2*i-1)*log(u[i]) + (2*n+1-2*i)*log(1-u[i])
    }
    
    A_2 = -n - mean(A_temp)
    W_2 = sum(W_temp) + 1/(12*n)
    W_star = W_2*(1+0.5/n)
    A_star = A_2*(1+0.75/n + 2.25/n^2)
    
    p = length(parameters)
    log.likelihood = -1*likelihood(parameters,data) 
    AICc = -2*log.likelihood + 2*p + 2*(p*(p+1))/(n-p-1)
    AIC  = -2*log.likelihood + 2*p
    BIC  = -2*log.likelihood + p*log(n)
    HQIC = -2*log.likelihood + 2*log(log(n))*p
    ks.testg = function(...) tryCatch(ks.test(...),
                                      warning = function(war) NA)
    KS = ks.test(x = data, y= "cdf", par = as.vector(parameters))
    
    result = (list("W" = W_star,"A" = A_star, "KS" = KS, "AIC" = AIC,
                      "CAIC" = AICc, "BIC" = BIC, "HQIC" = HQIC))
    class(result) <- "list" 
    return(result)
  }
}

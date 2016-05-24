function.EL <-
function(data, marker, status, tag.healthy = 0, CFN, CFP, control = control.gsym.point(), confidence.level)
{
  # Marker in the healthy population:
  X0 <- data[data[,status] == tag.healthy, marker]
  # Marker in the diseased population:
  X1 <- data[data[,status] != tag.healthy, marker]

  X0 = sort(X0)
  X1 = sort(X1)

  n0 <- length(X0)
  n1 <- length(X1)

  # Total sample size:
  n = n0 + n1

  EL_c_grid<- rep(NA,length = n)

  tocontrol=list("trace"=0)

  # Constants needed for the empirical likelihood method:
  # c_sampling for resampling,
  # c_F for estimating the distribution,
  # c_ELq for estimating the empirical likelihood function,
  # and c_R for estimating the ROC curve:
  constant1 = control$c_sampling
  constant2 = control$c_F
  constant3 = control$c_ELq
  constant4 = control$c_R

  # Compute the bandwidths b0 and b1 for resampling from healthy and diseased
  # populations, respectively
  b0 = constant1*sd(X0)*n0^(-0.2)
  b1 = constant1*sd(X1)*n1^(-0.2)

  # First resampling for finding the initial t:
  ele = 1
  count = 0

  B <- control$B

  t0hat_b <- rep(NA,length = B)
  s1hat_b <- rep(NA,length = B)

  chat_b <- rep(NA,length = B)
  spehat_b <- rep(NA,length = B)
  senhat_b <- rep(NA,length = B)
 
  rho = CFP/CFN

  lower = 0
  upper = min(c(1/rho,1))
  
  lower2 = 1-min(c(rho,1))
  upper2 = 1

  while (ele <= B)
  {

    X0b = resampling_gauss(b0, n0, X0)
    X1b = resampling_gauss(b1, n1, X1)

    X0b = sort(X0b)
    X1b = sort(X1b)

    wb = constant2*sd(X0b)*n0^(-1/3)
    x0b = relative_sample(X0b,X1b,wb)

    w1b = constant2*sd(X1b)*n1^(-1/3)
    x1b = relative_sample(X1b,X0b,w1b)

    w0b = constant4*sd(x0b)*n0^(-1/3)
    if (w0b > 0)
    {
	parameters = c(rho,w0b,x0b)
      f_lower = cOpt_gsym_kernel(lower, parameters)
  	f_upper = cOpt_gsym_kernel(upper, parameters)
	value1 = f_lower*f_upper
    }
    else
    {
      parameters = c(rho,x0b)
      f_lower = cOpt_gsym_empirical(lower, parameters)
      f_upper = cOpt_gsym_empirical(upper, parameters)
      value1 = f_lower*f_upper
    }
    w11b = constant4*sd(x1b)*n1^(-1/3)
    if (w11b > 0)
    {
	parameters = c(rho,w11b,x1b)
      f_lower = cOpt_gsym2_kernel(lower2, parameters)
	f_upper = cOpt_gsym2_kernel(upper2, parameters)
	value2 = f_lower*f_upper
    }
    else
    {
	parameters = c(rho,x1b)
        f_lower = cOpt_gsym2_empirical(lower2, parameters)
        f_upper = cOpt_gsym2_empirical(upper2, parameters)
        value2 = f_lower*f_upper
    }        

    # If the samples do overlap:    
    if (max(X0b)>min(X1b) & value1 < 0 & value2 < 0)
    { 
	   w0b = constant4*sd(x0b)*n0^(-1/3)
    	   if (w0b > 0)
    	   {
   	  	parameters = c(rho,w0b,x0b)
                t0hat_b[ele] <- uniroot(cOpt_gsym_kernel, interval = c(0,min(c(1/rho,1))), parameters = parameters, tol = 1e-12)$root          
           }
           # If the bandwith is equal to zero, we use the empirical:  
    	   else
           {
      	  	parameters = c(rho,x0b)
                t0hat_b[ele] <- uniroot(cOpt_gsym_empirical,interval = c(0,min(c(1/rho,1))), parameters = parameters,tol = 1e-12)$root
           }

           w11b = constant4*sd(x1b)*n1^(-1/3)
          
           if (w11b > 0)
           {
   	  	parameters = c(rho,w11b,x1b)
                s1hat_b[ele] <- uniroot(cOpt_gsym2_kernel, interval =  c(1-min(c(rho,1)),1), parameters = parameters, tol = 1e-12)$root    
           }
    
           # If the bandwith is equal to zero, we use the empirical:  
	   else
           {
          	parameters = c(rho,x1b)
                s1hat_b[ele] <- uniroot(cOpt_gsym2_empirical, interval = c(1-min(c(rho,1)),1), parameters = parameters, tol = 1e-12)$root
    	   }

           spehat_b[ele] = 1-(0.5*t0hat_b[ele]+0.5*((1-s1hat_b[ele])/rho))     
           senhat_b[ele] = 1-rho*(0.5*t0hat_b[ele]+0.5*((1-s1hat_b[ele])/rho))
            
           bb0 = constant3*sd(X0b)*n0^(-0.5)
           bb1 = constant3*sd(X1b)*n1^(-0.5)

           c_grid = sort(c(X0b,X1b))
           for (eme in 1:n)
           {
          	EL_c_grid[eme] = ELquantiles_gsym(c_grid[eme], spehat_b[ele], X0b, X1b, n0, n1, bb0, bb1, rho)          
           }

           ind_c_ini = which(EL_c_grid == min(EL_c_grid))
           c_ini = mean(c_grid[ind_c_ini])
     
           chat_b[ele] <- solnp (pars = c_ini, fun = ELquantiles_gsym, control=tocontrol, q0 = spehat_b[ele], X0 = X0b, X1 = X1b, n0 = n0, n1 = n1, h0 = bb0, h1 = bb1, rho = rho)$pars
           ele = ele+1
     }

    
     # If the samples do not overlap:
     else
     {
     	count = count+1
     }
  }

  # Point estimates of cutpoint c, spe and sen:
  w = constant2*sd(X0)*n0^(-1/3)
  ww = constant2*sd(X1)*n1^(-1/3)

  x0 = relative_sample(X0,X1,w)
  x1 = relative_sample(X1,X0,ww)

  bb0 = constant3*sd(X0)*n0^(-0.5)
  bb1 = constant3*sd(X1)*n1^(-0.5)

  w0 = constant4*sd(x0)*n0^(-1/3)
  w1 = constant4*sd(x1)*n1^(-1/3)
  
  if (w0 > 0)
  {
      parameters = c(rho,w0,x0)      
      t0tilde <- uniroot (cOpt_gsym_kernel, interval = c(0,min(c(1/rho,1))), tol = 1e-12, parameters = parameters)$root      
  }
  
  # If the bandwidth is equal to zero, we use the empirical:
  else
  {
      parameters = c(rho,x0)      
      t0tilde <- uniroot (cOpt_gsym_empirical, interval = c(0,min(c(1/rho,1))), tol = 1e-12, parameters = parameters)$root      
  }

  
  if (w1 > 0)
  {
      parameters = c(rho,w1,x1)
      s1tilde <- uniroot (cOpt_gsym2_kernel, interval = c(1-min(c(rho,1)),1), tol = 1e-12, parameters = parameters)$root      
  }

  # If the bandwdidth is equal to zero, we use the empirical:
  else
  {
      parameters = c(rho,x1)
      s1tilde <- uniroot (cOpt_gsym2_empirical, interval = c(1-min(c(rho,1)),1), tol = 1e-12, parameters = parameters)$root      
      
  }

  spehat = 1-(0.5*t0tilde+0.5*((1-s1tilde)/rho))  
  senhat = 1-rho*(0.5*t0tilde+0.5*((1-s1tilde)/rho))  

  c_grid = sort(c(X0,X1))
  for (eme in 1:n)
  {
      EL_c_grid[eme] = ELquantiles_gsym(c_grid[eme], spehat, X0, X1, n0, n1, bb0, bb1, rho)      
  }

  ind_c_ini = which(EL_c_grid == min(EL_c_grid))
  c_ini = mean(c_grid[ind_c_ini])

  chat <- solnp (pars = c_ini, fun = ELquantiles_gsym, control = tocontrol, q0 = spehat, X0 = X0, X1 = X1, n0 = n0, n1 = n1, h0 = bb0, h1 = bb1, rho = rho)$pars  

  CIEL_c <- numeric(length = 2)
  CIEL_spe <- matrix(NA, nrow = 1, ncol = 2)
  CIEL_sen <- matrix(NA, nrow = 1, ncol = 2)

  chat_b = sort(chat_b)
  CIEL_c[1] = chat_b[(B+1)*(1-confidence.level)/2]
  CIEL_c[2] = chat_b[(B+1)*(1-(1-confidence.level)/2)]
  
  spehat_b = sort(spehat_b)
  CIEL_spe[1,1:2] = c(spehat_b[(B+1)*(1-confidence.level)/2],spehat_b[(B+1)*(1-(1-confidence.level)/2)])
  
  senhat_b = sort(senhat_b)
  CIEL_sen[1,1:2] = c(senhat_b[(B+1)*(1-confidence.level)/2],senhat_b[(B+1)*(1-(1-confidence.level)/2)])
  
  optimal.cutoff <- matrix(ncol=3, nrow = length(chat))

  optimal.cutoff [,1] <- chat
  optimal.cutoff [,-1] <- CIEL_c

  Sp <- matrix(ncol=3, nrow = length(chat))
  Sp[,1] <- spehat
  Sp[,-1] <- CIEL_spe

  Se <- matrix(ncol=3, nrow = length(chat))
  Se[,1] <- senhat
  Se[,-1] <- CIEL_sen

  optimal.result <- list(cutoff = optimal.cutoff, Specificity = Sp, Sensitivity = Se)
  AUC <- calculate.empirical.AUC(data, marker, status, tag.healthy, confidence.level)

  GPQ <- function.GPQ(data, marker, status, tag.healthy = 0, CFN, CFP, control = control.gsym.point(), confidence.level)

  # If original data are not normally distributed:
  if ("lambda" %in% names(GPQ))     
  { 
  	res <- list (optimal.result = optimal.result, AUC = AUC, rho = rho, lambda = GPQ$lambda, normality.transformed= GPQ$normality.transformed, pvalue.healthy = GPQ$pvalue.healthy, pvalue.diseased = GPQ$pvalue.diseased, pvalue.healthy.transformed = GPQ$pvalue.healthy.transformed, pvalue.diseased.transformed = GPQ$pvalue.diseased.transformed)  
  } 

  # If original data are normally distributed:
  else
  {
  	res <- list (optimal.result = optimal.result, AUC = AUC, rho = rho, pvalue.healthy = GPQ$pvalue.healthy, pvalue.diseased = GPQ$pvalue.diseased)  
  }  

  return (res)
}

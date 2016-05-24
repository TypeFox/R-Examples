function.GPQ <-
function(data, marker, status, tag.healthy = 0, CFN, CFP, control = control.gsym.point(), confidence.level)

{
  # Marker in the healthy population:
  X0 <- data[data[,status] == tag.healthy, marker]
  # Marker in the diseased population:
  X1 <- data[data[,status] != tag.healthy, marker]

  n0 <- length(X0)
  n1 <- length(X1)


  model ="norm"

  # flag =  0 normal assumption for original data
  # flag =  1 a BoxCox transformation is done with lambda estimated via ML
  # flag = -1 a BoxCox transformation is done with lambda estimated via ML,
  #           besides the sign of the transformed biomarkers is changed
  #           due to the violation of the assumed criterium that
  #           larger values of the biomarker are associated with the disease

  # normality.transformed = "yes" normal assumption for Box-Cox transformed data 

  pvalue0 <- shapiro.test(X0)$p.value  # p-value for normality assumption healthy data
  pvalue1 <- shapiro.test(X1)$p.value  # p-value for normality assumption diseased data     

  # If original data are not normally distributed:
  if((pvalue0 < 0.05) | (pvalue1 < 0.05))
  {

	lambda_sol = BoxCox_binormal_MLestimate(X0,X1,n0,n1)

        if (abs(lambda_sol) < 10^(-15))
  	{
      		X0 = log(X0)
      		X1 = log(X1)
       	}

  	else
  	{
      		X0 =((X0^(lambda_sol))-1)/lambda_sol
      		X1 =((X1^(lambda_sol))-1)/lambda_sol
  	}

  	if (mean(X0) > mean(X1))
  	{
    		flag = -1
    		X0 = -X0
    		X1 = -X1
  	}

  	else
  	{
    		flag = 1
  	}

   }

   # Original data are normally distributed:
   else
   {
  	flag = 0	    
   }


   X0 = sort(X0)
   X1 = sort(X1)

   # If original data are not normally distributed:
   pvalue0.transformed <- shapiro.test(X0)$p.value
   pvalue1.transformed <- shapiro.test(X1)$p.value

   if(flag != 0)
   { 
	if((pvalue0.transformed>=0.05) & (pvalue1.transformed>=0.05))
  	{
        	normality.transformed = "yes"
  	}
  	
	else 
	{
		normality.transformed = "no"
        }
   }


   m0 = mean(X0)
   m1 = mean(X1)

   s0 = sd(X0)
   s1 = sd(X1)
   a = (m1-m0)/s1
   b = s0/s1
  
  rho = CFP/CFN
  
  parameters = c(a,b,rho)


  if (rho == 1)   # the symmetry point
  {

  	# Exact solution:
  	that <- 1-pnorm(a/(1+b))  
  	chat = (s1*m0+s0*m1)/(s0+s1)

  }

  else
  {
	that <- uniroot(f = cOpt_gsym, interval =  c(0,min(c(1/rho,1))), tol = 1e-12, parameters = parameters, model = model)$root
        chat = qnorm(1-that, m0, s0)
  }

  spehat = 1-that
  senhat = 1-rho*that

  t <- numeric(length = control$I)
  c <- numeric(length = control$I)
  spe <- numeric(length = control$I)
  sen <- numeric(length = control$I)

  for (k in 1:control$I)
  {
    t0 = rt(1, df = n0-1)
    t1 = rt(1, df = n1-1)
    chi2_0 = rchisq(1, df = n0-1)
    chi2_1 = rchisq(1, df = n1-1)
    Rm0 = m0-(t0*s0/sqrt(n0))
    Rm1 = m1-(t1*s1/sqrt(n1))
    Rs0 = sqrt((n0-1)*(s0^2)/chi2_0)
    Rs1 = sqrt((n1-1)*(s1^2)/chi2_1)
    Ra = (Rm1-Rm0)/Rs1
    Rb = Rs0/Rs1
    parameters = c(Ra,Rb,rho)

    if (rho == 1)   # the symmetry point
    {
	p<-Ra/(1+Rb)    
    	t[k] <- 1-pnorm(p)

	# The optimal cutpoint is computed:
    	c[k] = (Rs1*Rm0+Rs0*Rm1)/(Rs0+Rs1)

    }

    else
    {
	t[k] <- uniroot(f = cOpt_gsym, interval = c(0,min(c(1/rho,1))), tol = 1e-12, parameters = parameters, model = model)$root

    	# The optimal cutpoint is computed:
    	c[k] = qnorm(1-t[k], mean = Rm0, sd = Rs0)
    }
    

    # The optimal Specificity is computed:
    spe[k] = 1-t[k]

    # The optimal Sensitivity is computed:
    sen[k] = 1-rho*t[k]
  }

  c = sort(c)
  CIgp_c <- numeric(length = 2)
  CIgp_c[1] = c[ceiling(0.5*(1-confidence.level)*control$I)]
  CIgp_c[2] = c[ceiling((1-0.5*(1-confidence.level))*control$I)]
  
  sen = sort(sen)
  CIgp_sen <- numeric(length = 2)
  CIgp_sen[1] = sen[ceiling(0.5*(1-confidence.level)*control$I)]
  CIgp_sen[2] = sen[ceiling((1-0.5*(1-confidence.level))*control$I)]
  
  spe = sort(spe)
  CIgp_spe <- numeric(length = 2)
  CIgp_spe[1] = spe[ceiling(0.5*(1-confidence.level)*control$I)]
  CIgp_spe[2] = spe[ceiling((1-0.5*(1-confidence.level))*control$I)]
  

  if (flag == -1)
  {
      CIgp_c = (-CIgp_c)
      chat = (-chat)
      if (abs(lambda_sol) < 10^(-15))
      {
          CIgp_c = sort(exp(CIgp_c))
          chat = exp(chat)
      }
      else
      {
          CIgp_c = sort((lambda_sol*CIgp_c+1)^(1/lambda_sol))
          chat = (lambda_sol*chat+1)^(1/lambda_sol)
      }
      
   }

   if (flag == 1)
   {
      if (abs(lambda_sol) < 10^(-15))
      {
          CIgp_c = sort(exp(CIgp_c))
          chat = exp(chat)
      }
      else
      {
        CIgp_c = (lambda_sol*CIgp_c+1)^(1/lambda_sol)
        chat = (lambda_sol*chat+1)^(1/lambda_sol)
      }
      
    }

   c.names <- c("Value", "ll", "ul")

   optimal.cutoff <- matrix(ncol=3, nrow = length(chat), dimnames = list(1:length(chat), c.names))
   
   optimal.cutoff [,1] <- chat
   optimal.cutoff [,-1] <- CIgp_c

   Sp <- matrix(ncol=3, nrow = length(chat), dimnames = list(1:length(chat), c.names))

   Sp[,1] <- spehat
   Sp[,-1] <- CIgp_spe

   Se <- matrix(ncol=3, nrow = length(chat), dimnames = list(1:length(chat), c.names))
   Se[,1] <- senhat
   Se[,-1] <- CIgp_sen

   optimal.result <- list(cutoff = optimal.cutoff, Specificity = Sp, Sensitivity = Se)
   AUC <- calculate.empirical.AUC(data, marker, status, tag.healthy, confidence.level)

   # If original data are normally distributed:
   if (flag == 0)
   {
   	res <- list (optimal.result = optimal.result, AUC = AUC, rho = rho, pvalue.healthy = pvalue0, pvalue.diseased = pvalue1)
   }

   # If original data are not normally distributed:
   else
   {
   	res <- list (optimal.result = optimal.result, AUC = AUC, rho = rho, lambda = lambda_sol, normality.transformed = normality.transformed, pvalue.healthy = pvalue0, pvalue.diseased = pvalue1, pvalue.healthy.transformed = pvalue0.transformed, pvalue.diseased.transformed = pvalue1.transformed)
   }
   return(res)

}

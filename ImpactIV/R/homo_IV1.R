homo_IV1 <-
function(Z, A, M, Y, X)
{
           ##Sample size
           N = length(Z)


           ###############################Step 1
           ##estimate p
           p_hat = mean(Z)


           ###############################Step 2
           ##Fit Logit model
           ##Transform the data: (A,M)->0,1,2,3
           Logit = 1*(A==0&M==0) + 2*(A==1&M==0) + 3*(A==0&M==1) + 4*(A==1&M==1)
           Logit = factor(Logit)

           ##prepare the covariates
           ##X contains a colum of 1s
           X = cbind(rep(1, N), X)
           ##interactions between X and Z
           X_data = data.frame(cbind(X, X*Z))

           ##Logit model
           logit_fit = multinom(Logit ~ . + 0, data = X_data)
           ##the coefficients
           gamma = coef(logit_fit)

           ##P^(a,m)|z: z=0 and z=1
           Z1 = rep(1, N)
           Z0 = rep(0, N)
           X1 = cbind(X, X*Z1)
           X0 = cbind(X, X*Z0)
           ##fitted values
           P10_1 = X1%*%gamma[1,]
           P10_0 = X0%*%gamma[1,]
           P01_1 = X1%*%gamma[2,]
           P01_0 = X0%*%gamma[2,]
           P11_1 = X1%*%gamma[3,]
           P11_0 = X0%*%gamma[3,]
           
           ##probability
           total_1 = 1 + exp(P10_1) + exp(P01_1) + exp(P11_1)
           total_0 = 1 + exp(P10_0) + exp(P01_0) + exp(P11_0)
           P10_1 = exp(P10_1)/total_1
           P10_0 = exp(P10_0)/total_0
           P01_1 = exp(P01_1)/total_1
           P01_0 = exp(P01_0)/total_0
           P11_1 = exp(P11_1)/total_1
           P11_0 = exp(P11_0)/total_0
           

           #################################Step 3
           ##instrumental variables
           Cons = rep(1, N)
           EA1_0 = (P10_1 + P11_1) - (P10_0 + P11_0)           
           EM1_0 = (P01_1 + P11_1) - (P01_0 + P11_0)
           EAM1_0 = P11_1 - P11_0
           ##N * 4 matrix
           G = cbind(Cons, EA1_0, EM1_0, EAM1_0)


           #################################Step 4
           ##Estimatinng Equation
           ##W matrix: N * 4
           W = cbind(Z, A, M, A*M)
           WZ_p = W*(Z - p_hat)
           YZ_p = Y*(Z - p_hat)
           ##beta_IV1
           beta = solve(t(G)%*%WZ_p)%*%t(G)%*%YZ_p



           #################################Estimation of the variance
           ##Omega
           g = Y - W%*%beta
           Omega_hat = mean(g^2)
           ##Asymptotic Variance-Covariance matrix      
           COV = solve(t(G)%*%G)*Omega_hat/p_hat/(1-p_hat)
           ##variance
           var = diag(COV)
           ##standard error
           se = sqrt(var)
           ##z value
           zvalue = beta/se
           ##p value
           pvalue = pnorm(abs(zvalue), mean = 0, sd = 1, lower.tail = FALSE)*2
           ##confidence intervals
           CI = cbind(beta + qnorm(0.05/2)*se, beta - qnorm(0.05/2)*se)

           ########A robust version of variance estimator
           A = solve(t(G)%*%WZ_p)
           B = t(G)%*%G
           COVr = A%*%B%*%t(A)*Omega_hat*p_hat*(1-p_hat) 
           ##variance
           varr = diag(COVr)
           ##standard error
           ser = sqrt(varr)
           ##z value
           zvaluer = beta/ser
           ##p value
           pvaluer = pnorm(abs(zvaluer), mean = 0, sd = 1, lower.tail = FALSE)*2
           ##confidence intervals
           CIr = cbind(beta + qnorm(0.05/2)*ser, beta - qnorm(0.05/2)*ser)

           
           #################################Output
           output = list(beta = beta, phat = p_hat, residual = g,
                         se = se, zvalue = zvalue, pvalue = pvalue, 
                         CI = CI,  COV = COV,
                         ser = ser, zvaluer = zvaluer, pvaluer = pvaluer, 
                         CIr = CIr,  COVr = COVr, 
                         samplesize = N, G = G, W = W, Omega = Omega_hat)
           return(output)

}


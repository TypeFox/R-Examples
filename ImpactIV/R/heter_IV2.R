heter_IV2 <-
function(Z, A, M, Y, X, polydegree = 2, step1 = NULL, truncate = 0.25, select = NULL)
{
         
           ###################################Step 1            
           if(is.null(step1))
           {
                 step1 = homo_IV1(Z, A, M, Y, X)  
           }
           ##residuals
           residual = step1$residual              
           ##sample size           
           N = step1$samplesize
           ##p hat
           p_hat = step1$phat


           ###################################Step 2
           ##square of the residuals
           res2 = residual^2
           ##using polynomial regression
           if(polydegree > 1)
           {
                  polynomials = poly(X, degree = polydegree, raw = TRUE)
           }
           if(polydegree == 1)
           {
                  polynomials = X
           }
           polynomials = data.frame(polynomials)
           Omegalm=lm(res2 ~ ., data = polynomials)
           ##choose the best model
           if(is.null(select) == 0) 
           {          
                if(select == "AIC")
                {
                        Omegalm = step(Omegalm, trace = 0)
                }
                if(select == "BIC")
                {
                        Omegalm = step(Omegalm, trace = 0, k = log(N))
                }
           }
           ##the fitted value
           Omega = Omegalm$fitted
           ##Truncated at 0
           Omega = pmax(Omega, truncate)
           ##But the estimated Omega can be close to 0!
           Omega_inv = 1/Omega


           #################################Step 3
           ##instrumental variables
           G = step1$G
           W = step1$W


           #################################Step 4
           ##Estimatinng Equation
           W_Omega_Z_p = W*(Z - p_hat)*Omega_inv
           Y_Omega_Z_p = Y*(Z - p_hat)*Omega_inv
           ##beta_IV2
           beta = solve(t(G)%*%W_Omega_Z_p)%*%t(G)%*%Y_Omega_Z_p


           #################################Variance Estimation
           ##Weighted G
           ##But the fitted value of Omega can be negative!
           G_weight = G*sqrt(Omega_inv)
           ##Asymptotic Variance-Covariance matrix      
           COV = solve(t(G_weight)%*%G_weight)/p_hat/(1-p_hat)
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

           #######A Robust version of variance estimation
           A = solve(t(G)%*%W_Omega_Z_p)
           B = t(G_weight)%*%G_weight
           COVr = A%*%B%*%t(A)*p_hat*(1 - p_hat)
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
           output = list(beta = beta, phat = p_hat, residual = Y - W%*%beta,
                         se = se, zvalue = zvalue, pvalue = pvalue, 
                         CI = CI,  COV = COV,
                         ser = ser, zvaluer = zvaluer, pvaluer = pvaluer, 
                         CIr = CIr,  COVr = COVr, 
                         samplesize = N, G = G, W = W, Omega = Omega)
           return(output)

}


`fplsrPI` <- function(X, Y, Xtest, order, method = c("delta", "boota"), alpha = 0.05, B = 1000, weight, beta = 0.1, 
                      adjust = c(FALSE, TRUE), backh = 10)
{
  p = dim(Y)[1]
  n = ncol(Y)
  method = match.arg(method)
  if (method == "delta"){
      output = plsrPIs(X, Y, Xtest, ncomp1 = 1, ncomp = order, alpha, weight = weight, beta = beta)
  }
  else {
     if (adjust == FALSE){
         output = plsrPIsboot(X, Y, Xtest, ncomp = order, B, alpha, weight = weight, beta = beta)
     }
     if (adjust == TRUE){
           q = matrix(, p, backh)
           p = dim(Y)[1]
           if (weight == FALSE) {
               for (i in 1:backh) {
                    j = (n - backh - 1) + i
                    pred = t(Y[, 1:(j - 1)])
                    resp = t(Y[, 2:j])
                    mresp = apply(resp, 2, mean)
                    mod = plsr(pred, resp, Xtest = Y[, j], ncomp = order, 
                               type = "simpls")
                    q[, i] = mod$Ypred[, , ] 
               }
           }
           if (weight == TRUE){
               for (i in 1:backh){
                    j = (n - backh - 1) + i
                    pred = t(Y[, 1:(j - 1)])
                    resp = t(Y[, 2:j])
                    mresp = apply(resp, 2, mean)
                    mod = plsr(pred, resp, Xtest = Y[, j], ncomp = order, 
                          type = "simpls", weight = TRUE, beta = beta)
                    q[, i] = mod$Ypred[, , ]
               }
           } 
           up = up1 = down = down1 = matrix(, p, 1)
           for (i in 1:p) {
                up[i, ] = quantile(q[i, ], 1 - alpha/2, na.rm = TRUE)
                down[i, ] = quantile(q[i, ], alpha/2, na.rm = TRUE)
                up1[i, ] = quantile(Y[i, (n - backh):n], 1 - alpha/2, 
                           na.rm = TRUE)
                down1[i, ] = quantile(Y[i, (n - backh):n], alpha/2, 
                           na.rm = TRUE)
           }
           lbadj = down1/down
           ubadj = up1/up 
           oldoutput = plsrPIsboot(X, Y, Xtest, ncomp = order, B, alpha, weight = weight, beta = beta)
           lbnew = oldoutput$lb * lbadj
           ubnew = oldoutput$ub * ubadj
           output = list(lb = oldoutput$lb, ub = oldoutput$ub, lbadj = lbnew, ubadj = ubnew, 
                         lbadjfactor = lbadj, ubadjfactor = ubadj)
     }
  }
  return(output)
}

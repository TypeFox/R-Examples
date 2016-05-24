`plsrPIs` <- function(X, Y, Xnew, ncomp1 = 1, ncomp, alpha = .05, weight = c("FALSE", "TRUE"), beta)
{
      n = dim(X)[2]
      p = dim(X)[1]
      if (weight == FALSE){
          pred = t(X)
          mpred = apply(pred,2,"mean")
          decenpred = sweep(pred,2,mpred)
      }
      else {
           q = matrix(,n,1)
           for(i in 1:n){
               q[i,] = beta*(1-beta)^(i-1)
           }
           wt = diag(rev(q))
           pred = t(X)
           resp = t(Y)
           mpred2 = apply(pred,2,mean)
           mresp2 = apply(resp,2,mean)
           decenpred = sweep(pred,2,mpred2)
           decenresp = sweep(resp,2,mresp2)
           mpred1 = t(matrix(rep(mpred2,n),p,n))
           mresp1 = t(matrix(rep(mresp2,n),p,n))           
           pred = wt%*%decenpred + mpred1
           mpred = apply(pred,2,mean)
           decenpred = sweep(pred,2,mpred)
           Y = t(wt%*%decenresp + mresp1)
      }
      resi = matrix(NA,n,p)
      decenXnew = Xnew-mpred
      a = b = se = ub = lb = matrix(NA,p,1)
      imat = diag(1:n)
      for(i in 1:p){
          respv = Y[i,]
          mresp = mean(respv)
          decenresp = respv-mresp
          mod = plsr(pred, respv, Xtest = pred, ncomp = ncomp, type = "simpls",
                     unit.weights = FALSE)
          resi[,i] = mod$Ypred[,,] - respv
          Jaco = Jacob(decenpred, decenresp, ncomp1)$Jacomatrix
          b[i,] = sqrt((t(decenXnew) %*% Jaco %*% t(Jaco) %*% decenXnew + (n+1)/n))
          a[i,] = sqrt((t(resi[,i]) %*% resi[,i]) / (n-p-1))
          se[i,] = a[i,] * b[i,] * qnorm(1 - alpha/2)
      }
      mod2 = plsr(pred, t(Y), Xtest = t(Xnew), ncomp = ncomp, type = "simpls",
                  unit.weights = FALSE)
      predictvalue = matrix(mod2$Ypred[,,], p, 1)
      for(i in 1:p){
          lb[i,] = predictvalue[i,] - se[i,]
          ub[i,] = predictvalue[i,] + se[i,]
      }
      return(list(lb = lb,ub = ub))
}

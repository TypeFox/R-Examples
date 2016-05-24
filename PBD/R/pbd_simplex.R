pbd_simplex = function(trparsopt,idparsopt,trparsfix,idparsfix,exteq,parsfunc,pars2,brts,missnumspec)
{
numpar = length(trparsopt)
## Setting up initial simplex
v = t(matrix(rep(trparsopt,each = numpar + 1),nrow = numpar + 1))
for(i in 1:numpar)
{
    parsoptff = 1.05 * trparsopt[i]/(1 - trparsopt[i])
    trparsoptff = parsoptff/(1 + parsoptff)
    fac = trparsoptff/trparsopt[i]
    if(v[i,i + 1] == 0)
    {
       v[i,i + 1] = 0.00025
    } else {
       v[i,i + 1] = v[i,i + 1] * min(1.05,fac)
    }
}

fv = rep(0,numpar + 1)
for(i in 1:(numpar + 1))
{
   fv[i] = -pbd_loglik_choosepar(v[,i],trparsfix,idparsopt,idparsfix,exteq,parsfunc,pars2,brts,missnumspec)
}

how = "initial"
itercount = 1
string = itercount
for(i in 1:numpar)
{
   string = paste(string, v[i,1]/(1 - v[i,1]), sep = " ")
}
string = paste(string, -fv[1], how, "\n", sep = " ")
cat(string)
flush.console()

tmp = order(fv)
if(numpar == 1)
{
   v = matrix(v[tmp],nrow = 1,ncol = 2)
} else {
   v = v[,tmp]
}
fv = fv[tmp]

## Iterate until stopping criterion is reached
reltolx = as.numeric(pars2[8])
reltolf = as.numeric(pars2[9])
abstolx = as.numeric(pars2[10])
maxiter = as.numeric(pars2[11])
rh = 1
ch = 2
ps = 0.5
si = 0.5

v2 = t(matrix(rep(v[,1],each = numpar + 1),nrow = numpar + 1))

while(itercount <= maxiter & ( ( is.nan(max(abs(fv - fv[1]))) | (max(abs(fv - fv[1])) - reltolf * abs(fv[1]) > 0) ) + ( (max(abs(v - v2)) - reltolx * abs(v2) > 0) | (max(abs(v - v2)) - abstolx > 0) ) ) )
{ 
   ## Calculate reflection point

   if(numpar == 1)
   {
       xbar = v[1]
   } else {
       xbar = rowSums(v[,1:numpar])/numpar
   }
   xr = (1 + rh) * xbar - rh * v[,numpar + 1]
   fxr = -pbd_loglik_choosepar(xr,trparsfix,idparsopt,idparsfix,exteq,parsfunc,pars2,brts,missnumspec)
 
   if(fxr < fv[1])
   {
       ## Calculate expansion point
       xe = (1 + rh * ch) * xbar - rh * ch * v[,numpar + 1]
       fxe = -pbd_loglik_choosepar(xe,trparsfix,idparsopt,idparsfix,exteq,parsfunc,pars2,brts,missnumspec)
       if(fxe < fxr)
       {
           v[,numpar + 1] = xe
           fv[numpar + 1] = fxe
           how = "expand"
       } else {
           v[,numpar + 1] = xr
           fv[numpar + 1] = fxr
           how = "reflect"
       }
   } else {
       if(fxr < fv[numpar])
       {      
           v[,numpar + 1] = xr
           fv[numpar + 1] = fxr
           how = "reflect"
       } else {
           if(fxr < fv[numpar + 1])
           {
              ## Calculate outside contraction point
              xco = (1 + ps * rh) * xbar - ps * rh * v[,numpar + 1]
              fxco = -pbd_loglik_choosepar(xco,trparsfix,idparsopt,idparsfix,exteq,parsfunc,pars2,brts,missnumspec)
              if(fxco <= fxr)
              {
                 v[,numpar + 1] = xco
                 fv[numpar + 1] = fxco            
                 how = "contract outside"
              } else {
                 how = "shrink"
              }
           } else {
              ## Calculate inside contraction point
              xci = (1 - ps) * xbar + ps * v[,numpar + 1]
              fxci = -pbd_loglik_choosepar(xci,trparsfix,idparsopt,idparsfix,exteq,parsfunc,pars2,brts,missnumspec)
              if(fxci < fv[numpar + 1])
              {  
                 v[,numpar + 1] = xci
                 fv[numpar + 1] = fxci
                 how = "contract inside"
              } else {
                 how = "shrink"
              }
           }
           if(how == "shrink")
           {
               for(j in 2:(numpar + 1))
               {

                   v[,j] = v[,1] + si * (v[,j] - v[,1])
                   fv[j] = -pbd_loglik_choosepar(v[,j],trparsfix,idparsopt,idparsfix,exteq,parsfunc,pars2,brts,missnumspec)
               }
           }
       }
   }
   tmp = order(fv)
   if(numpar == 1)
   {
      v = matrix(v[tmp],nrow = 1,ncol = 2)
   } else {
      v = v[,tmp]
   }
   fv = fv[tmp]
   itercount = itercount + 1
   string = itercount;
   for(i in 1:numpar)
   {
       string = paste(string, v[i,1]/(1 - v[i,1]), sep = " ")
   }
   string = paste(string, -fv[1], how, "\n", sep = " ")
   cat(string)
   flush.console()
   v2 = t(matrix(rep(v[,1],each = numpar + 1),nrow = numpar + 1))
}
if(itercount < maxiter)
{
   cat("Optimization has terminated successfully.","\n")
} else {
   cat("Maximum number of iterations has been exceeded.","\n")
}
out = list(par = v[,1], fvalues = -fv[1], conv = as.numeric(itercount > maxiter))
invisible(out)
}
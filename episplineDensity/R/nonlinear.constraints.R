nonlinear.constraints <- function(epiparameters,softinfo,r)
{
# function [g,geq,grad_g,grad_geq] = cons_fun_file(epiparameters,softinfo,optim,r)
# make nonlinear constraints 

N <- epiparameters$Ndiscr
m0 <- epiparameters$m0
mN <- epiparameters$mN
p <- epiparameters$order
Delta <- (mN-m0)/N

w <- c(0.152753387, 0.152753387, 0.149172986, 0.149172986, 0.142096109,
      0.142096109, 0.131688638, 0.131688638, 0.118194532, 0.118194532,
      0.10193012,  0.10193012,  0.083276742, 0.083276742, 0.062672048,
      0.062672048, 0.04060143, 0.04060143, 0.017614007, 0.017614007)

xivec <- c(-0.076526521, 0.076526521, -0.227785851, 0.227785851, -0.373706089,
           0.373706089, 0.510867002, 0.510867002, -0.636053681, 0.636053681,
          -0.746331906, 0.746331906, -0.839116972, 0.839116972, -0.912234428,
           0.912234428, -0.963971927, 0.963971927, -0.993128599, 0.993128599)

geq <- 0; g1 <- 0; g2 <- 0
grad_g1 <- grad_g2 <- grad_geq <- rep (0, (p+2)* N + 1)

if (any (names (softinfo) == "upperbound1moment")) {
    for (k in 1:N) {
        dex <- (Delta/2)*xivec + ((k-1)*Delta + m0 + k*Delta + m0)/2
        g1 <- g1 + (Delta/2) * sum(w * 
               epispline (epiparameters, r, dex, exponentiate=TRUE, moment=1))
# gradient
        grad_g1 <- grad_g1 + (Delta/2) * 
               grad_expepispline (epiparameters, r, dex, moment=1) %*% w
    }
    g1 <- g1 - softinfo$upperbound1moment
      
}
    else g1 <- NULL  

if  (any (names (softinfo) == "upperbound2moment")) {
  
    for (k in 1:N) {
        dex <- (Delta/2)*xivec + ((k-1)*Delta + m0 + k*Delta + m0)/2
        g2 <- g2 + (Delta/2)*sum(w * 
              epispline(epiparameters,r, dex, exponentiate=TRUE, moment=2))
# gradient
        grad_g2 = grad_g2 + (Delta/2) * 
            grad_expepispline(epiparameters,r, dex, moment=2) %*% w
    }
    g2 <- g2 - softinfo$upperbound2moment
}
    else g2 <- NULL

if ( is.null (g1) &&  is.null (g2)) stop ("Called nonlinear; shouldn't have\n")
if (!is.null (g1) &&  is.null (g2))
   return (list (g = g1, grad_g = grad_g1))
if ( is.null (g1) && !is.null (g2))
   return (list (g = g2, grad_g = grad_g2))
if (!is.null (g1) && !is.null (g2))
    return (list (g = c(g1, g2), grad_g = cbind (grad_g1 = grad_g1, grad_g2 = grad_g2)))
}

integrate.to.one <- function(epiparameters, softinfo, r, force=FALSE)
{
#
# Impose the "integrate to one" constraint. This is an equality constraint!
#

N <- epiparameters$Ndiscr
m0 <- epiparameters$m0
mN <- epiparameters$mN
p <- epiparameters$order
Delta <- (mN-m0)/N

# Gaussian quadrature rule in each segment of order 5
# For [-1,1] integration domain; 
# For five points of evaluation
# Points of evaluation                 Weights
# 0                                     128/225
# (1/3)*sqrt(5-2*sqrt(10/7))            (322+13*sqrt(70))/900
# -(1/3)*sqrt(5-2*sqrt(10/7))           (322+13*sqrt(70))/900
# (1/3)*sqrt(5+2*sqrt(10/7))            (322-13*sqrt(70))/900
# -(1/3)*sqrt(5+2*sqrt(10/7))           (322-13*sqrt(70))/900
# Change of integration range from [-1,1] to [(k-1)*Delta + m0 , k*Delta + m0]
# int from (k-1)*Delta + m0 to k*Delta + m0 of expepispline(xi) dxi approx 
# (Delta/2)*sum i=1 to 5 w_i expepispline((Delta/2)*xi_i +
# ((k-1)*Delta + m0 + k*Delta + m0)/2)
# 
#w1 = (322+13*sqrt(70))/900;
#w2 = (322-13*sqrt(70))/900;
#w = [128/225; w1; w1; w2; w2]; 
#xi1 = (1/3)*sqrt(5-2*sqrt(10/7));
#xi2 = (1/3)*sqrt(5+2*sqrt(10/7));
#xivec = [0; xi1; -xi1; xi2; -xi2];

# For 20 points of evaluation:

w <- c(0.152753387, 0.152753387, 0.149172986, 0.149172986, 0.142096109,
      0.142096109, 0.131688638, 0.131688638, 0.118194532, 0.118194532,
      0.10193012,  0.10193012,  0.083276742, 0.083276742, 0.062672048,
      0.062672048, 0.04060143, 0.04060143, 0.017614007, 0.017614007)

xivec <- c(-0.076526521, 0.076526521, -0.227785851, 0.227785851, -0.373706089,
           0.373706089, 0.510867002, 0.510867002, -0.636053681, 0.636053681,
          -0.746331906, 0.746331906, -0.839116972, 0.839116972, -0.912234428,
           0.912234428, -0.963971927, 0.963971927, -0.993128599, 0.993128599)


geq <- 0
grad_geq <- rep (0, (p+2)* N + 1)



if  (force == TRUE || (any (names (softinfo) == "integrateToOne") &&
                 softinfo$integrateToOne == TRUE))
{

    for (k in 1:N) {
        dex <- (Delta/2)*xivec + ((k-1)*Delta + m0 + k*Delta + m0)/2
        geq <- geq + (Delta/2) * 
           sum(w * epispline(epiparameters, r, dex, exponentiate=TRUE))
# gradient 
        grad_geq = grad_geq + (Delta/2) * 
                   grad_expepispline(epiparameters, r, dex) %*% w
    }

    geq <- geq - 1

}
if (1 > 10) {
#    cat ("Inside integrate to one!\n")
#    cat ("Int-to-1: geq - 1 is ", geq, ", and the gradient is", grad_geq, "\n")
#    browser()
}
return (list (nonlinconst.eq.val = geq, nonlinconst.eq.grad = grad_geq))

}

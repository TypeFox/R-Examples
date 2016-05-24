# This file is intented to be called by calex_1d.R.  It specifies the
# true values for the parameters with the intent that we can try and
# estimate them from the generated dataset.

theta.TRUE <- 0.5

beta1.TRUE <- c(0,1,1)
names(beta1.TRUE) <- c("const" , "x", "A")

psi1.TRUE <- c(90,40,0.001)
names(psi1.TRUE) <- c("x","A","sigma1squared")

beta2.TRUE <- c(1,1)
names(beta2.TRUE) <- c("const","x.cubed")

psi2.TRUE <- c(80,0.001)
names(psi2.TRUE) <- c("x_rough","sigma2SQUARED")

lambda.TRUE <- 0.01

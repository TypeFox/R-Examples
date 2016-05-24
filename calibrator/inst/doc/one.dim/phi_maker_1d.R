# This file is intented to be called by calex_1d.R.  It generates a
# hyperparameter object phi.true.  Note that 

jj.psi1 <- 1:3
names(jj.psi1) <- c("x",  "A", "s1sq")

jj.psi2 <- 1:2
names(jj.psi2) <- c("x",         "s2sq")

jj.mean1 <- rep(1,3)
names(jj.mean1) <- names(jj.psi1)

jj.sigma1 <- diag(c(1.1, 1.1, 1.2))
rownames(jj.sigma1) <- names(jj.psi1)
colnames(jj.sigma1) <- names(jj.psi1)

jj.mean2 <- c(1,0.1,rep(1.1,2))
names(jj.mean2) <- c("rho","lambda",names(jj.psi2))
jj.sigma2 <- diag(c(1, 0.2, 0.2, 0.3))/10
rownames(jj.sigma2) <- names(jj.mean2)
colnames(jj.sigma2) <- names(jj.mean2)

jj.mean.th <- 1
names(jj.mean.th) <- c("A")
jj.sigma.th <- diag(c(1.5),nrow=1)
rownames(jj.sigma.th) <- names(jj.mean.th)
colnames(jj.sigma.th) <- names(jj.mean.th)


###################################################
### chunk number 13: 
###################################################
phi.1d <-
     phi.fun.1d(rho=1,
     lambda=1,
     psi1          = jj.psi1,
     psi2          = jj.psi2,
     psi1.apriori  = list(mean=jj.mean1,sigma=jj.sigma1),
     psi2.apriori  = list(mean=jj.mean2,sigma=jj.sigma2),
     theta.apriori = list(mean=jj.mean.th,sigma=jj.sigma.th)
                )

phi.1d2 <- phi.change(old.phi=phi.1d, phi.fun=phi.fun.1d, rho=3)

print(phi.1d2$rho)

phi.TRUE <- phi.change(phi.fun=phi.fun.1d, old.phi=phi.1d,
psi1=psi1.TRUE, psi2=psi2.TRUE,lambda=lambda.TRUE,rho=1)



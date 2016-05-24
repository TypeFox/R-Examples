### R code from vignette source 'calex.Rnw'

###################################################
### code chunk number 1: LoadLibrary
###################################################



###################################################
### code chunk number 2: calex.Rnw:162-163
###################################################
library("calibrator")


###################################################
### code chunk number 3: size_of_datasets
###################################################
do_from_scratch = FALSE
n1 <- 20
n2 <- 21


###################################################
### code chunk number 4: make_D1
###################################################
D1.int <- latin.hypercube(n1,4)
rownames(D1.int) <- paste("coderun",1:nrow(D1.int),sep=".")
colnames(D1.int) <- c("x", "y", "A", "B")
head(D1.int)


###################################################
### code chunk number 5: make_D2
###################################################
D2.int <- latin.hypercube(n2,2)
rownames(D2.int) <- paste("obs",1:nrow(D2.int),sep=".")
colnames(D2.int) <- c("x", "y")
head(D2.int)


###################################################
### code chunk number 6: make_extractor.int
###################################################
extractor.int <- function(D1){
return(list(x.star = D1[, 1:2, drop = FALSE], t.vec = D1[, 3:4, drop = FALSE]))
}


###################################################
### code chunk number 7: make_h1.int
###################################################
h1.int <- function(xin){
out <- c(1,xin[1],xin[3],xin[1]*xin[3], xin[2]^2*xin[4])
names(out) <- c("const" , "x", "A", "Ax" , "By.sq")
return(out)
}


###################################################
### code chunk number 8: make_H1.int
###################################################
H1.int <- function(D1){

    if (is.vector(D1)) {
        D1 <- t(D1)
    }
    out <- t(apply(D1, 1, h1.int))
    colnames(out) <- c("const" , "x", "A", "Ax" , "By.sq")
    return(out)
}


###################################################
### code chunk number 9: make_h2.int
###################################################
h2.int <- function(x){
out <- c(x[1],x[1]*x[2])
names(out) <- c("h2.x" ,"h2.xy")
return(out)
}

H2.int <- 
function (D2) 
{
    if (is.vector(D2)) {
        D2 <- t(D2)
    }
    out <- t(apply(D2, 1, h2.int))
    colnames(out) <- names(h2.int(D2[1, , drop = TRUE]))
    return(out)
}


###################################################
### code chunk number 10: make_pdm.maker.psi1
###################################################
    pdm.maker.psi1 <- function(psi1) {
        jj.omega_x <- diag(psi1[1:2])
        rownames(jj.omega_x) <- names(psi1[1:2])
        colnames(jj.omega_x) <- names(psi1[1:2])
        jj.omega_t <- diag(psi1[3:4],ncol=2)
        rownames(jj.omega_t) <- names(psi1[3:4])
        colnames(jj.omega_t) <- names(psi1[3:4])
        sigma1squared <- psi1[5]
        return(list(omega_x = jj.omega_x, omega_t = jj.omega_t, 
            sigma1squared = sigma1squared))
    }
    pdm.maker.psi2 <- function(psi2) {
        jj.omegastar_x <- diag(psi2[1:2],ncol=2)
        sigma2squared <- psi2[3]
        return(list(omegastar_x = jj.omegastar_x, sigma2squared = sigma2squared))
    }


###################################################
### code chunk number 11: phi.fun.int
###################################################
 phi.fun.int <- 
function (rho, lambda, psi1, psi1.apriori, psi2, psi2.apriori, 
    theta.apriori) 
{
    # CHANGES START
    pdm.maker.psi1 <- function(psi1) {
        jj.omega_x <- diag(psi1[1:2])
        rownames(jj.omega_x) <- names(psi1[1:2])
        colnames(jj.omega_x) <- names(psi1[1:2])
        jj.omega_t <- diag(psi1[3:4],ncol=2)
        rownames(jj.omega_t) <- names(psi1[3:4])
        colnames(jj.omega_t) <- names(psi1[3:4])
        sigma1squared <- psi1[5]
        return(list(omega_x = jj.omega_x, omega_t = jj.omega_t, 
            sigma1squared = sigma1squared))
    }
    pdm.maker.psi2 <- function(psi2) {
        jj.omegastar_x <- diag(psi2[1:2],ncol=2)
        sigma2squared <- psi2[3]
        return(list(omegastar_x = jj.omegastar_x, sigma2squared = sigma2squared))
    }
    # CHANGES END: remainder of function unaltered
    jj.mean <- theta.apriori$mean
    jj.V_theta <- theta.apriori$sigma
    jj.discard.psi1 <- pdm.maker.psi1(psi1)
    jj.omega_t <- jj.discard.psi1$omega_t
    jj.omega_x <- jj.discard.psi1$omega_x
    jj.sigma1squared <- jj.discard.psi1$sigma1squared
    jj.discard.psi2 <- pdm.maker.psi2(psi2)
    jj.omegastar_x <- jj.discard.psi2$omegastar_x
    jj.sigma2squared <- jj.discard.psi2$sigma2squared
    jj.omega_t.upper <- chol(jj.omega_t)
    jj.omega_t.lower <- t(jj.omega_t.upper)
    jj.omega_x.upper <- chol(jj.omega_x)
    jj.omega_x.lower <- t(jj.omega_x.upper)
    jj.a <- solve(solve(jj.V_theta) + 2 * jj.omega_t, solve(jj.V_theta, 
        jj.mean))
    jj.b <- t(2 * solve(solve(jj.V_theta) + 2 * jj.omega_t) %*% 
        jj.omega_t)
    jj.c <- jj.sigma1squared/sqrt(det(diag(nrow = nrow(jj.V_theta)) + 
        2 * jj.V_theta %*% jj.omega_t))
    names(jj.c) <- "ht.fun.precalc"
    jj.A <- solve(jj.V_theta + solve(jj.omega_t)/4)
    jj.A.upper <- chol(jj.A)
    jj.A.lower <- t(jj.A.upper)
    list(rho = rho, lambda = lambda, psi1 = psi1, psi1.apriori = psi1.apriori, 
        psi2 = psi2, psi2.apriori = psi2.apriori, theta.apriori = theta.apriori, 
        omega_x = jj.omega_x, omega_t = jj.omega_t, 
        omegastar_x = jj.omegastar_x, sigma1squared = jj.sigma1squared, 
        sigma2squared = jj.sigma2squared, omega_x.upper = jj.omega_x.upper, 
        omega_x.lower = jj.omega_x.lower, omega_t.upper = jj.omega_t.upper, 
        omega_t.lower = jj.omega_t.lower, a = jj.a, b = jj.b, 
        c = jj.c, A = jj.A, A.upper = jj.A.upper, A.lower = jj.A.lower)
}


###################################################
### code chunk number 12: make_jj.psi1
###################################################
jj.psi1 <- 1:5
names(jj.psi1) <- c("x", "y", "A","B", "s1sq")

jj.psi2 <- 1:3
names(jj.psi2) <- c("x", "y",          "s1sq")

jj.mean1 <- rep(1,5)
names(jj.mean1) <- names(jj.psi1)

jj.sigma1 <- diag(c(1.1, 1.1, 1.2, 1.3, 1.1))
rownames(jj.sigma1) <- names(jj.psi1)
colnames(jj.sigma1) <- names(jj.psi1)

jj.mean2 <- c(1,0.1,rep(1.1,3))
names(jj.mean2) <- c("rho","lambda",names(jj.psi2))
jj.sigma2 <- diag(c(1,0.2,1.1, 1.1, 1.2))/10
rownames(jj.sigma2) <- names(jj.mean2)
colnames(jj.sigma2) <- names(jj.mean2)

jj.mean.th <- 1:2
names(jj.mean.th) <- c("A","B")
jj.sigma.th <- diag(c(1.5,1.7))
rownames(jj.sigma.th) <- names(jj.mean.th)
colnames(jj.sigma.th) <- names(jj.mean.th)


###################################################
### code chunk number 13: make_phi.int
###################################################
phi.int <-
     phi.fun.int(rho=1,
     lambda=1,
     psi1          = jj.psi1,
     psi2          = jj.psi2,
     psi1.apriori  = list(mean=jj.mean1,sigma=jj.sigma1),
     psi2.apriori  = list(mean=jj.mean2,sigma=jj.sigma2),
     theta.apriori = list(mean=jj.mean.th,sigma=jj.sigma.th)
                 )



###################################################
### code chunk number 14: call_phi.change
###################################################
phi.int2 <- phi.change(old.phi=phi.int, phi.fun=phi.fun.int, rho=3)
print(phi.int2$rho)


###################################################
### code chunk number 15: examine_E.theta.toy
###################################################
E.theta.toy


###################################################
### code chunk number 16: examine_h1.int
###################################################
h1.int


###################################################
### code chunk number 17: make_E.theta.int
###################################################
E.theta.int <- function(D2 = NULL, H1 = NULL, x1 = NULL, x2 = NULL, phi, give.mean = TRUE)
{
    if (give.mean) {
        m_theta <- phi$theta.apriori$mean
        return(H1(D1.fun(D2, t.vec = m_theta)))
    }
    else {
        out <- matrix(0, 5, 5)
        out[3,3] <- phi$theta.apriori$sigma[1,1]
        out[3,4] <- phi$theta.apriori$sigma[1,1]*x1[1]
        out[4,3] <- phi$theta.apriori$sigma[1,1]*x2[1]
        out[4,4] <- phi$theta.apriori$sigma[1,1]*x1[1]*x2[1]
        out[5,5] <- phi$theta.apriori$sigma[2,2]*x1[2]*x2[2]
        return(out)
    }
}


###################################################
### code chunk number 18: make_Edash.theta.toy
###################################################
Edash.theta.int <- Edash.theta.toy


###################################################
### code chunk number 19: make_theta.TRUE
###################################################
theta.TRUE <- c(1,1)


###################################################
### code chunk number 20: make_beta1.TRUE
###################################################
beta1.TRUE <- c(0,0,0,1,1)
psi1.TRUE <- c(4,4,4,4,0.5)


###################################################
### code chunk number 21: make_two.designs
###################################################
two.designs <- rbind(D1.int,D1.fun(x.star=D2.int,t.vec=theta.TRUE))


###################################################
### code chunk number 22: examine_two.designs
###################################################
two.designs[c(1:3,(n1+n2):(n1+n2-2)),]


###################################################
### code chunk number 23: make_y.int
###################################################
jj.mean <- H1.int(two.designs)  %*%  beta1.TRUE
jj.sigma <- psi1.TRUE[5]*corr.matrix(two.designs, scales = psi1.TRUE[1:4])
code.and.obs <- as.vector(rmvnorm(n=1,mean=jj.mean,sigma=jj.sigma))
y.int <- code.and.obs[1:n1]
z.int <- code.and.obs[(n1+1):(n1+n2)]
names(y.int) <- rownames(D1.int)
head(y.int)


###################################################
### code chunk number 24: make_beta2.TRIE
###################################################
beta2.TRUE <- c(1,1)
psi2.TRUE <- c(3,3,0.6)


###################################################
### code chunk number 25: make_model.inadequacy
###################################################
jj.mean <-  drop(H2.int(D2.int) %*% beta2.TRUE)
jj.sigma <- corr.matrix(D2.int, scales=psi2.TRUE[1:2])*psi2.TRUE[3]
model.inadequacy <- rmvnorm(n=1, mean=jj.mean,sigma=jj.sigma)
z.int <- as.vector(z.int +  model.inadequacy) 
names(z.int) <- rownames(D2.int)


###################################################
### code chunk number 26: add_error
###################################################
lambda.TRUE <- 0.00
jj.obs.error <-  rnorm(n2)*lambda.TRUE
z.int <- z.int + jj.obs.error
head(z.int)


###################################################
### code chunk number 27: make_d.int
###################################################
d.int <- c(y.int , z.int)


###################################################
### code chunk number 28: make_psi1
###################################################
phi.true <- phi.change(phi.fun=phi.fun.int, old.phi=phi.int,
psi1=psi1.TRUE, psi2=psi2.TRUE,lambda=0.1,rho=1)


###################################################
### code chunk number 29: call_betahat.fun.koh
###################################################
betahat.fun.koh(theta=theta.TRUE, d=d.int, D1=D1.int, D2=D2.int, H1=H1.int, H2=H2.int, phi=phi.true)


###################################################
### code chunk number 30: call_betahat.fun.koh_otherparams
###################################################
betahat.fun.koh(theta=c(-5,5), d=d.int, D1=D1.int, D2=D2.int, H1=H1.int, H2=H2.int, phi=phi.true)


###################################################
### code chunk number 31: call_p.eqn8.supp
###################################################
p.eqn8.supp(theta=c(1,1), D1=D1.int, D2=D2.int, H1=H1.int, H2=H2.int, d=d.int, phi=phi.true)


###################################################
### code chunk number 32: call_p.eqn8.supp_difftheta
###################################################
p.eqn8.supp(theta=c(5,-6), D1=D1.int, D2=D2.int, H1=H1.int, H2=H2.int, d=d.int, phi=phi.true)


###################################################
### code chunk number 33: call_phi.stage1
###################################################
phi.stage1 <- stage1(D1=D1.int, y=y.int, H1=H1.int, maxit=10,
 method="SANN", trace=0, do.print=FALSE, phi.fun=phi.fun.int,
 phi=phi.int)


###################################################
### code chunk number 34: examine_stage1
###################################################
phi.stage1$psi1


###################################################
### code chunk number 35: examine_stage1_cont
###################################################
phi.stage1$sigma1squared
phi.stage1$omega_x
phi.stage1$omega_t


###################################################
### code chunk number 36: call_stage2
###################################################
  use1 <- 1:10
  use2 <- 1:11
  phi.stage2 <- stage2(D1=D1.int[use1,], D2=D2.int[use2,], H1=H1.int, H2=H2.int,
      y=y.int[use1], z=z.int[use2], extractor=extractor.int,
     phi.fun=phi.fun.int, E.theta=E.theta.int, Edash.theta=Edash.theta.int,
     maxit=1, method="SANN", phi=phi.stage1)


###################################################
### code chunk number 37: examine_stage2
###################################################
phi.stage2$rho
phi.stage2$lambda
phi.stage2$sigma2squared
phi.stage2$psi2


###################################################
### code chunk number 38: make_Xdist
###################################################
jj.xdist.mean <- rep(0.5,2)
names(jj.xdist.mean) <- c("x","y")
jj.xdist.var <- 0.05+diag(c(0.1,0.1))
rownames(jj.xdist.var) <- c("x","y")
colnames(jj.xdist.var) <- c("x","y")

X.dist.int <- list(mean=jj.xdist.mean,var=jj.xdist.var)


###################################################
### code chunk number 39: make_hbar.fun.int
###################################################
hbar.fun.int <- 
function (theta, X.dist, phi) 
{
    if (is.vector(theta)) {
        theta <- t(theta)
    }
    first.bit <- phi$rho * H1.int(D1.fun(X.dist$mean, theta))
    first.bit[,5] <- first.bit[,5] + theta[,2]*X.dist$var[2,2]
    second.bit <- H2.int(X.dist$mean)
    jj.names <- colnames(second.bit)
    second.bit <- kronecker(second.bit, rep(1, nrow(first.bit)))
    colnames(second.bit) <- jj.names
    return(t(cbind(first.bit, second.bit)))
}


###################################################
### code chunk number 40: call_EK.eqn10.theta
###################################################
  if(do_from_scratch){
    jj <- EK.eqn10.supp(X.dist=X.dist.int, D1=D1.int, D2=D2.int,
                        H1=H1.int, H2=H2.int, d=d.int, hbar.fun=hbar.fun.int,
                        lower.theta=c(-3,-3), upper.theta=c(3,3),
                        extractor=extractor.int,
                        phi=phi.stage2,eps=0.8)
  } else {
    jj <- 1.95633284138600
  }


###################################################
### code chunk number 41: printjj
###################################################
jj


###################################################
### code chunk number 42: defineMHwrapper
###################################################
p <- function(theta){
p.eqn8.supp(theta=theta, D1=D1.int, D2=D2.int, H1=H1.int, H2=H2.int, d=d.int, phi=phi.true)
}
p(c(1,1))
p(c(10,10))


###################################################
### code chunk number 43: sample.from.MH
###################################################
ss <- MH(n=10,start=c(1,1),sigma=diag(2),pi=p)


###################################################
### code chunk number 44: scatterplot
###################################################
plot(jitter(ss),xlab="A",ylab="B",main="Sample from posterior: try with more points")


###################################################
### code chunk number 45: SetTheBib
###################################################
bib <- system.file( "doc", "uncertainty.bib", package = "emulator" )
bib <- sub('.bib$','',bib)


###################################################
### code chunk number 46: usethebib
###################################################
cat( "\\bibliography{",bib,"}\n",sep='')



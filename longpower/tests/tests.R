library(longpower)

n = 3
t = c(0,2,5)
u = list(u1 = t, u2 = rep(0,n))
v = list(v1 = cbind(1,1,rep(0,n)),
         v2 = cbind(1,0,t))         
rho = c(0.2, 0.5, 0.8)
sigma2 = c(100, 200, 300)
tab = outer(rho, sigma2, 
      Vectorize(function(rho, sigma2){
        ceiling(diggle.linear.power(
          d=0.5,
          t=t,
          sigma2=sigma2,
          R=rho,
          alternative="one.sided",
          power=0.80)$n)}))
colnames(tab) = paste("sigma2 =", sigma2)
rownames(tab) = paste("rho =", rho)
tab

# reproduce Table 1 from Lu, Luo, & Chen (2008)
phi1 <- c(rep(1, 6), 2, 2)
phi2 <- c(1, 1, rep(2, 6))
lambda <- c(1, 2, sqrt(1/2), 1/2, 1, 2, 1, 2)
ztest <- ttest1 <- c()
for(i in 1:8){
  Na <- (phi1[i] + lambda[i] * phi2[i])*(qnorm(0.05/2) + qnorm(1-0.90))^2*(0.5^-2)
  Nb <- Na/lambda[i]
  ztest <- c(ztest, Na + Nb)
  v <- Na + Nb - 2
  Na <- (phi1[i] + lambda[i] * phi2[i])*(qt(0.05/2, df = v) + qt(1-0.90, df = v))^2*(0.5^-2)
  Nb <- Na/lambda[i]
  ttest1 <- c(ttest1, Na + Nb)
}
data.frame(phi1, phi2, lambda, ztest, ttest1)

# reproduce Table 2 from Lu, Luo, & Chen (2008)
tab <- c()
for(J in c(2,4))
for(aJ in (1:4)/10)
for(p1J in c(0, c(1, 3, 5, 7, 9)/10)){
  rJ <- 1-aJ
  r <- seq(1, rJ, length = J)
  # p1J = p^(J-1)
  tab <- c(tab, power.mmrm.ar1(rho = p1J^(1/(J-1)), ra = r, sigmaa = 1,
    lambda = 1, times = 1:J,
    delta = 1, sig.level = 0.05, power = 0.80)$phi1)
}
matrix(tab, ncol = 6, byrow = TRUE)

lmmpower(delta=1.5, t = seq(0,1.5,0.25),
	sig2.i = 55, sig2.s = 24, sig2.e = 10, cov.s.i=0.8*sqrt(55)*sqrt(24), power = 0.80)
lmmpower(n=208, t = seq(0,1.5,0.25),
	sig2.i = 55, sig2.s = 24, sig2.e = 10, cov.s.i=0.8*sqrt(55)*sqrt(24), power = 0.80)
lmmpower(beta = 5, pct.change = 0.30, t = seq(0,1.5,0.25),
	sig2.i = 55, sig2.s = 24, sig2.e = 10, cov.s.i=0.8*sqrt(55)*sqrt(24), power = 0.80)

library(lme4)
fm1 <- lmer(Reaction ~ Days + (Days|Subject), sleepstudy)
lmmpower(fm1, pct.change = 0.30, t = seq(0,9,1), power = 0.80)

library(nlme)
fm2 <- lme(Reaction ~ Days, random=~Days|Subject, sleepstudy)
lmmpower(fm2, pct.change = 0.30, t = seq(0,9,1), power = 0.80)

# random intercept only
fm3 <- lme(Reaction ~ Days, random=~1|Subject, sleepstudy)
lmmpower(fm3, pct.change = 0.30, t = seq(0,9,1), power = 0.80)

library(gee)
fm4 <- gee(Reaction ~ Days, id = Subject,
            data = sleepstudy,
            corstr = "exchangeable")
lmmpower(fm4, pct.change = 0.30, t = seq(0,9,1), power = 0.80)

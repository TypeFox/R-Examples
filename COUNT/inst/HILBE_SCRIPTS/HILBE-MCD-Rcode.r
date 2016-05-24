# Modeling Count Data, Cambridge University Press (2014)
# Joseph M Hilbe  hilbe@asu.edu

# ===============================================================
#           CHAPTER 1  Varieties of Count Data
# ===============================================================
#
# 1.2.1
# ------
#
# R CODE
sbp    <- c(131,132,122,119,123,115)
male   <- c(1,1,1,0,0,0)
smoker <- c(1,1,0,0,1,0)
age    <- c(34,36,30,32,26,23)
summary(reg1 <- lm(sbp~ male+smoker+age))
mu <- predict(reg1)
mu
cof <- reg1$coef
cof
xb <- cof[1] + cof[2]*male + cof[3]*smoker + cof[4]*age
xb
diff <- sbp - mu
diff

# 1.2.3
# --------

# Table 1.2a R: Code for Figure 1.2a
# =====================================================
obs <- 15; mu <- 4; y <- (0:140)/10; alpha <- .5
amu <- mu*alpha; layout(1)
all.lines <- vector(mode = 'list', length = 5)
for (i in 1:length(mu)) {
    yp = exp(-mu[i])*(mu[i]^y)/factorial(y)
    ynb1 = exp(  log(gamma(mu[i]/alpha + y))
               - log(gamma(y+1)) 
               - log(gamma(mu[i]/alpha)) 
               + (mu[i]/alpha)*log(1/(1+alpha)) 
               + y*log(1-1/(1+alpha)))
    ynb2 = exp(  y*log(amu[i]/(1+amu[i])) 
               - (1/alpha)*log(1+amu[i]) 
               + log( gamma(y +1/alpha) )
               - log( gamma(y+1) ) 
               - log( gamma(1/alpha) ))
    ypig = exp( (-(y-mu[i])^2)/(alpha*2*y*mu[i]^2)) * (sqrt(1/(alpha*2*pi*y^3)))
    ygp = exp( log((1-alpha)*mu[i]) 
               + (y-1)*log((1-alpha) * mu[i]+alpha*y)
               - (1-alpha)*mu[i] 
               - alpha*y 
               - log(gamma(y+1)))
    all.lines = list(yp = yp, ynb1 = ynb1, ynb2 = ynb2, ypig = ypig, ygp = ygp)
    ymax = max(unlist(all.lines), na.rm=TRUE)
    cols = c("red","blue","black","green","purple")
    plot(y, all.lines[[1]], ylim = 
         c(0, ymax), type = "n", main="5 Count Distributions: mean=4; alpha=0.5")
    for (j in 1:5)
        lines(y, all.lines[[j]], ylim = c(0, ymax), col=cols[j],type='b',pch=19, lty=j) 
         legend("topright",cex = 1.5, pch=19,
         legend=c("NB2","POI","PIG","NB1","GP"),
         col = c(1,2,3,4,5),
         lty = c(1,1,1,1,1),
         lwd = c(1,1,1,1,3))
}
# =======================================================

# 1.4.2
# -------

# Table 1.4 R: Poisson probabilities for y from 0 through 4
# ===========================================================
y <- c(4, 2, 0,3, 1, 2)
y0 <- exp(-2)* (2^0)/factorial(0)
y1 <- exp(-2)* (2^1)/factorial(1)
y2 <- exp(-2)* (2^2)/factorial(2)
y3 <- exp(-2)* (2^3)/factorial(3)
y4 <- exp(-2)* (2^4)/factorial(4)
poisProb <- c(y0, y1, y2, y3, y4); poisProb

# OR
dpois(0:4, lambda=2)

# CUMULATIVE
ppois(0:4, lambda=2)

# to plot a histogram
py <- 0:4
plot(poisProb ~ py, xlim=c(0,4), type="o", main="Poisson Prob 0-4: Mean=2")
# ============================================================


# Table 1.5 R : Code for Figure 1.3
# ===============================
m<- c(0.5,1,3,5)           #Poisson means 
y<- 0:11                   #Observed counts 
layout(1) 
for (i in 1:length(m)) { 
  p<- dpois(y, m[i])       #poisson pdf 
  if (i==1) { 
  plot(y, p, col=i, type='l', lty=i) 
  } else { 
  lines(y, p, col=i, lty=i) 
  } 
}
# ===============================


# ===============================================================
#              CHAPTER 2  Poisson Regression
# ===============================================================
#
# 2.3
# ----
#
# Table 2.4  R:  Synthetic Poisson Model
# =============================================
library(MASS); library(COUNT); set.seed(4590); nobs <- 50000
x1 <- runif(nobs); x2 <- runif(nobs); x3 <- runif(nobs)
py <- rpois(nobs, exp(1 + 0.75*x1 - 1.25*x2 + .5*x3))
cnt <- table(py)
dataf <- data.frame(prop.table(table(py) ) )
dataf$cumulative <- cumsum(dataf$Freq)
datafall <- data.frame(cnt, dataf$Freq*100, dataf$cumulative * 100)
datafall; summary(py)
summary(py1 <- glm(py ~ x1 + x2 + x3, family=poisson))
confint.default(py1); py1$aic/(py1$df.null+1)
r <- resid(py1, type = "pearson")
pchi2 <- sum(residuals(py1, type="pearson")^2)
disp <- pchi2/py1$df.residual; pchi2; disp
# =====================================================



# Table 2.6 R: Monte Carlo Poisson
# ==================================================
 mysim <- function()
 {
  nobs <- 50000
  x1 <- runif(nobs)
  x2 <- runif(nobs)
  x3 <- runif(nobs)
  py <- rpois(nobs, exp(2 + .75*x1 - 1.25*x2 + .5*x3))       
  poi <- glm(py ~ x1 + x2 + x3, family=poisson)  
    pr <- sum(residuals(poi, type="pearson")^2)
    prdisp <- pr/poi$df.residual
    beta <- poi$coef
    list(beta,prdisp)
}
B <- replicate(100, mysim())   
apply(matrix(unlist(B[1,]),4,100),1,mean)

# ===================================================
# Dispersion
mean(unlist(B[2,]))

# 2.4
# -----
# Table 2.7 R:  Example Poisson Model and Associated Statistics
# ===========================================================
library(COUNT)
data(rwm5yr); rwm1984 <- subset(rwm5yr, year==1984)
cage <- rwm1984$age - mean(rwm1984$age)
summary(poic <- glm(docvis ~ outwork + cage, family=poisson, data=rwm1984))
pr <- sum(residuals(poic, type="pearson")^2)  # Pearson Chi2
pr/poic$df.residual                           # dispersion statistic
modelfit(poic)
cnt <- table(rwm1984$docvis)
dataf <- data.frame(prop.table(table(rwm1984$docvis) ) )
dataf$cumulative <- cumsum(dataf$Freq)
datafall <- data.frame(cnt, dataf$Freq*100, dataf$cumulative * 100)
datafall
# ===========================================================

# Table 2.8 R: Change Levels in Categorical Predictor
# ======================================================
levels(rwm1984$edlevel)              # levels of edlevel
elevel <- rwm1984$edlevel            # new variable
levels(elevel)[2] <- "Not HS grad"   # assign level 1 to 2
levels(elevel)[1] <- "HS"            # rename level 1 to "HS"
levels(elevel)                       # levels of elevel
summary(tst2 <- glm(docvis ~ outwork + cage + female + married + kids 
    + factor(elevel), family=poisson, data=rwm1984))
# =======================================================


# 2.5.1
# -------

summary(pyq <- glm(docvis ~ outwork + age, family=poisson, data=rwm1984))

# Likelihood Profiling of SE
confint(pyq)

# Traditional Model-based SE
confint.default(pyq)


# 2.5.2
# ------

# Table 2.11  R: Poisson Model  - Rate Ratio Parameterization
# ============================================================
library(COUNT)
data(rwm5yr); rwm1984 <- subset(rwm5yr, year==1984)
summary(poi1 <- glm(docvis ~ outwork + age, family=poisson, data=rwm1984))
pr <- sum(residuals(poi1, type="pearson")^2)  # Pearson Chi2
pr/poi1$df.residual                    # dispersion statistic
poi1$aic / (poi1$df.null+1)            # AIC/n
exp(coef(poi1))                        # IRR
exp(coef(poi1))*sqrt(diag(vcov(poi1))) # delta method 
exp(confint.default(poi1))             # CI of IRR
# ============================================================


# 2.6
# ----

# Table 2.12 R:  Poisson with Exposure
# ===========================================
data(fasttrakg)
summary(fast <- glm(die ~ anterior + hcabg + factor(killip), 
                        family=poisson, 
                        offset=log(cases), 
                        data=fasttrakg))
exp(coef(fast))
exp(coef(fast))*sqrt(diag(vcov(fast))) 
exp(confint.default(fast))             
modelfit(fast)
# ============================================


# 2.7
# ----

# R  prediction
myglm <- glm(docvis ~ outwork + age, family=poisson, data=rwm1984)
lpred <- predict(myglm, newdata=rwm1984, type="link", se.fit=TRUE)
up <- lpred$fit + 1.96*lpred$se.fit; lo <- lpred$fit - 1.96*lpred$se.fit
eta <- lpred$fit ; mu <- myglm$family$linkinv(eta)
upci <- get(myglm$family$link)(up); loci <- get(myglm$family$link)(lo)


# 2.8.1
# -----

# Table 2.14 R:  Marginal Effects at Mean
# ===========================================================
library(COUNT)
data(rwm5yr); rwm1984 <- subset(rwm5yr, year==1984)
summary(pmem <- glm(docvis ~ outwork + age, family=poisson, data=rwm1984))
mout <- mean(rwm1984$outwork); mage <- mean(rwm1984$age)
xb <- coef(pmem)[1] + coef(pmem)[2]*mout + coef(pmem)[3]*mage
dfdxb <- exp(xb) * coef(pmem)[3]
mean(dfdxb)
# ===========================================================


# 2.8.2
# -------

# R CODE
mean(rwm1984$docvis) * coef(pmem)[3]

# 2.8.3
# -------

# R CODE  discrete change
summary(pmem <- glm(docvis ~ outwork + age, family=poisson, data=rwm1984))
mu0 <- exp(pmem$coef[1] + pmem$coef[3]*mage)
mu1 <- exp(pmem$coef[1] + pmem$coef[2] + pmem$coef[3]*mage)
pe <- mu1 - mu0
mean(pe)

# R CODE  avg partial effects
summary(pmem <- glm(docvis ~ outwork + age, family=poisson, data=rwm1984))
bout = coef(pmem)[2]
mu = fitted.values(pmem)
xb = pmem$linear.predictors
pe_out = 0
pe_out = ifelse(rwm1984$outwork == 0, exp(xb + bout)-exp(xb), NA)
pe_out = ifelse(rwm1984$outwork == 1, exp(xb)-exp(xb-bout),pe_out)
mean(pe_out) 


# ===============================================================
#              CHAPTER 3 Testing Overdispersion
# ===============================================================
#
# 3.1
# ----


# Table 3.1  R: Deviance Goodness-of-Fit Test
# ================================================================
library(COUNT); data(rwm5yr); rwm1984 <- subset(rwm5yr, year==1984)
mymod <-glm(docvis ~ outwork + age, family=poisson, data=rwm1984)
mymod

dev<-deviance(mymod); df<-df.residual(mymod)
p_value<-1-pchisq(dev,df)
print(matrix(c("Deviance GOF","D","df","p-value", " ", 
   round(dev,4),df, p_value), ncol=2))

# ==================================================================



# R CODE
mymod <-glm(docvis ~ outwork + age, family=poisson, data=rwm1984)
pr <- sum(residuals(mymod, type="pearson")^2)   # get Pearson Chi2
pchisq(pr, mymod$df.residual, lower=F)          # calc p-value 
pchisq(mymod$deviance, mymod$df.residual, lower= F)  # calc p-vl

# P__disp is now available when the COUNT pacakge is loaded. 
# Table 3.2  R: Function to Calculate Pearson Chi2 and Dispersion Statistics
# =============================================================
P__disp <- function(x) {
   pr <- sum(residuals(x, type="pearson")^2)
   dispersion <- pr/x$df.residual
   cat("\n Pearson Chi2 = ", pr , 
       "\n Dispersion   = ", dispersion, "\n")
}
# =============================================================

P__disp(mymod)

# R CODE
library(COUNT)
data(rwm5yr)
rwm1984 <- subset(rwm5yr, year==1984)
mymod <-glm(docvis ~ outwork + age, family=poisson, data=rwm1984)
P__disp(mymod)

# 3.3.1
# ----

# Table 3.4 R:  Z-Score Test
# ===========================================================
library(COUNT); data(rwm5yr); rwm1984 <- subset(rwm5yr, year==1984)
summary(poi <- glm(docvis ~ outwork + age, family=poisson, data=rwm1984))
mu <-predict(poi, type="response")
z <- ((rwm1984$docvis - mu)^2 - rwm1984$docvis)/ (mu * sqrt(2))
summary(zscore <- lm(z ~ 1))
# ==========================================================


# 3.3.2
# -----
# Table 3.5 R: Lagrange Multiplier Test
# ===============================================
obs <- nrow(rwm1984)   # continue from Table 3.2
mmu <- mean(mu); nybar <- obs*mmu; musq <- mu*mu
mu2 <- mean(musq)*obs
chival <- (mu2 - nybar)^2/(2*mu2); chival 
pchisq(chival,1,lower.tail = FALSE)
# ===============================================

# 3.3.3
# -----

# Table 3.7  R:  Poisson Model with Ancillary Statistics
# ============================================================
library(COUNT); data(rwm5yr); rwm1984 <- subset(rwm5yr, year==1984)
summary(poi1 <- glm(docvis ~ outwork + age, family=poisson, data=rwm1984))
pr <- sum(residuals(poi1, type="pearson")^2)  # Pearson Chi2
pr/poi1$df.residual                     # dispersion statistic
poi1$aic / (poi1$df.null+1)             # AIC/n
exp(coef(poi1))                         # IRR
exp(coef(poi1))*sqrt(diag(vcov(poi1)))  # delta method 
exp(confint.default(poi1))              # CI of IRR
modelfit(poi1)                          # same as Stata abic
sd(rwm1984$docvis)^2                    # observed variance
xbp <- predict(poi1)                    # xb, linear predictor
mup <- exp(xbp)                         # mu, fitted Poisson 
mean(mup)                               # expected variance: mean=variance
# Table of observed vs expected counts
rbind(obs=table(rwm1984$docvis)[1:18], 
     exp = round(sapply(0:17, function(x)sum(dpois(x, fitted(poi1))))))
meany <- mean(rwm1984$docvis)           # mean docvis
expect0 <- exp(-meany)*meany^0 / exp(log(factorial(0))) # expected prob of 0
zerodays <- (poi1$df.null+1) *expect0   # expected zero days
obs=table(rwm1984$docvis)[1:18]         # observed number values in each count 0-17
exp = round(sapply(0:17, function(x)sum(dpois(x, fitted(poi1))))) # expected each count
chisq.test(obs, exp)                    # ChiSq test if obs & exp from same pop
# ============================================================


# 3.4.1
# -----

data(medpar)

# R:  quasipoisson
# =============================================================
summary(poiql <- glm(los ~ hmo + white + hmo + factor(type), 
                           family=quasipoisson, data=medpar))
# =============================================================


# Table 3.8  R:  Scaling SE  Medpar Data
# =====================================================
library(COUNT); data(medpar); attach(medpar)
summary(poi <- glm(los ~ hmo + white + factor(type), family=poisson, 
    data=medpar)) 
confint(poi)                                    # profile confidence interval
pr <- sum(residuals(poi,type="pearson")^2)      # Pearson statistic
dispersion <- pr/poi$df.residual; dispersion          # dispersion
sse <- sqrt(diag(vcov(poi))) * sqrt(dispersion); sse  # model SE
# OR
poiQL <- glm(los ~ hmo + white + factor(type), family=quasipoisson, 
   data=medpar)
coef(poiQL); confint(poiQL)                     # coeff & scaled SEs
modelfit(poi)                                   # AIC,BIC statistics
# ======================================================


# 3.4.2
# ------

# Table 3.10  R: Quasi-likelihood Poisson Standard Errors
# ===============================================
poiQL <- glm(los ~ hmo+white+type2+type3, family=poisson, data=medpar)
summary(poiQL)
pr <-sum(residuals(poiQL, type="pearson")^2 ) 
disp <- pr/poiQL$df.residual        # Pearson dispersion
se <-sqrt(diag(vcov(poiQL)))
QLse <- se/sqrt(disp); QLse
# ===============================================


# 3.4.3
# ------

# Table 3.12  R:  Robust Standard Errors of medpar Model
# ====================================================
library(sandwich)
poi <- glm(los ~ hmo + white + factor(type), family=poisson, data=medpar)
vcovHC(poi)
sqrt(diag(vcovHC(poi, type="HC0")))                 # final HC0 = H-C-zero
# Clustering
poi <- glm(los ~ hmo + white + factor(type), family=poisson, data=medpar)
library(haplo.ccs)
sandcov(poi, medpar$provnum)
sqrt(diag(sandcov(poi, medpar$provnum)))
# ====================================================

summary(poi1 <- glm(los ~ hmo+white+factor(type), family=poisson, data=medpar))

# 3.4.4
# ------

# Table 3.13  R: Bootstrap Standard Errors
# ======================================================
library(COUNT); library(boot); data(medpar)
poi <- glm(los ~ hmo + white + factor(type), family=poisson, data=medpar)
summary(poi)
t <- function (x, i) {
xx <- x[i,]
 bsglm <- glm( los ~ hmo + white + factor(type), family=poisson, data=medpar)
 return(sqrt(diag(vcov(bsglm))))
 }
bse <- boot(medpar, t, R=1000)
sqrt(diag(vcov(poi))); apply(bse$t,2, mean)
# =======================================================

# ===============================================================
#              CHAPTER 4  Assessment of Fit
# ===============================================================
#
# 4.1

# R  Pearson Chi2 and statistic and graph
summary(pexp <- glm(docvis ~ outwork + cage, family=poisson, data=rwm1984))
presid <- residuals(pexp, type="pearson")
pchi2 <- sum(residuals(pexp, type="pearson")^2)  # Pearson Chi2
summary(rwm <- glm(docvis ~ outwork + age, family=poisson, data=rwm1984))
P__disp(rwm)
mu <- predict(rwm)
grd <- par(mfrow = c(2,2))
plot(x=mu, y= rwm$docvis, main = "Response residuals")
plot(x=mu, y= presid, main = "Pearson residuals")

# 4.2.1 
# -------
# Table 4.2 R: Likelihood Ratio Test
# ====================================================
library(COUNT); library(lmtest); data(rwm5yr)
rwm1984 <- subset(rwm5yr, year==1984)
poi1 <- glm(docvis ~ outwork + age, family=poisson, data=rwm1984)
poi1a <- glm(docvis ~ outwork, family=poisson, data=rwm1984)
lrtest(poi1, poi1a)
drop1(poi1, test="Chisq")
# ====================================================

# 4.2.2
# ------
# Chi2 test boundary LR 
pchisq(2.705,1, lower.tail=FALSE)/2

# 4.3.2
# ------

# Table 4.4 R:  Version of Stata User Command, abic 
# ======================================================== 
modelfit <- function(x) { 
obs <- x$df.null + 1 
aic <- x$aic 
xvars <- x$df.null - x$df.residual + 1
rdof <- x$df.residual 
aic_n <- aic/obs 
ll <- xvars - aic/2 
bic_r <- x$deviance - (rdof * log(obs)) 
bic_l <- -2*ll + xvars * log(obs) 
bic_qh <- -2*(ll - xvars * log(xvars))/obs 
c(AICn=aic_n, AIC=aic, BICqh=bic_qh, BICl=bic_l) 
} 
# modelfit(x) # substitute fitted model name for x 
# ======================================================== 
library(COUNT) 
data(medpar) 
mymodel <- glm(los ~ hmo + white + factor(type), family=poisson, data=medpar) 
modelfit(mymodel)

# ===============================================================
#              CHAPTER 5  Negative Binomial Regression
# ===============================================================
#
# 5.3.1

# Table 5.4  R: rwm1984 Modeling Example
# ============================================================
# make certain the appropriate packages are loaded
library(COUNT); data(rwm5yr); rwm1984 <- subset(rwm5yr, year==1984)
# USING glm.nb
summary(nbx <- glm.nb(docvis ~ outwork + age + married + female + 
       edlevel2 + edlevel3 + edlevel4, data=rwm1984))
exp(coef(nbx)); exp(coef(nbx))*sqrt(diag(vcov(nbx))) 
exp(confint.default(nbx))
alpha <- 1/nbx$theta; alpha; P__disp(nbx)
modelfit(nbx) 
xbnb <- predict(poi1); munb <- exp(xbnb)
# expected variance of NB model (using alpha where alpha=1/theta)
mean(munb)+ (1/nbx$theta)*mean(munb)^2
round(sqrt(rbind(diag(vcov(nbx)), diag(sandwich(nbx)))), digits=4)
# USING nbinomial
nb1 <- nbinomial(docvis ~ outwork + age + married + female + 
       edlevel2 + edlevel3 + edlevel4, data=rwm1984)
summary(nb1)
modelfit(nb1)
# ============================================================

# Table 5.5 R: rwm1984 Poisson and NB2 Models
# ============================================================
library(COUNT);library(msme)
data(rwm5yr);rwm1984 <- subset(rwm5yr, year==1984)
# POISSON
poi <- glm(docvis ~ outwork + age + married + female + 
                    edlevel2 + edlevel3 + edlevel4, 
                    family = poisson, data = rwm1984)
summary(poi)
#NB2
nb1 <- nbinomial(docvis ~ outwork + age + married + female + 
        edlevel2 + edlevel3 + edlevel4, data=rwm1984)
summary(nb1)
# NB1
library(gamlss)
summary(gamlss(formula = docvis ~ outwork + age + married + female +  
    edlevel2 + edlevel3 + edlevel4, family = NBII, data = rwm1984))
# ============================================================


# 5.4.3
# ------

# R CODE 
# POISSON
library(COUNT) ; data(nuts)
nut <- subset(nuts, dbh<.6)
sntrees <- scale(nut$ntrees)
sheight <- scale(nut$height)
scover  <- scale(nut$cover)
summary(PO <- glm(cones ~ sntrees + sheight + scover, family=quasipoisson, data=nut))
table(nut$cones)
summary(nut$cones)

# NEGATIVE BINOMIAL
library(msme)
NB <- nbinomial(cones ~ sntrees + sheight + scover, data=nut)
summary(NB)

# HETEROGENEOUS NEGATIVE BINOMIAL
summary(HNB <- nbinomial(cones ~ sntrees + sheight + scover, 
     formula2 =~ sntrees + sheight + scover, data=nut, family = "negBinomial", 
     scale.link = "log_s"))
exp(coef(HNB))


# ===============================================================
#        CHAPTER 6  Poisson Inverse Gaussian Regression
# ===============================================================
#
# 6.2.2
# ------

# Table 6.3  R: Poisson Inverse Gaussian - rwm1984
# ========================================================
library(gamlss); library(COUNT); library(msme); library(sandwich)
data(rwm5yr); rwm1984 <- subset(rwm5yr, year==1984)
summary(nbmod <- glm.nb(docvis ~ outwork + age, data=rwm1984))
vcovHC(nbmod)
sqrt(diag(vcovHC(nbmod, type="HC0")))  
pigmod <- gamlss(docvis ~ outwork + age, data=rwm1984, family=PIG)
summary(pigmod)
exp(coef(pigmod))
# =========================================================
exp(1.344)

# Table 6.5  R: Poisson Inverse Gaussian - medpar
# ========================================================
library(gamlss); library(COUNT); library(msme); library(sandwich)
data(medpar)
rwm1984 <- subset(rwm5yr, year==1984)
summary(nbmod1 <- glm.nb(los ~ hmo + white + factor(type), data=medpar))
vcovHC(nbmod1)
sqrt(diag(vcovHC(nbmod1, type="HC0")))  
pigmod1 <- gamlss(los ~ hmo + white + factor(type), data=medpar, family=PIG)
summary(pigmod1)
exp(coef(pigmod1))
# =========================================================


# ===============================================================
#        CHAPTER 7  Probleems with Zero Counts
# ===============================================================
#
# 
# R
# ================================================
exp(-3) * 3^0 / exp(log(factorial(0)))
100* (exp(-3) * 3^0 / exp(log(factorial(0))))
# ================================================

# 7.1.1
# ------

# Table 7.1 R:  Poisson and Zero-truncated Poisson 
# ====================================================
library(msme)
library(gamlss.tr)
data(medpar)
poi <- glm(los ~ white + hmo + factor(type), family=poisson, data=medpar)
summary(poi)
ztp <- gamlss(los ~ white + hmo + factor(type), data=medpar, family="PO")
gen.trun(0, "PO", type="left", name = "lefttr")
lt0poi <- gamlss(los~white+hmo+ factor(type), data=medpar, family="POlefttr")
summary(lt0poi)
# ====================================================

# 7.1.2
# ------

# Table 7.2   Zero-Truncated Negative Binomial
# ==========================================================
library(msme); library(gamlss.tr)
data(medpar)
nb <- nbinomial(los~ white + hmo + factor(type), data=medpar)
summary(nb)
ztnb <- gamlss(los~ white + hmo + factor(type),data=medpar, family="NBI")
gen.trun(0, "NBI", type="left", name = "lefttr")
lt0nb <- gamlss(los~white+hmo+ factor(type), data=medpar, family="NBIlefttr")
summary(lt0nb)
# ==========================================================

# R:  Calculate NB2 expected 0's for given ? and ?
# ==========================================
a <- 1; mu <- 2 ; y <- 0
 exp(y*log(a*mu/(1+a*mu))-(1/a)*log(1+a*mu)+ 
 log(gamma(y +1/a))-log(gamma(y+1))-log( gamma(1/a)))
# ==========================================

# R: Proof that sum of y probabilities from 0 to 100 is 1
# ================================================
a <- 1 ; mu <- 2 ; y <- 0:100
ff <- exp(y*log(a*mu/(1+a*mu))-(1/a)*log(1+a*mu)+ 
 log(gamma(y +1/a))-log(gamma(y+1))-log( gamma(1/a)))  
sum(ff)
# ================================================

# 7.2.1
# ------

# Table 7.3 R:  Poisson-Logit Hurdle
# ====================================================
library(pscl); library(COUNT) 
data(rwm5yr); rwm1984 <- subset(rwm5yr, year==1984)
hpl <- hurdle(docvis ~ outwork + age, dist="poisson", data=rwm1984, 
                       zero.dist="binomial", link="logit")
summary(hpl); AIC(hpl)
# ====================================================


# Table 7.4: R Components to Poisson-Logit Hurdle 
# ===================================================
visit <- ifelse(rwm1984$docvis >0, 1, 0)
table(visit)
logis <- glm(visit ~ outwork + age, data=rwm1984, 
                     family=binomial(link="logit"))
summary(logis)
library(pscl)
hpl2 <- hurdle(docvis ~ outwork + age, data=rwm1984,
    dist = "poisson", zero.dist="binomial", link="logit")
summary(hpl2)
logit <- glm(visit ~ outwork + age, data=rwm1984, 
             family=binomial(link="logit"))
summary(logit)
# ===================================================

# Table 7.5 R  NB2-logit Hurdle <Assume Model  from 7.3 Loaded>
# ======================================================
hnbl <- hurdle(docvis ~ outwork + age, dist="poisson", data=rwm1984, 
                       zero.dist="binomial", link="logit")
summary(hnbl); AIC(hnbl)
alpha <- 1/hnbl$theta ; alpha
exp(coef(hnbl))
predhnbl <- hnbl$fitted.values
# ======================================================

# Table 7.6.  R - ZIP
# ====================================================================
library(pscl); library(COUNT)
data(rwm5yr) ; rwm1984 <- subset(rwm5yr, year==1984)
poi <- glm(docvis ~ outwork + age, data=rwm1984, family=poisson)
zip <- zeroinfl(docvis ~ outwork + age | outwork + age, data=rwm1984, dist="poisson")
summary(zip)
print(vuong(zip,poi))
exp(coef(zip))
round(colSums(predict(zip, type="prob")[,1:17]))  # expected counts
rbind(obs=table(rwm1984$docvis)[1:18])            # observed counts
# =================================================================

# 7.3.5
# -----

# R CODE
pred <- round(colSums(predict(zip, type="prob") [,1:13]))
obs <- table(rwm1984$docvis)[1:13]
rbind(obs, pred)

# 7.3.6
# ------

# Table 7.7.  R - ZINB
# ===============================================================
library(pscl); library(COUNT)
data(rwm5yr); rwm1984 <- subset(rwm5yr, year==1984)
nb2 <- glm.nb(docvis ~ outwork + age, data=rwm1984)
zinb <- zeroinfl(docvis ~ outwork + age | outwork + age, data=rwm1984, dist="negbin")
summary(zinb)
print(vuong(zinb,nb2))
exp(coef(zinb))
pred <- round(colSums(predict(zinb, type="prob")[,1:13])) # expected counts
obs <- table(rwm1984$docvis)[1:13]                        # observed counts        
rbind(obs, pred)
# ====================================================================

# 7.3.7
# -----

# R  ZIPIG
# ======================================================
library(gamlss) ; data(rwm1984); attach(rwm1984)
zpig <- gamlss(docvis ~ outwork + age, sigma.fo= ~ -1, 
    family="ZIPIG", data=rwm1984)
# summary(zpig) Throws error; uncomment to see
# code for calculating vuong, LR test, etc on book's website
# ======================================================================

# ===============================================================
#        CHAPTER 8  Generalized Poisson
# ===============================================================
# ===============================================================
#        CHAPTER 9  More Advanced Models
# ===============================================================
#
# 9.1  exact models

# R  
# ==============================================
library(COUNT)
data(azcabgptca); attach(azcabgptca)
table(los); table(procedure, type); table(los, procedure)
summary(los)
summary(c91a <- glm(los ~ procedure+ type, family=poisson, data=azcabgptca))
modelfit(c91a)
summary(c91b <- glm(los ~ procedure+ type, family=quasipoisson, data=azcabgptca))
modelfit(c91b)
library(sandwich); sqrt(diag(vcovHC(c91a, type="HC0")))
# ==========================================================

# R
# ==============================================================
library(gamlss.tr)
gen.trun(0,"PO", type="left", name="leftr")
summary(c91c <- gamlss(los~ procedure+type, data=azcabgptca, family="POleftr"))
# ==============================================================================

# 9.2.1  Truncated models
# -----

# Table 9.2 R Left-Truncated at 3 Poisson
# ==================================================
library(COUNT); data(rwm5yr); rwm1984 <- subset(rwm5yr, year==1984)
summary(plt <- gamlss(docvis~outwork + age,data=rwm1984,family="PO"))
library(gamlss); library(gamlss.tr); pltvis<-subset(rwm1984, rwm1984$docvis>3)
summary(ltpo <- gamlss(docvis~outwork+age, family=trun(3, "PO", "left"), data=pltvis))
# -----------------
pltvis<-subset(rwm1984, rwm1984$docvis>3)    # alternative method
gen.trun(3, "PO", "left")                    # saved globally for session
summary(lt3po <- gamlss(docvis~outwork+age, family=POleft, data=pltvis))
# ==================================================

# Table 9.3 R: Right-Truncated Poisson : cut=10
# ==================================================
rtp<-subset(rwm1984, rwm1984$docvis<10)                 
summary(rtpo <- gamlss(docvis~outwork + age, data=rtp, 
   family=trun(10, "PO", type="right"))) # LT, not LE
# ==================================================


# 9.2.2
# ------

# Table 9.4 R  Left Censored Poisson at Cut=3
# =================================================
library(gamlss.cens); library(survival); library(COUNT)
data(rwm5yr); rwm1984 <- subset(rwm5yr, year==1984)
lcvis <- rwm1984
cy <- with(lcvis, ifelse(docvis<3, 3, docvis))  
ci <- with(lcvis, ifelse(docvis<=3, 0, 1))
Surv(cy,ci, type="left")[1:100]
cbind(Surv(cy,ci, type="left")[1:50], rwm1984$docvis[1:50])
lcmdvis <- data.frame(lcvis, cy, ci )
rm(cy,ci); gen.cens("PO",type="left")
lcat30<-gamlss(Surv(cy, ci, type="left") ~ outwork + age, 
    data=lcmdvis, family=POlc)
summary(lcat30)
# ================================================


# Table 9.5 R Right censored Poisson at 10
# =========================================
library(gamlss.cens); library(survival); rcvis <- rwm1984 
cy <- with(rcvis, ifelse(docvis>=10, 9, docvis)) 
ci <- with(rcvis, ifelse(docvis>=10, 0, 1)) 
rcvis <- data.frame(rcvis, cy, ci ) 
rm(cy,ci) ; gen.cens("PO",type="right") 
summary(rcat30<-gamlss(Surv(cy, ci) ~ outwork + age, 
   data=rcvis, family=POrc, n.cyc=100)) 
# ==========================================

# 9.3
# ------

# Table 9.6 R: Poisson-Poisson Finite Mixture Model
# =================================================
library(COUNT)
library(flexmix)
data(fishing)
attach(fishing)
## FIXME the following code gives an error that the model argument should
## be one of "gaussian", "binomial", "poisson", "Gamma" 

fmm_pg <- flexmix(totabund~meandepth + offset(log(sweptarea)),
                  data=rwm1984, k=2,
                  model=list(FLXMRglm(totabund~., family="poisson"),
                    FLXMRglm(totabund~., family="poisson")))
parameters(fmm_pg, component=1, model=1)
parameters(fmm_pg, component=2, model=1)
summary(fmm_pg)

# =================================================

# 9.4
# ------

#Table 9.7 R: GAM 
# ===========================================================
library(COUNT); library(mgcv)
data(rwm5yr); rwm1984 <- subset(rwm5yr, year==1984)
summary(pglm <- glm(docvis ~ outwork + age + female + married + 
     edlevel2 + edlevel3 + edlevel4, family=poisson, data=rwm1984))
summary(pgam <- gam(docvis ~ outwork + s(age) + female + married + 
     edlevel2 + edlevel3 + edlevel4, family=poisson, data=rwm1984))
plot(pgam)
# ===========================================================

# 9.6.1
#-----

# Table 9.8 R: GEE
# ================================================
library(COUNT); library(gee); data(medpar)
summary(pgee <- gee(los ~ hmo + white + age80 + type2 + type3, 
                  data=medpar,   id=medpar$provnum,
                  corstr='exchangeable', family=poisson))
# =============================================

# 9.6.2
# -----

# Table 9.9 R  Random Intercept Poisson
# ================================================
library(gamlss.mx)
summary(rip <- gamlssNP(los ~ hmo + white + type2 + type3, 
                          random=~1|provnum, data=medpar, 
                          family="PO", mixture="gq", K=20))
# ================================================

# 9.7
# ----

# Table 9.10  R  Generalized Waring Regression
# ========================================================
library(COUNT); library(GWRM)  
data(rwm5yr); rwm1984 <- subset(rwm5yr, year==1984)
war <- GWRM.fit(docvis ~ outwork + age + female + married, data=rwm1984)
GWRM.display(war)
# =========================================================

# 9.8
# ----

# Table 9.11 R: Bayesian Poisson  MCMC
# ======================================================
library(COUNT); library(MCMCpack); data(medpar)
summary(poi <- glm(los ~ hmo + white + type2 + type3, 
                   family=poisson, data=medpar))
confint.default(poi)
summary(poibayes <- MCMCpoisson(los ~ hmo + white + type2 + type3,
                burnin = 5000, mcmc = 100000, data=medpar))
# ======================================================

#-----------------------------------------------------------------------------#
#   Chapter 4 
#-----------------------------------------------------------------------------#

data(wool)

mod1 <- lm(log(y) ~ . , data = wool)
par(mfrow=c(1,2))
plot(fitted(mod1), residuals(mod1))
qqnorm(studres(mod1))
par(mfrow=c(1,1))

# a separate forward-search for each lambda value

lambda <- seq(-1,1,by=0.5)
mod2 <- fwdsco(y ~ x1 + x2 + x3, data = wool, lambda=lambda, nsamp="exact")
plot(mod2, plot.mle=F)
plot(mod2, plot.Sco=F, plot.Lik=T)

# computes the score test statistic using the data ordered by a search on untransformed data

mod3 <- fwdlm(y ~ x1 + x2 + x3, data = wool, nsamp="exact")
inc <- fwd.included(mod3)
y <- wool$y
x <- cbind(1,wool[,1:3])
score <- matrix(NA, length(inc), length(lambda))
for (i in 1:length(inc))
    for (j in 1:length(lambda))
        score[i,j] <- score.s(x[inc[[i]],], y[inc[[i]]], lambda[j])$Score

plot(0, 0, type="n", 
     ylim=range(score, na.rm=T), xlim=c(ncol(x),nrow(x)+1), 
     xlab="Subset Size", ylab="Score Test Statistic")
for (j in 1:length(lambda))
    { lines(ncol(x):nrow(x), score[,j], lty=j)
      text(nrow(x)+1, score[nrow(score), j], lambda[j])
    }
abline(h=c(1,-1)*2.58, col="lightgrey")

# forward plot with 95% likelihood intervals:

y <- wool$y
x <- cbind(1,wool[,1:3])
lambda <- seq(-2, 2, by=0.1)

lin.approx <- function(x,y,c)
{ approx <- spline(x, y, n=1000)
  j <- which(approx$y > c)
  c(mean(approx$x[c(min(j), min(j)-1)]), 
    lik.approx$x[which.max(approx$y)],
    mean(approx$x[c(max(j), max(j)+1)]))
}

# lambda = 1
inc <- fwd.included(fwdlm(y ~ x1 + x2 + x3, data = wool, nsamp="exact"))
mle <- matrix(NA, nrow=length(inc), 3)
for (j in 2:length(inc))
    { lik <- sapply(lambda, function(l, x, y, i) 
                                    score.s(x[i,], y[i], l)$Lik, 
                    x = x, y = y, i = inc[[j]] )
      mle[j,] <- lin.approx(lambda, lik, max(lik) - 1/2 * qchisq(0.95,1))
    }
plot(ncol(x):nrow(x), mle[,2], type="l",
     ylim=c(-1.1, 1.1), xlim=c(ncol(x),nrow(x)+1), 
     xlab="Subset Size", ylab="MLE of lambda (response = y)")
lines(ncol(x):nrow(x), mle[,1], lty=2)
lines(ncol(x):nrow(x), mle[,3], lty=2)
 
# lambda = 0
inc <- fwd.included(fwdlm(log(y) ~ x1 + x2 + x3, data = wool, nsamp="exact"))
mle <- matrix(NA, nrow=length(inc), 3)
for (j in 2:length(inc))
    { lik <- sapply(lambda, function(l, x, y, i) 
                                    score.s(x[i,], y[i], l)$Lik, 
                    x = x, y = y, i = inc[[j]] )
      mle[j,] <- lin.approx(lambda, lik, max(lik) - 1/2 * qchisq(0.95,1))
    }
plot(ncol(x):nrow(x), mle[,2], type="l",
     ylim=c(-1.1, 1.1), xlim=c(ncol(x),nrow(x)+1), 
     xlab="Subset Size", ylab="MLE of lambda (response = log(y) )")
lines(ncol(x):nrow(x), mle[,1], lty=2)
lines(ncol(x):nrow(x), mle[,3], lty=2)

#-----------------------------------------------------------------------------#

data(poison)
poison$poison <- as.factor(poison$poison)

mod <- lm(time ~ poison + treat, data = poison, x=T, y=T)
summary(mod)

mod1 <- fwdsco(time ~ poison + treat, data = poison, lambda=seq(-1,1,by=0.5))
plot(mod1, plot.mle=F)

sapply(mod1$Unit, function(x) x[(nrow(x)-5):nrow(x),1])


# Modified poison data

poison$time[8] <- 0.13
mod2 <- fwdsco(time ~ poison + treat, data = poison, lambda=seq(-1,1,by=0.5))
plot(mod2, plot.mle=F)

# search on untransformed data

mod3 <- fwdlm(time ~ poison + treat, data = poison)
inc <- fwd.included(mod3)
y <- poison$time
x <- mod$x
lambda <- seq(-1,1,by=0.5)
score <- matrix(NA, length(inc), length(lambda))
for (i in 1:length(inc))
    for (j in 1:length(lambda))
        score[i,j] <- score.s(x[inc[[i]],], y[inc[[i]]], lambda[j])$Score

plot(0, 0, type="n", 
     ylim=range(score, na.rm=T), xlim=c(ncol(x),nrow(x)+1), 
     xlab="Subset Size", ylab="Score Test Statistic")
for (j in 1:length(lambda))
    { lines(ncol(x):nrow(x), score[,j], lty=j)
      text(nrow(x)+1, score[nrow(score), j], lambda[j])
    }
abline(h=c(1,-1)*2.58, col="lightgrey")


# Doubly modified poison data

poison$time[38] <- 0.14

# single-deletion Score test statistic 
y <- poison$time
x <- mod$x
lambda <- c(-1,-0.5,0,1)
score <- matrix(NA, length(y), length(lambda))
for (i in 1:length(y))
    for (j in 1:length(lambda))
        score[i,j] <- score.s(x[-i,], y[-i], lambda[j])$Score
# full data Score test statistic (see also below)
score.all <- rep(NA, length(lambda))
for (j in 1:length(lambda))
    score.all[j] <- score.s(x, y, lambda[j])$Score
# diagnostic plots 
par(mfrow=c(2,2))
for (j in 1:length(lambda))
    { plot(0,0, type="n",
           ylim=c(min(score[,j])-0.3*diff(range(score[,j])),
                  max(score[,j])+0.3*diff(range(score[,j]))),
           xlim=c(1,nrow(score)), 
           xlab = "Case number", ylab = "Score test statistic",
           main = paste("lambda = ", lambda[j], sep=""))
      for (i in 1:nrow(score))
          segments(i, score.all[j], i, score[i,j], lwd=2)
      abline(h=score.all[j])
      abline(h=c(2.58, -2.58), col="lightgrey", lty=2)
    }
par(mfrow=c(1,1))

mod4 <- fwdsco(time ~ poison + treat, data = poison, lambda=seq(-1,1,by=0.5))
summary(mod4)
plot(mod4, plot.mle=F)

mod4$ScoreTest[nrow(mod4$ScoreTest),]


# Multiply modified poison data

poison <- read.table("data/poison.txt", header=T)
poison$poison <- as.factor(poison$poison)

poison$time[c(6,9,10,11)] <- c(0.14, 0.08, 0.07, 0.06)

# single-deletion Score test statistic 

y <- poison$time
x <- mod$x
lambda <- c(0, 1/3, 0.5)
score <- matrix(NA, length(y), length(lambda))
for (i in 1:length(y))
    for (j in 1:length(lambda))
        score[i,j] <- score.s(x[-i,], y[-i], lambda[j])$Score

# full data Score test statistic (see also below)

score.all <- rep(NA, length(lambda))
for (j in 1:length(lambda))
    score.all[j] <- score.s(x, y, lambda[j])$Score

# diagnostic plots 

par(mfrow=c(3,1))
for (j in 1:length(lambda))
    { plot(0,0, type="n",
           ylim=c(min(score[,j])-0.3*diff(range(score[,j])),
                  max(score[,j])+0.3*diff(range(score[,j]))),
           xlim=c(1,nrow(score)), 
           xlab = "Case number", ylab = "Score test statistic",
           main = paste("lambda = ", lambda[j], sep=""))
      for (i in 1:nrow(score))
          segments(i, score.all[j], i, score[i,j], lwd=2)
      abline(h=score.all[j])
      abline(h=c(1.96, -1.96), col="lightgrey", lty=2)
    }
par(mfrow=c(1,1))

mod5 <- fwdlm(time^(1/3) ~ poison + treat, data = poison) 
par(mfrow=c(1,3))
qqnorm(mod5$Residuals[,"m=6"], main="m=6"); qqline(mod5$Residuals[,"m=6"])
qqnorm(mod5$Residuals[,"m=27"], main="m=27"); qqline(mod5$Residuals[,"m=27"])
qqnorm(mod5$Residuals[,"m=48"], main="m=48"); qqline(mod5$Residuals[,"m=48"])
par(mfrow=c(1,1))

mod6 <- fwdsco(time ~ poison + treat, data = poison, lambda=seq(-1,1,by=0.5))
summary(mod6)
plot(mod6, plot.mle=F, ylim=c(-10,23))

mod6$ScoreTest[nrow(mod6$ScoreTest),]

mod7 <- fwdlm(time^(-1) ~ poison + treat, data = poison) 
plot(mod7,1)

#-----------------------------------------------------------------------------#

data(ozone)
ozone$Time <- 1:nrow(ozone)

mod1 <- fwdsco(y ~ Time + x2 + x5 + x6 + x8, data=ozone, lambda=seq(-1,1,by=0.5)) 
summary(mod1)
plot(mod1, plot.mle=F)

#-----------------------------------------------------------------------------#

data(stackloss)

mod1 <- fwdsco(Loss ~ Air + Temp + Conc, data=stackloss, nsamp="exact")
plot(mod1, plot.mle=F, ylim=c(-4,6))

mod2 <- fwdsco(Loss ~ Air + Temp, data=stackloss, nsamp="exact")
plot(mod2, plot.mle=F)
summary(mod2)

mod3 <- fwdlm(sqrt(Loss) ~ Air + Temp, data=stackloss, nsamp="exact")
par(mfrow=c(1,3))
#
qq <- qqnorm(mod3$Residuals[,"m=12"], main="m=12", xlim=c(-2.5,2.5))
j <- which(ifelse(qq$y < -2 | qq$y > 2, T, F))
text(qq$x[j], qq$y[j], j, pos=4)
qqline(mod3$Residuals[,"m=12"])
#
qq <- qqnorm(mod3$Residuals[,"m=19"], main="m=19", xlim=c(-2.5,2.5))
j <- which(ifelse(qq$y < -2 | qq$y > 2, T, F))
text(qq$x[j], qq$y[j], j, pos=4)
qqline(mod3$Residuals[,"m=19"])
#
qq <- qqnorm(mod3$Residuals[,"m=21"], main="m=21", xlim=c(-2.5,2.5))
j <- which(ifelse(qq$y < -2 | qq$y > 2, T, F))
text(qq$x[j], qq$y[j], j, pos=4)
qqline(mod3$Residuals[,"m=21"])
par(mfrow=c(1,1))



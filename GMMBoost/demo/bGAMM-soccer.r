library(GMMBoost)
data("soccer")
## generalized additive mixed model
## grid for the smoothing parameter


smoothpa <- seq(1,10001,by=1000)
AIC_vec <- rep(0,length(smoothpa))

## determination of optimal smoothing parameter
for (j in 1:length(smoothpa))
{
gamm <- bGAMM(points ~ ball.possession + tackles,
        ~ transfer.spendings + transfer.receits + unfair.score + ave.attend + sold.out,
        rnd = list(team=~1), data = soccer, lambda = smoothpa[j], family = poisson(link = log),
        control = list(overdispersion=TRUE,start=c(1,rep(0,25))))

AIC_vec[j] <- gamm$IC_sel[gamm$opt]
}

## final fit

gamm1 <- bGAMM(points ~ ball.possession + tackles,
         ~ transfer.spendings + transfer.receits + unfair.score + ave.attend + sold.out,
         rnd = list(team=~1), data = soccer, lambda = smoothpa[match(min(AIC_vec),AIC_vec)],
         family = poisson(link = log), control = list(overdispersion=TRUE,start=c(1,rep(0,25))))

plot(gamm1)

## generalized additive mixed model with random slopes and different control variables
## grid for the smoothing parameter
BIC_vec <- rep(0,length(smoothpa))

## determination of optimal smoothing parameter
for (j in 1:length(smoothpa))
{
gamm <- bGAMM(points ~ tackles + as.factor(yellow.red.card),
        ~ transfer.spendings + transfer.receits + sold.out + ave.attend,
        rnd = list(team=~1+tackles), data = soccer, family = poisson(link=log),lambda = smoothpa[j],
        control = list(nbasis=15,spline.degree=2,diff.ord=1,add.fix="ave.attend",
        overdispersion=TRUE,start=c(5,rep(0,52)),sel.method="bic",method="REML"))

BIC_vec[j] <- gamm$IC_sel[gamm$opt]
}

## final fit

gamm2 <- bGAMM(points ~ tackles + as.factor(yellow.red.card),
        ~ transfer.spendings + transfer.receits + sold.out + ave.attend,
        rnd = list(team=~1+tackles), data = soccer, family = poisson(link=log),
        lambda = smoothpa[match(min(BIC_vec),BIC_vec)],
        control = list(nbasis=15,spline.degree=2,diff.ord=1,add.fix="ave.attend",
        overdispersion=TRUE,start=c(5,rep(0,52)),sel.method="bic",method="REML"))

plot(gamm2)

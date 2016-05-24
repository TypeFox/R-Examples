## ----setup, include=FALSE------------------------------------------------
library(knitr)

## ----package, echo=FALSE, message=FALSE----------------------------------
library(mcprofile)

## ----, ctadata-----------------------------------------------------------
data(cta)
str(cta)
cta$concf <- factor(cta$conc, levels=unique(cta$conc))

ggplot(cta, aes(y=foci, x=concf)) + 
   geom_boxplot() +
   geom_dotplot(binaxis = "y", stackdir = "center", binwidth = 0.2) +    
   xlab("concentration")

## ----ctam1, message=FALSE------------------------------------------------
(m1 <- glm(foci ~ concf, data=cta, family=poisson(link="log")))
# confidence intervals for beta_j
confint(m1, parm=2:8)

## ----m1profile-----------------------------------------------------------
cm1 <- cbind(0, diag(7))
rownames(cm1) <- paste(levels(cta$concf)[-1], levels(cta$concf)[1], sep=" - ")
print(cm1)

p1 <- mcprofile(m1, cm1)
confint(p1, adjust="none")

## ----ctam2---------------------------------------------------------------
(m2 <- glm(foci ~ concf-1, data=cta, family=poisson(link="log")))

# m1 and m2 are the same models with different parameterization
all.equal(deviance(m1), deviance(m2))

## ----m2profile, message=FALSE--------------------------------------------
library(multcomp)
(cm2 <- contrMat(table(cta$concf), type="Dunnett"))

p2 <- mcprofile(m2, cm2)
(ci <- confint(p2, adjust="none"))

## ----plot1, warning=FALSE------------------------------------------------
plot(p2) + ylim(c(-5,5))

## ----plot2---------------------------------------------------------------
plot(ci)

## ----plot3---------------------------------------------------------------
plot(exp(ci)) + 
  coord_trans(x="log") +
  geom_vline(xintercept=1, lty=2)

## ----grid, warning=FALSE-------------------------------------------------
cgrid <- sapply(1:nrow(cm2), function(i) seq(-3, 5, length=15))
head(cgrid)

p3 <- mcprofile(m2, cm2, grid=cgrid)
plot(p3) 

## ----ci1-----------------------------------------------------------------
(ci2 <- confint(p2, adjust="single-step"))

## ----test1---------------------------------------------------------------
summary(p2, margin=1, alternative="greater")

## ----wald, warning=FALSE-------------------------------------------------
w2 <- wald(p2)
plot(w2) + ylim(c(-5,5))

## ----multcomp------------------------------------------------------------
confint(glht(m2, cm2))
confint(w2)

## ----hoa, warning=FALSE--------------------------------------------------
h2 <- hoa(p2)
confint(h2)

## ----toxin---------------------------------------------------------------
data(toxinLD)
toxinLD$logdose <- log(toxinLD$dose)
ggplot(toxinLD, aes(x=logdose, y=dead/(dead+alive))) + 
  geom_point(size=3)

## ----logisticregression--------------------------------------------------
mlr <- glm(cbind(dead, alive) ~ logdose, data=toxinLD, family=binomial(link="logit"))

## ----lrcm, message=FALSE, warning=FALSE----------------------------------
pdose <- seq(-1,2.3, length=7)
cmlr <- model.matrix(~ pdose)
rownames(cmlr) <- round(pdose, 2)
head(cmlr)

plr <- mcprofile(mlr, cmlr)
cilr <- confint(plr)

plot(cilr) +
  xlab("logit(mu)") +
  ylab("log dose level")

## ----lrciplot, message=FALSE, warning=FALSE------------------------------
pdat <- data.frame(logdose=pdose)
pdat$estimate <- mlr$family$linkinv(cilr$estimate$Estimate)
pdat$lower <- mlr$family$linkinv(cilr$confint$lower)
pdat$upper <- mlr$family$linkinv(cilr$confint$upper)

ggplot(toxinLD, aes(x=logdose, y=dead/(dead+alive))) + 
  geom_ribbon(data=pdat, aes(y=estimate, ymin=lower, ymax=upper), size=0.95, fill="lightblue") +
  geom_line(data=pdat, aes(y=estimate), size=0.95) +
  geom_point(size=3) +
  ylab("Mortality rate")


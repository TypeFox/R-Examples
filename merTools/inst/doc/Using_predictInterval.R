## ----setup, echo = FALSE, message=FALSE, warning=FALSE, results='hide'----
knitr::opts_chunk$set(
  cache=FALSE,
  comment="#>",
  collapse=TRUE, 
  echo = TRUE
)
library(knitr); library(merTools)

## ----Prep, message=FALSE, warning=FALSE----------------------------------
set.seed(271828)
data(sleepstudy)
fm1 <- lmer(Reaction ~ Days + (Days|Subject), data=sleepstudy)
display(fm1)

## ----predInt-------------------------------------------------------------
PI.time <- system.time(
  PI <- predictInterval(merMod = fm1, newdata = sleepstudy, 
                        level = 0.95, n.sims = 1000,
                        stat = "median", type="linear.prediction",
                        include.resid.var = TRUE)
)

## ----Inspect predInt, results="asis", echo=FALSE-------------------------
kable(head(PI))

## ----Inspect predInt 2, fig.width=7, fig.align="center"------------------
library(ggplot2);
ggplot(aes(x=1:30, y=fit, ymin=lwr, ymax=upr), data=PI[1:30,]) +
  geom_point() + 
  geom_linerange() +
  labs(x="Index", y="Prediction w/ 95% PI") + theme_bw()

## ----arm.Sim, fig.width=7, fig.height=4, fig.align="center"--------------
PI.arm.time <- system.time(
  PI.arm.sims <- arm::sim(fm1, 1000)
)

PI.arm <- data.frame(
  fit=apply(fitted(PI.arm.sims, fm1), 1, function(x) quantile(x, 0.500)),
  upr=apply(fitted(PI.arm.sims, fm1), 1, function(x) quantile(x, 0.975)),
  lwr=apply(fitted(PI.arm.sims, fm1), 1, function(x) quantile(x, 0.025))
)

comp.data <- rbind(data.frame(Predict.Method="predictInterval()", x=(1:nrow(PI))-0.1, PI),
                   data.frame(Predict.Method="arm::sim()", x=(1:nrow(PI.arm))+0.1, PI.arm))

ggplot(aes(x=x, y=fit, ymin=lwr, ymax=upr, color=Predict.Method), data=comp.data[c(1:30,181:210),]) +
  geom_point() + 
  geom_linerange() +
  labs(x="Index", y="Prediction w/ 95% PI") +
  theme_bw() +  theme(legend.position="bottom") +
  scale_color_brewer(type = "qual", palette = 2)
          
              

## ----bootMer.1, fig.width=7, fig.height=4, fig.align="center"------------
##Functions for bootMer() and objects
####Return predicted values from bootstrap
mySumm <- function(.) {
  predict(., newdata=sleepstudy, re.form=NULL)
}
####Collapse bootstrap into median, 95% PI
sumBoot <- function(merBoot) {
  return(
    data.frame(fit = apply(merBoot$t, 2, function(x) as.numeric(quantile(x, probs=.5, na.rm=TRUE))),
               lwr = apply(merBoot$t, 2, function(x) as.numeric(quantile(x, probs=.025, na.rm=TRUE))),
               upr = apply(merBoot$t, 2, function(x) as.numeric(quantile(x, probs=.975, na.rm=TRUE)))
    )
  )
}

##lme4::bootMer() method 1
PI.boot1.time <- system.time(
  boot1 <- lme4::bootMer(fm1, mySumm, nsim=1000, use.u=FALSE, type="parametric")
)

PI.boot1 <- sumBoot(boot1)

comp.data <- rbind(data.frame(Predict.Method="predictInterval()", x=(1:nrow(PI))-0.1, PI),
                   data.frame(Predict.Method="lme4::bootMer() - Method 1", x=(1:nrow(PI.boot1))+0.1, PI.boot1))

ggplot(aes(x=x, y=fit, ymin=lwr, ymax=upr, color=Predict.Method), data=comp.data[c(1:30,181:210),]) +
  geom_point() + 
  geom_linerange() +
  labs(x="Index", y="Prediction w/ 95% PI") +
  theme_bw() +  theme(legend.position="bottom") +
  scale_color_brewer(type = "qual", palette = 2)
  

## ----bootMer.2, fig.width=7, fig.height=4, fig.align="center"------------
##lme4::bootMer() method 2
PI.boot2.time <- system.time(
  boot2 <- lme4::bootMer(fm1, mySumm, nsim=1000, use.u=TRUE, type="parametric")
)

PI.boot2 <- sumBoot(boot2)

comp.data <- rbind(data.frame(Predict.Method="predictInterval()", x=(1:nrow(PI))-0.1, PI),
                   data.frame(Predict.Method="lme4::bootMer() - Method 2", x=(1:nrow(PI.boot2))+0.1, PI.boot2))

ggplot(aes(x=x, y=fit, ymin=lwr, ymax=upr, color=Predict.Method), data=comp.data[c(1:30,181:210),]) +
  geom_point() + 
  geom_linerange() +
  labs(x="Index", y="Prediction w/ 95% PI") +
  theme_bw() +  theme(legend.position="bottom") +
  scale_color_brewer(type = "qual", palette = 2)
  
  

## ----bootMer.3, fig.width=7, fig.height=4, fig.align="center"------------
##lme4::bootMer() method 3
PI.boot3.time <- system.time(
  boot3 <- lme4::bootMer(fm1, mySumm, nsim=1000, use.u=TRUE, type="semiparametric")
)

PI.boot3 <- sumBoot(boot3)

comp.data <- rbind(data.frame(Predict.Method="predictInterval()", x=(1:nrow(PI))-0.1, PI),
                   data.frame(Predict.Method="lme4::bootMer() - Method 3", x=(1:nrow(PI.boot3))+0.1, PI.boot3))

ggplot(aes(x=x, y=fit, ymin=lwr, ymax=upr, color=Predict.Method), data=comp.data[c(1:30,181:210),]) +
  geom_point() + 
  geom_linerange() +
  labs(x="Index", y="Prediction w/ 95% PI") +
  theme_bw() +  theme(legend.position="bottom") +
  scale_color_brewer(type = "qual", palette = 2)
  
  

## ---- echo=FALSE---------------------------------------------------------
times <- rbind(PI.time, PI.arm.time, PI.boot1.time, PI.boot2.time, PI.boot3.time)[,1:3]
rownames(times) <- c("predictInterval()", "arm::sim()", "lme4::bootMer()-Method 1", "lme4::bootMer()-Method 2", "lme4::bootMer()-Method 3")
kable(times)


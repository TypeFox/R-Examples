#' Simulate COST alike data
#' 
#' Simulate data alike the data from the Copenhagen stroke study (COST)
#' 
#' This uses functionality of the lava package.
#' 
#' @param N Sample size
#' @return Data frame
#' @author Thomas Alexander Gerds
#' @export
simCost <- function(N){
  requireNamespace("lava")
  ## psmT <- psm(Surv(time,status)~ age + sex + hypTen + prevStroke + othDisease + alcohol + diabetes + smoke + atrialFib + hemor + strokeScore + cholest,data=cost)
  ## psmCens <- psm(Surv(time,1-status)~1, data=cost)
  psmCens.scale <- 0.04330608
  psmT.scale <- 1.120334
  psmCens.coef <- c("(Intercept)"=8.297624)
  psmT.coef <- c("(Intercept)"=11.4457,"age"=-0.0572,"sex=male"=-0.4433,"hypTen=yes"=-0.2346,"prevStroke=yes"=-0.1873,"othDisease=yes"=-0.1327,"alcohol=yes"=0.1001,"diabetes=yes"=-0.333,"smoke=yes"=-0.3277,"atrialFib=yes"=-0.399,"hemor=yes"=0.0725,"strokeScore"=0.0250,"cholest"=-0.0049)
  m <- lava::lvm(~ age + sex + hypTen + prevStroke + othDisease + alcohol + diabetes + smoke + atrialFib + hemor + strokeScore + cholest + T + C)
  ## lava::distribution(m, ~age) <- lava::normal.lvm(mean=mean(cost$age),sd=sd(cost$age))
  lava::distribution(m, ~age) <- lava::normal.lvm(mean=73.29151,sd=11.32405)
  ## lava::distribution(m, ~sex) <- lava::binomial.lvm(p=mean(as.numeric(cost$sex=="male")))
  lava::distribution(m, ~sex) <- lava::binomial.lvm(p=0.465251)
  ## lava::distribution(m, ~strokeScore) <- normal.lvm(mean=mean(cost$strokeScore),sd=sd(cost$strokeScore))
  lava::distribution(m, ~strokeScore) <- lava::normal.lvm(mean=43.52896,sd=13.01235)
  ## lava::distribution(m, ~hypTen) <- lava::binomial.lvm(p=mean(as.numeric(cost$hypTen=="yes")))
  lava::distribution(m, ~hypTen) <- lava::binomial.lvm(p=0.3301158)
  ## lava::distribution(m, ~prevStroke) <- lava::binomial.lvm(p=mean(as.numeric(cost$prevStroke=="yes")))
  lava::distribution(m, ~prevStroke) <- lava::binomial.lvm(p=0.1833977)
  ## lava::distribution(m, ~othDisease) <- lava::binomial.lvm(p=mean(as.numeric(cost$othDisease=="yes")))
  lava::distribution(m, ~othDisease) <- lava::binomial.lvm(p=0.1621622)
  ## lava::distribution(m, ~alcohol) <- lava::binomial.lvm(p=mean(as.numeric(cost$alcohol=="yes")))
  lava::distribution(m, ~alcohol) <- lava::binomial.lvm(p=0.3166023)
  ## lava::distribution(m, ~diabetes) <- lava::binomial.lvm(p=mean(as.numeric(cost$diabetes=="yes")))
  lava::distribution(m, ~diabetes) <- lava::binomial.lvm(p=0.1409266)
  ## lava::distribution(m, ~smoke) <- lava::binomial.lvm(p=mean(as.numeric(cost$smoke=="yes")))
  lava::distribution(m, ~smoke) <- lava::binomial.lvm(p=0.4555985)
  ## lava::distribution(m, ~hemor) <- lava::binomial.lvm(p=mean(as.numeric(cost$hemor=="yes")))
  lava::distribution(m, ~hemor) <- lava::binomial.lvm(p=0.05019305)
  ## lava::distribution(m, ~atrialFib) <- lava::binomial.lvm(p=mean(as.numeric(cost$atrialFib=="yes")))
  lava::distribution(m, ~atrialFib) <- lava::binomial.lvm(p=0.1254826)
  lava::distribution(m,~T) <- eval(call("coxWeibull.lvm",scale=exp(-psmT.coef["(Intercept)"]/psmT.scale),shape=1/psmT.scale))
  lava::distribution(m,~C) <- eval(call("coxWeibull.lvm",scale=exp(-psmCens.coef["(Intercept)"]/psmCens.scale),shape=1/psmCens.scale))
  TCoef <- -psmT.coef[-1]/psmT.scale
  lava::regression(m) <- formula(paste("T ~ f(strokeScore,",TCoef[["strokeScore"]],") + f(age,",TCoef[["age"]],") + f(sex,",TCoef[["sex=male"]],") + f(hypTen,",TCoef[["hypTen=yes"]],") + f(prevStroke,",TCoef[["prevStroke=yes"]],") + f(othDisease,",TCoef[["othDisease=yes"]],") + f(alcohol,",TCoef[["alcohol=yes"]],") + f(diabetes,",TCoef[["diabetes=yes"]],") + f(smoke,",TCoef[["smoke=yes"]],") + f(atrialFib,",TCoef[["atrialFib=yes"]],") + f(hemor,",TCoef[["hemor=yes"]],")",sep=""))
  m <- lava::eventTime(m,time~min(T=1,C=0),"status")
  d <- lava::sim(m,N)
  d <- d[,-match(c("T","C"),names(d))]
  d
}

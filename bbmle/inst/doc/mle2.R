## ----knitropts,echo=FALSE,message=FALSE----------------------------------
if (require("knitr")) opts_chunk$set(fig.width=5,fig.height=5,tidy=FALSE,warning=FALSE,error=TRUE)

## ----setup,results="hide",echo=FALSE,message=FALSE-----------------------
library(Hmisc)

## ----emdbook,message=FALSE-----------------------------------------------
library(emdbook)

## ----bbsim---------------------------------------------------------------
set.seed(1001)
x1 <- rbetabinom(n=1000,prob=0.1,size=50,theta=10)

## ----bbmle,message=FALSE-------------------------------------------------
library("bbmle")

## ----likfun1-------------------------------------------------------------
mtmp <- function(prob,size,theta) {
  -sum(dbetabinom(x1,prob,size,theta,log=TRUE))
}

## ----fit1,warning=FALSE--------------------------------------------------
(m0 <- mle2(mtmp,start=list(prob=0.2,theta=9),data=list(size=50)))

## ----sum1----------------------------------------------------------------
summary(m0)

## ----confint1,warning=FALSE----------------------------------------------
confint(p0)
confint(m0,method="quad")
confint(m0,method="uniroot")

## ----profplot1,fig.height=5,fig.width=10,out.width="\\textwidth"---------
par(mfrow=c(1,2))
plot(p0,plot.confstr=TRUE)

## ----fit2,warning=FALSE--------------------------------------------------
m0f <- mle2(x1~dbetabinom(prob,size=50,theta),
            start=list(prob=0.2,theta=9),data=data.frame(x1))

## ----fit2f---------------------------------------------------------------
m0cf <- mle2(x1~dbetabinom(prob=plogis(lprob),size=50,theta=exp(ltheta)),
            start=list(lprob=0,ltheta=2),data=data.frame(x1))
confint(m0cf,method="uniroot")
confint(m0cf,method="spline")

## ----orobdata------------------------------------------------------------
load(system.file("vignetteData","orob1.rda",package="bbmle"))
summary(orob1)

## ----aodlikfun-----------------------------------------------------------
ML1 <- function(prob1,prob2,prob3,theta,x) {
  prob <- c(prob1,prob2,prob3)[as.numeric(x$dilution)]
  size <- x$n
  -sum(dbetabinom(x$m,prob,size,theta,log=TRUE))
}

## ----crowdertab,echo=FALSE,results="asis"--------------------------------
crowder.results <- matrix(c(0.132,0.871,0.839,78.424,0.027,0.028,0.032,-34.991,
                            rep(NA,7),-34.829,
                            rep(NA,7),-56.258),
                          dimnames=list(c("prop diffs","full model","homog model"),
                            c("prob1","prob2","prob3","theta","sd.prob1","sd.prob2","sd.prob3","NLL")),
                          byrow=TRUE,nrow=3)
latex(crowder.results,file="",table.env=FALSE,title="model")

## ----eval=FALSE----------------------------------------------------------
#  ## would prefer ~dilution-1, but problems with starting values ...
#  (m1B <- mle2(m~dbetabinom(prob,size=n,theta),
#               param=list(prob~dilution),
#               start=list(prob=0.5,theta=1),
#      data=orob1))

## ----suppWarn,echo=FALSE-------------------------------------------------
opts_chunk$set(warning=FALSE)

## ----aodstderr-----------------------------------------------------------
round(stdEr(m2),3)

## ----aodvar--------------------------------------------------------------
sqrt(1/(1+coef(m2)["theta"]))

## ----deltavar------------------------------------------------------------
sqrt(deltavar(sqrt(1/(1+theta)),meanval=coef(m2)["theta"],
         vars="theta",Sigma=vcov(m2)[4,4]))

## ----sigma3--------------------------------------------------------------
m2b <- mle2(m~dbetabinom(prob,size=n,theta=1/sigma^2-1),
            data=orob1,
            parameters=list(prob~dilution,sigma~1),
            start=list(prob=0.5,sigma=0.1))
## ignore warnings (we haven't bothered to bound sigma<1)
round(stdEr(m2b)["sigma"],3)
p2b <- profile(m2b,prof.lower=c(-Inf,-Inf,-Inf,0))

## ----compquad------------------------------------------------------------
r1 <- rbind(confint(p2)["theta",],
            confint(m2,method="quad")["theta",])
rownames(r1) <- c("spline","quad")
r1

## ----profplottheta-------------------------------------------------------
plot(p2,which="theta",plot.confstr=TRUE)

## ----profplotsigma-------------------------------------------------------
plot(p2b,which="sigma",plot.confstr=TRUE,
     show.points=TRUE)

## ----homogmodel----------------------------------------------------------
ml0 <- function(prob,theta,x) {
  size <- x$n
  -sum(dbetabinom(x$m,prob,size,theta,log=TRUE))
}
m0 <- mle2(ml0,start=list(prob=0.5,theta=100),
          data=list(x=orob1))

## ----logLikcomp----------------------------------------------------------
logLik(m0)

## ----formulafit----------------------------------------------------------
m0f <- mle2(m~dbetabinom(prob,size=n,theta),
            parameters=list(prob~1,theta~1),
            data=orob1,
            start=list(prob=0.5,theta=100))
m2f <- update(m0f,
              parameters=list(prob~dilution,theta~1),
              start=list(prob=0.5,theta=78.424))
m3f <- update(m0f,
              parameters=list(prob~dilution,theta~dilution),
              start=list(prob=0.5,theta=78.424))

## ----anovafit------------------------------------------------------------
anova(m0f,m2f,m3f)

## ----ICtabfit------------------------------------------------------------
AICtab(m0f,m2f,m3f,weights=TRUE,delta=TRUE,sort=TRUE)
BICtab(m0f,m2f,m3f,delta=TRUE,nobs=nrow(orob1),sort=TRUE,weights=TRUE)
AICctab(m0f,m2f,m3f,delta=TRUE,nobs=nrow(orob1),sort=TRUE,weights=TRUE)

## ----reWarn,echo=FALSE---------------------------------------------------
opts_chunk$set(warning=FALSE)

## ----frogsetup-----------------------------------------------------------
frogdat <- data.frame(
  size=rep(c(9,12,21,25,37),each=3),
  killed=c(0,2,1,3,4,5,rep(0,4),1,rep(0,4)))
frogdat$initial <- rep(10,nrow(frogdat))

## ----getgg---------------------------------------------------------------
library(ggplot2)

## ----gg1-----------------------------------------------------------------
gg1 <- ggplot(frogdat,aes(x=size,y=killed))+geom_point()+
      stat_sum(aes(size=factor(..n..)))+
      labs(size="#")+scale_x_continuous(limits=c(0,40))

## ----gg1plot-------------------------------------------------------------
gg1 + geom_line(data=pdat1,colour="red")+
      geom_line(data=pdat2,colour="blue")

## ----basegraphprofplot---------------------------------------------------
plot(prof4)

## ----latticeprof,fig.height=5,fig.width=10,out.width="\\textwidth"-------
prof4_df <- as.data.frame(prof4)
library(lattice)
xyplot(abs(z)~focal|param,data=prof4_df,
      subset=abs(z)<3,
       type="b",
       xlab="",
       ylab=expression(paste(abs(z),
           " (square root of ",Delta," deviance)")),
       scale=list(x=list(relation="free")),
             layout=c(3,1))

## ----ggplotprof,fig.height=5,fig.width=10--------------------------------
ss <-subset(prof4_df,abs(z)<3)
ggplot(ss,
       aes(x=focal,y=abs(z)))+geom_line()+
      geom_point()+
      facet_grid(.~param,scale="free_x")

## ----oldargs,eval=FALSE--------------------------------------------------
#  function (x, levels, conf = c(99, 95, 90, 80, 50)/100, nseg = 50,
#            absVal = TRUE, ...) {}

## ----newargs,eval=FALSE--------------------------------------------------
#  function (x, levels, which=1:p, conf = c(99, 95, 90, 80, 50)/100, nseg = 50,
#            plot.confstr = FALSE, confstr = NULL, absVal = TRUE, add = FALSE,
#            col.minval="green", lty.minval=2,
#            col.conf="magenta", lty.conf=2,
#            col.prof="blue", lty.prof=1,
#            xlabs=nm, ylab="score",
#            show.points=FALSE,
#            main, xlim, ylim, ...) {}


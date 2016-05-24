## ----setup,include=FALSE-------------------------------------------------
library(rotations)
options(digits=3)
library(knitr)
opts_chunk$set(fig.path='figure/', fig.width=5, fig.height=5, fig.show='hold',fig.align='center')

## ----setup2,include=FALSE------------------------------------------------
options(digits=3)

## ----ex2-----------------------------------------------------------------
r <- pi/2
U <- c(0, 1 ,0)
W <- U*r
R <- as.SO3(W)
R
identical(R, as.SO3(U, r))

## ----ex3-----------------------------------------------------------------
mis.angle(R)*2/pi
mis.axis(R)

## ----ex4-----------------------------------------------------------------
as.Q4(U, r)
as.Q4(as.SO3(U, r))

## ----ex5-----------------------------------------------------------------
Rs <- ruars(n = 20, rangle = rcayley, kappa = 1, space = 'SO3')
Qs <- ruars(n = 20, rangle = rcayley, kappa = 1, space = 'Q4')
Rs <- ruars(n = 20, rangle = rcayley, nu = 1, space = 'SO3')
Qs <- ruars(n = 20, rangle = rcayley, nu = 1, space = 'Q4')
head(Rs,3)

## ----gridsearch----------------------------------------------------------
# error function definition
L1.error <- function(sample, Shat) {
    sum(rot.dist(sample, Shat, method = "intrinsic", p = 1))
}
cayley.sample <- ruars(n = 50, rangle = rcayley, nu = 1, space = "SO3")
# gradient based optimization
system.time(SL1 <- gradient.search(cayley.sample, L1.error))
# in-built function
system.time(S <- median(cayley.sample, type = "geometric"))
rot.dist(S, SL1$Shat)

## ----ex6-----------------------------------------------------------------
Rs <- ruars(50, rcayley, kappa = 10)
region(Rs, method="direct", type="asymptotic", estimator="mean", alp=0.05)
region(Rs, method="direct", type="bootstrap", estimator="mean", alp=0.05, m=300)
region(Rs, method="direct", type="asymptotic", estimator="median", alp=0.05)
region(Rs, method="direct", type="bootstrap", estimator="median", alp=0.05, m=300)

## ----ex7,fig.cap="The $x$-axis of a random sample from the Cayley-UARS distribution with $\\kappa=1$, $n=50$.  All for point estimates are displayed on the left and all three region methods along with the projected mean are on the right.",fig.lp="figure:eye1",out.width=".4\\textwidth",fig.pos="h",dev='png'----
plot(Rs, center = mean(Rs), show_estimates = "all")
plot(Rs, center = mean(Rs), show_estimates = "proj.mean", 
 mean_regions = "all",  alp = .05)

## ----ex1-----------------------------------------------------------------
data(drill)
head(drill)
Subj1Wrist<-subset(drill, Subject == '1' & Joint == 'Wrist')
Subj1Wdata <- as.Q4(Subj1Wrist[,5:8])
mean(Subj1Wdata)

## ----ex21----------------------------------------------------------------
data(nickel)
head(nickel[,1:6])
Location1<-subset(nickel, location==1)
Loc1data<-as.SO3(Location1[,5:13])
mean(Loc1data)

## ----summary-------------------------------------------------------------
 Qs<-ruars(20, rcayley, space='Q4')
 Rs<-as.SO3(Qs)
 suppressMessages(require(onion))
 onionQs <- as.quaternion(t(Qs))
 suppressMessages(require(orientlib))
 orientRs <- rotvector(matrix(Rs, ncol = 9))


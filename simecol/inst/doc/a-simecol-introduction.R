### R code from vignette source 'a-simecol-introduction.Rnw'

###################################################
### code chunk number 1: init
###################################################
library("simecol")
data(lv, package = "simecol")
options("width"=72)
options("prompt" = "R> ", "continue" = "+  ")


###################################################
### code chunk number 2: out (eval = FALSE)
###################################################
## library("simecol")
## data(lv, package = "simecol")
## lv <- sim(lv)
## plot(lv)
## o1 <- out(lv)


###################################################
### code chunk number 3: parm1
###################################################
parms(lv)
parms(lv)["k1"] <- 0.4


###################################################
### code chunk number 4: parm2
###################################################
parms(lv)["a"] <- 1
parms(lv)


###################################################
### code chunk number 5: parms3
###################################################
parms(lv) <- parms(lv)[-4]
parms(lv)


###################################################
### code chunk number 6: CA (eval = FALSE)
###################################################
## library("simecol")
## data(CA, package="simecol")
## CA <- sim(CA)
## plot(CA)


###################################################
### code chunk number 7: CA (eval = FALSE)
###################################################
## times(CA)
## times(CA) <- c(to=100)
## CA <- sim(CA)
## plot(CA)


###################################################
### code chunk number 8: lv1
###################################################
library("simecol")
data(lv, package="simecol")
lv1 <- lv2 <- lv


###################################################
### code chunk number 9: lv2
###################################################
init(lv1)
parms(lv1)
parms(lv2)["k3"] <- 0.1
lv1 <- sim(lv1)
lv2 <- sim(lv2)


###################################################
### code chunk number 10: lv12
###################################################
o1 <- out(lv1); o2 <- out(lv2)
par(mfrow=c(1,2), las=1)
plot(o1$time, o1$prey, ylab="Prey, Predator",
  xlab="Time", col="blue", type="l", ylim=c(0,2), main="Scenario 1")
lines(o1$time, o1$predator, col="red", lty="dashed")
plot(o2$time, o2$prey, ylab="Prey, Predator",
  xlab="Time", col="blue", type="l", ylim=c(0,2), main="Scenario 2")
lines(o2$time, o2$predator, col="red", lty="dashed")
legend(10, 2, c("Prey", "Predator"), col=c("blue", "red"),
  lty=c("solid", "dashed"), lwd=1)


###################################################
### code chunk number 11: lv_range
###################################################
sapply(o1[c("predator", "prey")], range)
sapply(o2[c("predator", "prey")], range)


###################################################
### code chunk number 12: lv_spectrum
###################################################
tlv <- times(lv1)
ots <- ts(o1[c("predator", "prey")], start=tlv["from"],
          end=tlv["to"], deltat=tlv["by"])
sp  <- spectrum(ots, spans=c(3,3), log="no")
1/sp$freq[sp$spec[,1] == max(sp$spec[,1])]


###################################################
### code chunk number 13: lv_ef1
###################################################
lv_ef <- lv

main(lv_ef) <-  function (time, init, parms, ...) {
  x <- init
  p <- parms
  S <- approxTime1(inputs, time, rule=2)["s.in"]
  dx1 <-   S * p["k1"] * x[1] - p["k2"] * x[1] * x[2]
  dx2 <- - p["k3"] * x[2] + p["k2"] * x[1] * x[2]
  list(c(dx1, dx2))
}


###################################################
### code chunk number 14: lv_ef2
###################################################
inputs(lv_ef) <-  as.matrix(data.frame(
  time = c(0, 30, 30.1, 100),
  s.in = c(0,  0,  .5,     .5)
))


###################################################
### code chunk number 15: lv_eff3
###################################################
o <- out(sim(lv_ef))
matplot(o$time,o[2:3], xlab="Time",
  ylab="Substrate, Prey, Predator", type="l",
  lty=c("solid", "dashed"), col=c("blue", "red"), las=1)
inp <- as.data.frame(inputs(lv_ef))
lines(inp$time, inp$s.in, col="darkgreen", lwd=2, lty="11")


###################################################
### code chunk number 16: a-simecol-introduction.Rnw:818-819
###################################################
options("prompt" = " ", "continue" = " ")


###################################################
### code chunk number 17: lv_efr
###################################################
lv_efr <- lv_ef
tt     <- fromtoby(times(lv_efr))
o      <- matrix(0, nrow=length(tt), ncol=10)
initfunc(lv_efr) <- function(obj) {
  tt <- fromtoby(times(obj))
  inputs(obj) <- as.matrix(data.frame(
    time = tt,
    s.in = pmax(rnorm(tt, mean=1, sd=0.5), 0)
  ))
  obj
}
for (i in 1:10) {
  lv_efr <- initialize(lv_efr)
  lv_efr <- sim(lv_efr)
  o[,i] <- out(lv_efr)$prey
}
matplot(tt, o, xlab="Time", ylab="Prey", las=1, type="l")


###################################################
### code chunk number 18: a-simecol-introduction.Rnw:843-844
###################################################
options("prompt" = "R> ", "continue" = "+ ")


###################################################
### code chunk number 19: ibm_class (eval = FALSE)
###################################################
## setClass("indbasedModel",
##   representation(
##     parms  = "list",
##     init   = "data.frame"
##   ),
##   contains = "simObj"
## )


###################################################
### code chunk number 20: ibm_daphnia
###################################################
source("ibm_daphnia.R")


###################################################
### code chunk number 21: sim_ibm_daphnia (eval = FALSE)
###################################################
## solver(ibm_daphnia) <- "iteration"
## ibm_daphnia <- sim(ibm_daphnia)


###################################################
### code chunk number 22: ibm_plot (eval = FALSE)
###################################################
## setMethod("plot", c("indbasedModel", "missing"), function(x, y, ...) {
##   o <- out(x)
##   par(mfrow=c(2, 2))
##   plot(o$times, o$meanage,  type="l", xlab="Time", ylab="Mean age (d)")
##   plot(o$times, o$meaneggs, type="l", xlab="Time", ylab="Eggs per individual")
##   plot(o$times, o$number,   type="l", xlab="Time", ylab="Abundance")
##   plot(o$times, o$number,   type="l", xlab="Time", ylab="Abundance", log="y")
## })


###################################################
### code chunk number 23: daphnia1
###################################################
solver(ibm_daphnia) <- "myiteration"
ibm_daphnia <- sim(ibm_daphnia)
plot(ibm_daphnia)


###################################################
### code chunk number 24: ibm_init
###################################################
Sc0                                   <- ibm_daphnia
times(Sc0)                            <- c(from=0, to=30, by=0.2)
parms(Sc0)[c("temp", "food", "mort")] <- c(15, 0.4, 0.1)
init(Sc0) <- data.frame(age=rep(10, 50), size = rep(2.5, 50),
                        eggs=rep(5, 50), eggage=runif(50, 0, 4))


###################################################
### code chunk number 25: ibm_eqations
###################################################
equations(Sc0)$survive = function(inds, parms) {
  abundance <- nrow(inds)
  rnd       <- runif(abundance)
  mort      <- fmort(parms$mort, inds$size) * parms$DELTAT
  subset(inds, rnd > mort)
}


###################################################
### code chunk number 26: ibm_mort
###################################################
Sc1 <- Sc2 <- Sc3 <- Sc0
equations(Sc0)$fmort <- function(mort, x) 0
equations(Sc1)$fmort <- function(mort, x) mort
equations(Sc2)$fmort <- function(mort, x){
  mort * 2 * rank(-x) / (length(x) + 1)
}
equations(Sc3)$fmort <- function(mort, x){
  mort * 2 * rank(x) / (length(x) + 1)
}


###################################################
### code chunk number 27: ibm_seed
###################################################
set.seed(123)


###################################################
### code chunk number 28: ibm_sim_sc
###################################################
sc <- lapply(list(Sc0=Sc0, Sc1=Sc1, Sc2=Sc2, Sc3=Sc3), sim)


###################################################
### code chunk number 29: daphnia2
###################################################
growthrate <- function(obj) {
  o <- subset(out(obj), times > 10)
  m <- lm(log(o$number) ~ o$times)
  as.vector(coef(m)[2])
}

abundplot <- function(ref, sc1, sc2, sc3){
  par(las=1)
  plot(out(ref)[c("times", "number")],
      type="l", xlab="Time", ylab="Abundance",
      log="y", ylim=c(10, 5000),
      main="",
      font=2, font.lab=2, lwd=2
  )
  lines(out(sc1)[c("times", "number")], col="blue", lty="6222", lwd=2)
  lines(out(sc2)[c("times", "number")], col="darkgreen", lty="22", lwd=2)
  lines(out(sc3)[c("times", "number")], col="red", lty="66", lwd=2)
  legend("topleft",
    c(paste("Sc0: no mortality,     r =", round(growthrate(ref),2), collapse=" "),
      paste("Sc1: random mort.   r =", round(growthrate(sc1),2), collapse=" "),
      paste("Sc2: small pref.        r =", round(growthrate(sc2),2), collapse=" "),
      paste("Sc3: large pref.         r =", round(growthrate(sc3),2), collapse=" ")
      ),
    col=c("black","blue","darkgreen","red"), lty=c("solid", "6222", "22", "66"),
    bty="n", lwd=2
  )
}
abundplot(sc$Sc0, sc$Sc1, sc$Sc2, sc$Sc3)


###################################################
### code chunk number 30: cleanup
###################################################
options("prompt" = "> ", "continue" = "+ ")



###################################################
### Initialization and graphics settings
###################################################
opar   <- par(no.readonly=TRUE)
oask   <- devAskNewPage(dev.interactive(orNone = TRUE))
defpar <- par(no.readonly = TRUE)

library("simecol")
pkgdir <- system.file(package="simecol")

###################################################
### First example
###################################################

data(lv, package = "simecol")
lv <- sim(lv)
plot(lv)
o1 <- out(lv)


###################################################
### Parameter manipulation
###################################################
parms(lv)
parms(lv)["k1"] <- 0.4

parms(lv)["a"] <- 1
parms(lv)

parms(lv) <- parms(lv)[-4]
## or even more general with:
#parms(lv) <- parms(lv)[-which(names(parms(lv)) == "a")]
parms(lv)

###################################################
readline("Press Return to Continue ")
par(opar)
devAskNewPage(FALSE) # enable animated graphics

###################################################
### Simple cellular automaton
###################################################

source(paste(pkgdir, "/doc/examples/", "CA.R", sep=""))
times(CA)
times(CA) <- c(to=100)
CA <- sim(CA)
plot(CA)

par(defpar) # interactive on
###################################################
### Figure 3
###################################################

set.seed(345)
times(CA) <- c(to=50)
CA <- sim(CA)

library(lattice)
tcol <- (terrain.colors(13))[-13]
x <- out(CA, last=TRUE)
x <- ifelse(x == 0, NA, x)
print(levelplot(x,
             cuts = 11,
             col.regions = tcol,
             colorkey=list(at=seq(0,55,5))
))

###################################################
### Two L&V Scenarios
###################################################
data(lv, package="simecol")
plot(sim(lv))

lv1 <- lv2 <- lv
parms(lv1)
init(lv1)
parms(lv2) <- c(k1=0.2, k2=0.2, k3=0.1)
lv1 <- sim(lv1)
lv2 <- sim(lv2)


###################################################
### Figure 4
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

par(defpar)

###################################################
### Model result statistics example
###################################################
sapply(o1[c("predator", "prey")], range)
sapply(o2[c("predator", "prey")], range)

tlv <- times(lv1)
ots        <- ts(o1[c("predator", "prey")], start=tlv["from"],
                 end=tlv["to"], deltat=tlv["by"])
sp         <- spectrum(ots, spans=c(3,3), log="no")
1/sp$freq[sp$spec[,1] == max(sp$spec[,1])]
Sys.sleep(2)
acf(ots, lag.max=100)

###################################################
### Non-autonomous system (with inputs)
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

inputs(lv_ef) <-  as.matrix(data.frame(
  time = c(0, 30, 30.1, 100),
  s.in = c(0,  0,  .5,     .5)
))

###################################################
### Figure 5
###################################################
o <- out(sim(lv_ef))
matplot(o$time,o[2:3], xlab="Time",
  ylab="Substrate, Prey, Predator", type="l",
  lty=c("solid", "dashed"), col=c("blue", "red"), las=1)
inp <- as.data.frame(inputs(lv_ef))
lines(inp$time, inp$s.in, col="darkgreen", lwd=2, lty="dotted")


###################################################
### Demo of initfunc (Monte-Carlo example)
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

###################################################
### Figure 6
###################################################
matplot(tt, o, xlab="time", ylab="Prey", las=1, type="l")


###################################################
### Daphnia example
### please see details in appendix
### or external file
### Figure 7
###################################################
source(paste(pkgdir, "/doc/examples/", "ibm_daphnia.R", sep=""))

## Show the contents of the example file
#edit(file=paste(pkgdir, "/examples/", "ibm_daphnia.R", sep=""))

solver(ibm_daphnia) <- "myiteration"
ibm_daphnia <- sim(ibm_daphnia)
plot(ibm_daphnia)


###################################################
### Selective predation example
###################################################
Sc0                                   <- ibm_daphnia
times(Sc0)                            <- c(from=0, to=30, by=0.2)
parms(Sc0)[c("temp", "food", "mort")] <- c(15, 0.4, 0.1)
init(Sc0) <- data.frame(age=rep(10, 50), size = rep(2.5, 50),
                        eggs=rep(5, 50), eggage=runif(50, 0, 4))

equations(Sc0)$survive = function(inds, parms) {
  abundance <- nrow(inds)
  rnd       <- runif(abundance)
  mort      <- fmort(parms$mort, inds$size) * parms$DELTAT
  subset(inds, rnd > mort)
}

Sc1 <- Sc2 <- Sc3 <- Sc0
equations(Sc0)$fmort <- function(mort, x) 0
equations(Sc1)$fmort <- function(mort, x) mort
equations(Sc2)$fmort <- function(mort, x){
  mort * 2 * rank(-x) / (length(x) + 1)
}
equations(Sc3)$fmort <- function(mort, x){
  mort * 2 * rank(x) / (length(x) + 1)
}

set.seed(123)
sc <- lapply(list(Sc0=Sc0, Sc1=Sc1, Sc2=Sc2, Sc3=Sc3), sim)


###################################################
### Figure 8
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

## clean up
par(defpar)
devAskNewPage(oask)

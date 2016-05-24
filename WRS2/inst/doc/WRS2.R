## ----echo=FALSE, message=FALSE, results='hide'---------------------------
require(WRS2)
require(beanplot)
require(car)
require(MASS)
require(colorspace)
require(reshape)
require(ez)

## ----soccer-plot, eval=FALSE, echo = FALSE-------------------------------
#  SpainGer <- subset(eurosoccer, League == "Spain" | League == "Germany")
#  SpainGer <- droplevels(SpainGer)
#  op <- par(mfrow = c(1,2))
#  boxplot(GoalsGame ~ League, data = SpainGer, main = "Boxplot Goals Scored per Game")
#  points(1:2, tapply(SpainGer$GoalsGame, SpainGer$League, mean, trim = 0.2), pch = 19, col = "red")
#  beanplot(GoalsGame ~ League, data = SpainGer, log = "", main = "Beanplot Goals Scored per Game",
#           col = "coral")
#  par(op)

## ----soccer-plot1, echo=FALSE, fig.height = 6, fig.width = 12, dev='postscript'----
SpainGer <- subset(eurosoccer, League == "Spain" | League == "Germany")
SpainGer <- droplevels(SpainGer)
op <- par(mfrow = c(1,2))
boxplot(GoalsGame ~ League, data = SpainGer, main = "Boxplot Goals Scored per Game")
points(1:2, tapply(SpainGer$GoalsGame, SpainGer$League, mean, trim = 0.2), pch = 19, col = "red")
beanplot(GoalsGame ~ League, data = SpainGer, log = "", main = "Beanplot Goals Scored per Game", 
         col = "coral")
par(op)

## ------------------------------------------------------------------------
yuen(GoalsGame ~ League, data = SpainGer)

## ------------------------------------------------------------------------
pb2gen(GoalsGame ~ League, data = SpainGer, est = "median")
pb2gen(GoalsGame ~ League, data = SpainGer, est = "onestep")

## ----soccer2-plot, eval=FALSE, echo=FALSE--------------------------------
#  op <- par(mfrow = c(2,1))
#  boxplot(GoalsGame ~ League, data = eurosoccer, main = "Boxplot Goals Scored per Game")
#  beanplot(GoalsGame ~ League, data = eurosoccer, log = "", main = "Beanplot Goals Scored per Game",
#           col = "coral")
#  par(op)

## ----soccer2-plot1, echo=FALSE, fig.height = 10, fig.width = 12, dev='postscript'----
op <- par(mfrow = c(2,1))
boxplot(GoalsGame ~ League, data = eurosoccer, main = "Boxplot Goals Scored per Game")
beanplot(GoalsGame ~ League, data = eurosoccer, log = "", main = "Beanplot Goals Scored per Game", 
         col = "coral")
par(op)

## ------------------------------------------------------------------------
t1way(GoalsGame ~ League, data = eurosoccer)
med1way(GoalsGame ~ League, data = eurosoccer)

## ------------------------------------------------------------------------
lincon(GoalsGame ~ League, data = eurosoccer)

## ----goggles-plot, eval=FALSE, echo=FALSE--------------------------------
#  attach(goggles)
#  op <- par(mfrow = c(1,2))
#  interaction.plot(gender, alcohol, attractiveness, fun = median, ylab = "Attractiveness", xlab = "Gender", type = "b", pch = 20,
#                   lty = 1, col = c("black", "cadetblue", "coral"), main = "Interaction Plot Alcohol|Gender", legend = FALSE)
#  legend("right", legend = c("None", "2 Pints","4 Pints"), col = c("black", "cadetblue", "coral"), lty = 1, cex = 0.8)
#  interaction.plot(alcohol, gender, attractiveness, fun = median, ylab = "Attractiveness", xlab = "Gender", type = "b", pch = 20,
#                   lty = 1, col = c("coral", "black"), main = "Interaction Plot Gender|Alcohol", legend = FALSE)
#  legend("bottomleft", legend = c("female", "male"), col = c("coral", "black"), lty = 1, cex = 0.8)
#  par(op)
#  detach(goggles)

## ----goggles-plot1, echo=FALSE, fig.height = 6, fig.width = 12, dev='postscript'----
attach(goggles)
op <- par(mfrow = c(1,2))
interaction.plot(gender, alcohol, attractiveness, fun = median, ylab = "Attractiveness", xlab = "Gender", type = "b", pch = 20,
                 lty = 1, col = c("black", "cadetblue", "coral"), main = "Interaction Plot Alcohol|Gender", legend = FALSE)  
legend("right", legend = c("None", "2 Pints","4 Pints"), col = c("black", "cadetblue", "coral"), lty = 1, cex = 0.8)
interaction.plot(alcohol, gender, attractiveness, fun = median, ylab = "Attractiveness", xlab = "Gender", type = "b", pch = 20,
                 lty = 1, col = c("coral", "black"), main = "Interaction Plot Gender|Alcohol", legend = FALSE) 
legend("bottomleft", legend = c("female", "male"), col = c("coral", "black"), lty = 1, cex = 0.8)
par(op)
detach(goggles)

## ----echo=3:6------------------------------------------------------------
set.seed(123)
goggles$alcohol <- relevel(goggles$alcohol, ref = "None") #FIXME
t2way(attractiveness ~ gender*alcohol, data = goggles)
med2way(attractiveness ~ gender*alcohol, data = goggles)
pbad2way(attractiveness ~ gender*alcohol, data = goggles, est = "onestep")
summary(aov(attractiveness ~ gender*alcohol, data = goggles))

## ------------------------------------------------------------------------
mcp2atm(attractiveness ~ gender*alcohol, data = goggles)

## ----swim-plot, echo=FALSE, eval=FALSE-----------------------------------
#  tmean20 <- function(x) mean(x, trim = 0.20)
#  optpes.male <- subset(swimming, Sex == "Male")
#  optpes.female <- subset(swimming, Sex == "Female")
#  op <- par(mfrow = c(1,2))
#  interaction.plot(optpes.male$Event, optpes.male$Optim, optpes.male$Ratio, fun = tmean20,
#                   xlab = "Event", ylab = "Ratio", main = "Interaction Men",
#                   type = "b", pch = 20, lty = 1, col = 1:2, legend = FALSE)
#  legend("topleft", legend = c("Optimists", "Pessimists"), col = 1:2, lty = 1, cex = 0.8)
#  interaction.plot(optpes.female$Event, optpes.female$Optim, optpes.female$Ratio, fun = tmean20,
#                   xlab = "Event", ylab = "Ratio", main = "Interaction Women",
#                   type = "b", pch = 20, lty = 1, col = 1:2, legend = FALSE)
#  legend("topleft", legend = c("Optimists", "Pessimists"), col = 1:2, lty = 1, cex = 0.8)
#  par(op)

## ----swim-plot1, echo=FALSE, fig.height = 6, fig.width = 12, dev='postscript'----
tmean20 <- function(x) mean(x, trim = 0.20)
optpes.male <- subset(swimming, Sex == "Male")
optpes.female <- subset(swimming, Sex == "Female")
op <- par(mfrow = c(1,2))
interaction.plot(optpes.male$Event, optpes.male$Optim, optpes.male$Ratio, fun = tmean20, 
                 xlab = "Event", ylab = "Ratio", main = "Interaction Men", 
                 type = "b", pch = 20, lty = 1, col = 1:2, legend = FALSE)
legend("topleft", legend = c("Optimists", "Pessimists"), col = 1:2, lty = 1, cex = 0.8)
interaction.plot(optpes.female$Event, optpes.female$Optim, optpes.female$Ratio, fun = tmean20,
                 xlab = "Event", ylab = "Ratio", main = "Interaction Women", 
                 type = "b", pch = 20, lty = 1, col = 1:2, legend = FALSE)
legend("topleft", legend = c("Optimists", "Pessimists"), col = 1:2, lty = 1, cex = 0.8)
par(op)

## ------------------------------------------------------------------------
t3way(Ratio ~ Optim*Sex*Event, data = swimming)
fitaov_op <- aov(Ratio ~ Optim*Sex*Event, data = swimming)  
Anova(fitaov_op, type = "II")       

## ----swim2-plot, echo=FALSE, eval=FALSE----------------------------------
#  op <- par(mfrow = c(1,2))
#  interaction.plot(swimming$Optim, swimming$Sex, swimming$Ratio, fun = tmean20,
#                   xlab = "Optimist/Pessimist", ylab = "Ratio", main = "Interaction Sex|Optim",
#                   type = "b", pch = 19, lty = 1, col = 1:2, legend = FALSE)
#  legend("bottomleft", legend = c("Male", "Female"), col = 1:2, lty = 1, cex = 0.8)
#  interaction.plot(swimming$Sex, swimming$Optim, swimming$Ratio, fun = tmean20,
#                   xlab = "Sex", ylab = "Ratio", main = "Interaction Optim|Sex",
#                   type = "b", pch = 19, lty = 1, col = 1:2, legend = FALSE)
#  legend("bottomright", legend = c("Optimists", "Pessimists"), col = 1:2, lty = 1, cex = 0.8)
#  par(op)

## ----swim2-plot1, echo=FALSE, fig.height = 6, fig.width = 12, dev='postscript'----
op <- par(mfrow = c(1,2))
interaction.plot(swimming$Optim, swimming$Sex, swimming$Ratio, fun = tmean20, 
                 xlab = "Optimist/Pessimist", ylab = "Ratio", main = "Interaction Sex|Optim", 
                 type = "b", pch = 19, lty = 1, col = 1:2, legend = FALSE)
legend("bottomleft", legend = c("Male", "Female"), col = 1:2, lty = 1, cex = 0.8)
interaction.plot(swimming$Sex, swimming$Optim, swimming$Ratio, fun = tmean20, 
                 xlab = "Sex", ylab = "Ratio", main = "Interaction Optim|Sex", 
                 type = "b", pch = 19, lty = 1, col = 1:2, legend = FALSE)
legend("bottomright", legend = c("Optimists", "Pessimists"), col = 1:2, lty = 1, cex = 0.8)
par(op)

## ------------------------------------------------------------------------
anorexiaFT <- subset(anorexia, subset = Treat == "FT")
yuend(anorexiaFT$Prewt, anorexiaFT$Postwt)

## ----ano-plot, echo=FALSE, eval=FALSE------------------------------------
#  colpal <- c(rainbow_hcl(17, c = 50))
#  matplot(t(anorexiaFT[,2:3]), type = "b", cex = 0.8, main = "Weight Trajectories",
#          xaxt = "n", ylab = "Weight (lbs.)", lty = 1, col = colpal, pch = 20)
#  axis(1, at = 1:2, labels = c("Prior", "Post"))

## ----ano-plot1, echo=FALSE, fig.height = 5, fig.width = 8, dev='postscript'----
colpal <- c(rainbow_hcl(17, c = 50))
matplot(t(anorexiaFT[,2:3]), type = "b", cex = 0.8, main = "Weight Trajectories", 
        xaxt = "n", ylab = "Weight (lbs.)", lty = 1, col = colpal, pch = 20)
axis(1, at = 1:2, labels = c("Prior", "Post"))

## ----wine-plot, echo=FALSE, eval=FALSE-----------------------------------
#  WineTasting_wide <- reshape(WineTasting, idvar = "Taster", timevar = "Wine", direction = "wide")[-1]  ## wide format
#  colpal <- c(rainbow_hcl(22, c = 100))
#  #pal <- palette(colpal)
#  matplot(t(WineTasting_wide), pch = 20, type = "b", xaxt = "n", xlab = "Wines", ylab = "Score", lty = 1, col = colpal, main = "Wine Trajectories")
#  axis(1, at = 1:3, labels = levels(WineTasting$Wine))
#  #palette(pal)

## ----wine-plot1, echo=FALSE, fig.height = 5, fig.width = 8, dev='postscript'----
WineTasting_wide <- reshape(WineTasting, idvar = "Taster", timevar = "Wine", direction = "wide")[-1]  ## wide format
colpal <- c(rainbow_hcl(22, c = 100))
#pal <- palette(colpal)
matplot(t(WineTasting_wide), pch = 20, type = "b", xaxt = "n", xlab = "Wines", ylab = "Score", lty = 1, col = colpal, main = "Wine Trajectories")
axis(1, at = 1:3, labels = levels(WineTasting$Wine))
#palette(pal)

## ----echo=2:3------------------------------------------------------------
attach(WineTasting)
rmanova(y = Taste, groups = Wine, block = Taster)
rmmcp(y = Taste, groups = Wine, block = Taster)
detach(WineTasting)

## ----hang-plot, eval=FALSE, echo=FALSE-----------------------------------
#  ind <- rep(1:6, each = 20)
#  symlist <- split(hangover$symptoms, ind)[c(1,4,2,5,3,6)]
#  gtmeans <- sapply(symlist, mean, trim = 0.2)
#  plot(1:3, type = "n", ylim = c(0, max(hangover$symptoms) + 10), xaxt = "n", xlab = "Time Points",
#       ylab = "Number of Symptoms", main = "Hangover Data")
#  axis(1, at = 1:3, labels = paste("Time", 1:3))
#  for (i in 1:6) points(jitter(rep(ceiling(i/2), 20)), symlist[[i]], cex = 0.6, col = ((i %% 2) + 1))
#  legend("topleft", legend = c("control", "alcoholic"), lty = 1, col = 1:2)
#  lines(1:3, gtmeans[c(1, 3, 5)], col = 1, type = "b", pch = 19)
#  lines(1:3, gtmeans[c(2, 4, 6)], col = 2, type = "b", pch = 19)

## ----hang-plot1, echo=FALSE, fig.height = 5, fig.width = 8, dev='postscript'----
ind <- rep(1:6, each = 20)
symlist <- split(hangover$symptoms, ind)[c(1,4,2,5,3,6)]
gtmeans <- sapply(symlist, mean, trim = 0.2)
plot(1:3, type = "n", ylim = c(0, max(hangover$symptoms) + 10), xaxt = "n", xlab = "Time Points", 
     ylab = "Number of Symptoms", main = "Hangover Data")
axis(1, at = 1:3, labels = paste("Time", 1:3))
for (i in 1:6) points(jitter(rep(ceiling(i/2), 20)), symlist[[i]], cex = 0.6, col = ((i %% 2) + 1))
legend("topleft", legend = c("control", "alcoholic"), lty = 1, col = 1:2)
lines(1:3, gtmeans[c(1, 3, 5)], col = 1, type = "b", pch = 19)
lines(1:3, gtmeans[c(2, 4, 6)], col = 2, type = "b", pch = 19)

## ------------------------------------------------------------------------
bwtrim(symptoms ~ group*time, id = id, data = hangover)

## ----warning=FALSE-------------------------------------------------------
bwtrim(symptoms ~ group*time, id = id, data = hangover, tr = 0)
fitF <- ezANOVA(hangover, symptoms, between = group, within = time, wid = id)
fitF$ANOVA

## ----echo=2:4------------------------------------------------------------
set.seed(123)
sppba(symptoms ~ group*time, id, data = hangover)
sppbb(symptoms ~ group*time, id, data = hangover)
sppbi(symptoms ~ group*time, id, data = hangover)

## ----smooth-plot, eval=FALSE, echo=FALSE---------------------------------
#  colpal <- c(rainbow_hcl(5, c = 100))
#  pal <- palette(colpal)
#  attach(chile)
#  op <- par(mfrow = c(1,2))
#  plot(length, heat, pch = 20, col = "gray", main = "Chile Smoothing I", xlab = "Length", ylab = "Heat")
#  fitmean <- runmean(length, heat)
#  fitmest <- rungen(length, heat)
#  fitmed <- rungen(length, heat, est = "median")
#  fitbag <- runmbo(length, heat, est = "onestep")
#  orderx <- order(length)
#  lines(length[orderx], fitmean[orderx], lwd = 2, col = 1)
#  lines(length[orderx], fitmest[orderx], lwd = 2, col = 2)
#  lines(length[orderx], fitmed[orderx], lwd = 2, col = 3)
#  lines(length[orderx], fitbag[orderx], lwd = 2, col = 4)
#  legend("topright", legend = c("Trimmed Mean", "MOM", "Median", "Bagged Onestep"), col = 1:4, lty = 1)
#  plot(length, heat, pch = 20, col = "gray", main = "Chile Smoothing II", xlab = "Length", ylab = "Heat")
#  fitmean1 <- runmean(length, heat, fr = 0.2)
#  fitmean2 <- runmean(length, heat, fr = 0.5)
#  fitmean3 <- runmean(length, heat, fr = 1)
#  fitmean4 <- runmean(length, heat, fr = 5)
#  orderx <- order(length)
#  lines(length[orderx], fitmean1[orderx], lwd = 2, col = 1)
#  lines(length[orderx], fitmean2[orderx], lwd = 2, col = 2)
#  lines(length[orderx], fitmean3[orderx], lwd = 2, col = 3)
#  lines(length[orderx], fitmean4[orderx], lwd = 2, col = 4)
#  legend("topright", legend = c("f = 0.2", "f = 0.5", "f = 1", "f = 5"), col = 1:4, lty = 1)
#  par(op)
#  palette(pal)
#  detach(chile)

## ----smooth-plot1, echo=FALSE, fig.height = 6, fig.width = 12, dev='postscript'----
colpal <- c(rainbow_hcl(5, c = 100))
pal <- palette(colpal)
attach(chile)
op <- par(mfrow = c(1,2))
plot(length, heat, pch = 20, col = "gray", main = "Chile Smoothing I", xlab = "Length", ylab = "Heat")
fitmean <- runmean(length, heat)
fitmest <- rungen(length, heat)
fitmed <- rungen(length, heat, est = "median")
fitbag <- runmbo(length, heat, est = "onestep")
orderx <- order(length)
lines(length[orderx], fitmean[orderx], lwd = 2, col = 1)
lines(length[orderx], fitmest[orderx], lwd = 2, col = 2)
lines(length[orderx], fitmed[orderx], lwd = 2, col = 3)
lines(length[orderx], fitbag[orderx], lwd = 2, col = 4)
legend("topright", legend = c("Trimmed Mean", "MOM", "Median", "Bagged Onestep"), col = 1:4, lty = 1)
plot(length, heat, pch = 20, col = "gray", main = "Chile Smoothing II", xlab = "Length", ylab = "Heat")
fitmean1 <- runmean(length, heat, fr = 0.2)
fitmean2 <- runmean(length, heat, fr = 0.5)
fitmean3 <- runmean(length, heat, fr = 1)
fitmean4 <- runmean(length, heat, fr = 5)
orderx <- order(length)
lines(length[orderx], fitmean1[orderx], lwd = 2, col = 1)
lines(length[orderx], fitmean2[orderx], lwd = 2, col = 2)
lines(length[orderx], fitmean3[orderx], lwd = 2, col = 3)
lines(length[orderx], fitmean4[orderx], lwd = 2, col = 4)
legend("topright", legend = c("f = 0.2", "f = 0.5", "f = 1", "f = 5"), col = 1:4, lty = 1)
par(op)
palette(pal)
detach(chile)

## ----echo=2:5------------------------------------------------------------
comppts <- c(18, 70, 80, 90, 100, 110)
fitanc <- ancova(Posttest ~ Pretest + Group, fr1 = 0.3, fr2 = 0.3, 
                 data = electric, pts = comppts)
fitanc

## ----anc-plot, eval=FALSE, echo=FALSE------------------------------------
#  plot(electric$Pretest, electric$Posttest, col = rep(1:2, each = 96), pch = 1, cex = 0.8,
#        xlab = "Pretest Score", ylab = "Posttest Score", main = "TV Show ANCOVA")
#  eltr <- subset(electric, subset = Group == "treatment")
#  elct <- subset(electric, subset = Group == "control")
#  ordtr <- order(eltr$Pretest)
#  lines(eltr$Pretest[ordtr], fitanc$fitted.values$treatment[ordtr], col = 1, lwd = 2)
#  abline(lm(eltr$Posttest ~ eltr$Pretest), col = 1, lty = 2)
#  ordct <- order(elct$Pretest)
#  lines(elct$Pretest[ordct], fitanc$fitted.values$control[ordct], col = 2, lwd = 2)
#  abline(lm(elct$Posttest ~ elct$Pretest), col = 2, lty = 2)
#  abline(v = comppts, lty = 2, col = "gray")
#  legend(30, 120, legend = c("treatment", "control"), lty = 1, col = 1:2)

## ----anc-plot1, echo=FALSE, fig.height = 6, fig.width = 9, dev='postscript'----
plot(electric$Pretest, electric$Posttest, col = rep(1:2, each = 96), pch = 1, cex = 0.8, 
      xlab = "Pretest Score", ylab = "Posttest Score", main = "TV Show ANCOVA")
eltr <- subset(electric, subset = Group == "treatment")
elct <- subset(electric, subset = Group == "control")
ordtr <- order(eltr$Pretest)
lines(eltr$Pretest[ordtr], fitanc$fitted.values$treatment[ordtr], col = 1, lwd = 2)
abline(lm(eltr$Posttest ~ eltr$Pretest), col = 1, lty = 2)
ordct <- order(elct$Pretest)
lines(elct$Pretest[ordct], fitanc$fitted.values$control[ordct], col = 2, lwd = 2)
abline(lm(elct$Posttest ~ elct$Pretest), col = 2, lty = 2)
abline(v = comppts, lty = 2, col = "gray")
legend(30, 120, legend = c("treatment", "control"), lty = 1, col = 1:2)


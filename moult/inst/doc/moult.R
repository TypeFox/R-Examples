### R code from vignette source 'moult.Rnw'

###################################################
### code chunk number 1: preliminaries
###################################################
options(prompt = "R> ", continue = "+  ", width = 70, useFancyQuotes = FALSE)


###################################################
### code chunk number 2: sanderling1
###################################################
library("moult")

head(sanderlings)

plot(MIndex ~ Day, data = sanderlings, pch = 16, cex = 0.5, 
     xlab = "days since July 1",
     ylab = "proportion of feather mass grown")


###################################################
### code chunk number 3: sanderling2
###################################################
m1 <- moult(MIndex ~ Day, data = sanderlings, type = 1)        

m2 <- moult(MIndex ~ Day, data = sanderlings)                  
summary(m2)

m3 <- moult(MIndex ~ Day, data = sanderlings, type = 3)             
summary(m3)


###################################################
### code chunk number 4: sanderling4
###################################################
format(as.Date(coef(m2, "mean"), origin = "2012-06-30"), "%d %b")


###################################################
### code chunk number 5: moult.Rnw:242-249
###################################################
durationmean2ab <- function(duration, mean)
  { ab <- c(- mean / duration, 1 / duration)
    names(ab) <- c("intercept", "slope")
    return(ab)
  } 

uz1 <- durationmean2ab(coef(m1, "duration"), coef(m1, "mean"))


###################################################
### code chunk number 6: sanderling.predict
###################################################
day <- unique(sanderlings$Day)
p1 <- predict(m2, newdata = data.frame(day))    
nn <- as.numeric(table(sanderlings$Day))
p1$M * nn                   


###################################################
### code chunk number 7: moult.Rnw:267-270 (eval = FALSE)
###################################################
## days2 <- seq(70, 310, by = 10)
## p2 <- predict(m2, newdata = data.frame(days2))
## p2$M * 100                   


###################################################
### code chunk number 8: probit
###################################################
sanderlings$started <- ifelse (sanderlings$MIndex > 0, 1, 0)
m.start <- glm(started ~ Day, family = binomial(link = "probit"),
               data = sanderlings)
summary(m.start)


###################################################
### code chunk number 9: moult.Rnw:291-296
###################################################
cfs <- as.numeric(coef(m.start))
sigma.hat1 <- 1 / cfs[2]
mu.hat1 <- - cfs[1] / cfs[2]
c(mu.hat1, sigma.hat1)
format(as.Date(mu.hat1, origin = "2012-06-30"), "%d %b")


###################################################
### code chunk number 10: moult.Rnw:301-312
###################################################
sanderlings$finished <- ifelse (sanderlings$MIndex < 1, 0, 1)
m.end <- glm(finished ~ Day, family = binomial(link = "probit"),
             data = sanderlings)

cfs <- as.numeric(coef(m.end))
sigma.hat2 <- 1 / cfs[2]
mu.hat2 <- - cfs[1] / cfs[2]
duration <- mu.hat2 - mu.hat1
c(duration, mu.hat2, sigma.hat2)
format(as.Date(mu.hat2, origin = "2012-06-30"), "%d %b")
nr1 <- durationmean2ab(duration, mu.hat1)


###################################################
### code chunk number 11: probit3
###################################################
scores <- c(sanderlings$started, sanderlings$finished)
x <- c(rep(0, times = length(sanderlings$started)), 
       rep(1, times = length(sanderlings$finished)))
ddays <- rep(sanderlings$Day, times = 2)

m.both <- glm(scores ~ ddays + x, family = binomial(link = "probit"))
summary(m.both)

cfs <- as.numeric(coef(m.both))
sigma.hat <- 1 / cfs[2]
mu.hat <- - cfs[1] / cfs[2]
mean.finish <- - (cfs[1] + cfs[3]) / cfs[2]
duration <- - cfs[3] / cfs[2]
c(duration, mu.hat, mean.finish, sigma.hat)
format(as.Date(c(mu.hat, mean.finish), origin = "2012-06-30"), "%d %b")

nr2 <- durationmean2ab(duration, mu.hat)


###################################################
### code chunk number 12: probitplots
###################################################
plot(MIndex ~ Day, data = sanderlings, pch = 16, cex = 0.5, 
     xlab = "days since July 1", ylab = "moult index (PFMG)")

abline(uz1, lwd = 2, col = "grey")             
abline(nr1, lty = 2, lwd = 2)                 
abline(nr2, lty = 3, col = "red", lwd = 3)    
legend(220, 0.3, lty = c(1, 2, 3), lwd = c(2, 2, 3), 
       col = c("grey", "black", "red"), bty = "n", 
       legend = c("UZ type 1", "NR separate", "NR combined"))


###################################################
### code chunk number 13: moult.Rnw:360-367
###################################################
rm(list = ls())

durationmean2ab <- function(duration, mean)
  { ab <- c(- mean / duration, 1 / duration)
    names(ab) <- c("intercept", "slope")
    return(ab)
  } 


###################################################
### code chunk number 14: weavers1
###################################################
data("weavers")
head(weavers)

if (is.numeric(weavers$Moult)) {
  scores <- format(weavers$Moult, scientific = FALSE, trim = TRUE) 
} else {
  scores <- weavers$Moult 
}

mscores <- substr(scores, 1, 9)
  
feather.mass <- c(10.4, 10.8, 11.5, 12.8, 14.4, 15.6, 16.3, 15.7, 15.7) 

weavers$pfmg <- ms2pfmg(mscores, feather.mass)  
weavers$day <- date2days(weavers$RDate, dateformat = "yyyy-mm-dd", 
                         startmonth = 8)  


###################################################
### code chunk number 15: moult.Rnw:407-414
###################################################
m88.2 <- moult(pfmg ~ day, data = weavers)
summary(m88.2)

m88c <- moult(pfmg ~ day, data = weavers, type = 3)
summary(m88c)

uz1 <- durationmean2ab(coef(m88c, "duration"), coef(m88c, "mean"))


###################################################
### code chunk number 16: weavers2
###################################################
ssex <- ifelse(weavers$Sex == 1 | weavers$Sex == 3, 'male', 
        ifelse(weavers$Sex == 2 | weavers$Sex == 4, 'female', NA))

weavers$ssex <- as.factor(ssex)

mmf <- moult(pfmg ~ day | ssex | ssex, data = weavers, type = 3)
summary(mmf)


###################################################
### code chunk number 17: weaverplots2
###################################################
fstart <- coef(mmf, "mean")[1]
fduration <- coef(mmf, "duration")[1] 

mstart <- sum(coef(mmf, "mean")[1:2])
mduration <- sum(coef(mmf, "duration")[1:2])

female.traj <- durationmean2ab(fduration, fstart)
male.traj <- durationmean2ab(mduration, mstart)

plot(pfmg ~ day, pch = 16, cex = 0.5, xlab = "days since August 1",
     ylab = "proportion of feather mass grown", las = 1, col = "grey", 
     data = weavers)
abline(uz1, lwd = 2)
abline(male.traj, col = "blue", lty = 2, lwd = 3)
abline(female.traj, col = "red", lty = 4, lwd = 3)
legend(280, 0.2, lty = c(2, 4), lwd = 2, col = c("blue", "red"),
       legend = c("males", "females"), bty = "n")


###################################################
### code chunk number 18: weavers3
###################################################
weavers$year.f <- as.factor(weavers$Year)
m88y <- moult(pfmg ~ day | 1 | year.f, data = weavers, type = 3)
summary(m88y)

weavers$pfmg[weavers$Year == 1988]


###################################################
### code chunk number 19: weavers4
###################################################
weav89 <- weavers[weavers$Year >= 1989, ]
weav89$year.f <- as.factor(weav89$Year)
m89y <- moult(pfmg ~ day | 1 | year.f, data = weav89, type = 3)
summary(m89y)


###################################################
### code chunk number 20: moult.Rnw:490-492
###################################################
pred.year <- predict(m89y, predict.type = "start", 
                     newdata = data.frame(year.f = as.factor(1989:2005)))


###################################################
### code chunk number 21: weaversrain2
###################################################
rainSep <- c(44.8, 110.2, 24.3, 84.4, 73.1, 10.5, 32.4, 3.8, 91.2,
               8.1, 27.8, 112.4, 58.7, 111.5, 20, 66.3, 43, 14.5)
weavers$rain <- rainSep[match(weavers$Year, 1988:2005)]

m88r <- moult(pfmg ~ day | 1 | rain, data = weavers, type = 3)
m88r2 <- moult(pfmg ~ day | 1 | rain + I(rain^2), data = weavers, type = 3)
summary(m88r2)


###################################################
### code chunk number 22: AIC
###################################################
AIC(m88y, m88c, m88r, m88r2)               


###################################################
### code chunk number 23: weaverplots1
###################################################
par(mar = c(7, 9, 1, 1), mgp = c(5, 1, 0))
plot(1988:2005, c(NA, pred.year[, 1]), type = "b", lwd = 1, pch = 19, 
     xlab = "year", ylab = "mean start of moult", 
     ylim = c(70, 210), las = 1, yaxt = "n")
abline(h = coef(m88c, "mean"), col = "grey")
upp <- pred.year[, 1] + 1.96 * pred.year[, 2]  
lwr <- pred.year[, 1] - 1.96 * pred.year[, 2]  
arrows(x0 = 1989:2005, y0 = lwr, y1 = upp, code = 3, angle = 90, 
       length = 0.05)
ylab <- format(as.Date(seq(70, 210, by = 20), origin = "2012-07-31"), "%d %b")
axis(2, at = seq(70, 210, by = 20), labels = ylab, las = 1)


###################################################
### code chunk number 24: moult.Rnw:540-541
###################################################
rm(list = ls())



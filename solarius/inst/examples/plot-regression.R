### inc
library(ggplot2)
library(RColorBrewer)
library(gridExtra)

### models
mod <- solarPolygenic(trait2 ~ age + age^2 + sex, dat30)

dat <- within(dat30, {
  sex = factor(sex, labels = c("M", "F"))
})
lmod <- lm(trait2 ~ age + age^2 + sex, dat)

### simulated data
N <- 100

age.range <- range(dat$age, na.rm = TRUE)
simdat <- data.frame(trait = 1:N, age = seq(age.range[1], age.range[2], length = N),
  sex = factor(rep("M", N), levels = c("M", "F")))
  
simdat2 <- dat
simdat2$sex <- "F"

### predict
pdat <- predict(lmod, simdat,
  interval = "confidence", level = 0.95)

pdat <- data.frame(pdat, age = simdat$age)

pdat2 <- predict(lmod, simdat2,
  interval = "confidence", level = 0.95)
pdat2 <- data.frame(pdat2, age = simdat2$age)

### plot
psize <- 2
cols <- brewer.pal(3, "Set1")

p1 <- ggplot() +
  geom_point(data = subset(dat, sex == "M"), aes(age, trait2), col = cols[2], size = psize) +
  geom_point(data = subset(dat, sex == "F"), aes(age, trait2), col = cols[1], size = psize) +  
  # Male
  geom_line(data = pdat, aes(age, fit), col = cols[2]) +
  geom_ribbon(data = pdat, aes(x = age, ymin = lwr, ymax = upr), fill = cols[2], alpha = 0.5) +
  geom_line(data = pdat, aes(x = age, y = upr), col = cols[2], linetype = 2) +
  geom_line(data = pdat, aes(x = age, y = lwr), col = cols[2], linetype = 2) +
  # Female
  geom_line(data = pdat2, aes(age, fit), col = cols[1]) +
  geom_ribbon(data = pdat2, aes(x = age, ymin = lwr, ymax = upr), fill = cols[1], alpha = 0.5) +
  geom_line(data = pdat2, aes(x = age, y = upr), col = cols[1], linetype = 2) +
  geom_line(data = pdat2, aes(x = age, y = lwr), col = cols[1], linetype = 2) +
  labs(x = "age", y = "trait2")


### JetLagKnees

# Approximate Figure 15.1-1
xyplot(shift ~ treatment, data = JetLagKnees)
bwplot(shift ~ treatment, data = JetLagKnees)

# Table 15.1-1
if (require(plyr)){
  smry <- ddply(JetLagKnees, .(treatment),
                function(x)c(Mean = mean(x$shift),
                             s = var(x$shift),
                             n = length(x$shift)))
  print(smry)
  
  # Grand mean
  weighted.mean(smry$Mean, smry$n)
}

# Subset the three treatment groups
control <- subset(JetLagKnees, treatment == "control")$shift
knee <- subset(JetLagKnees, treatment == "knee")$shift
eyes <- subset(JetLagKnees, treatment == "eyes")$shift

# k is the number of groups
k <- length(unique(JetLagKnees$treatment))

# Calculate n
n <- length(JetLagKnees$shift)
control.n <- length(control)
knee.n <- length(knee)
eyes.n <- length(eyes)

# Calculate standard deviations
control.sd <- sd(control)
knee.sd <- sd(knee)
eyes.sd <- sd(eyes)

# Error mean square
(SS.error <- ((control.sd^2 * (control.n - 1)) +
  (knee.sd^2 * (knee.n - 1)) +
  (eyes.sd^2 * (eyes.n - 1))))
(MS.error <- SS.error / (n - k))

# Grand mean
(grand.mean <- (control.n * mean(control) + knee.n * mean(knee) + 
  eyes.n * mean(eyes)) / n)

# Group mean square
(SS.groups <- (control.n * (mean(control) - grand.mean)^2) +
  (knee.n * (mean(knee) - grand.mean)^2) +
  (eyes.n * (mean(eyes) - grand.mean)^2))
(MS.groups <- SS.groups / (k - 1))

# F
(F <- MS.groups / MS.error)

# P-value
pf(F, 2, 19, lower.tail = FALSE)

# Figure 15.1-3
(fcrit <- qf(0.05, 2, 19, lower.tail = FALSE))
curve(df(x, 2, 19), from = 0, to = 10,
      ylab = "Probability Density", 
      xlab = expression(F[paste("2,19")]),
      xaxs = "i", yaxs = "i")
x <- seq(fcrit, 10, length = 100)
y <- df(x, 2, 19)
polygon(c(x[1], x, x[100]), c(0, y, df(10, 2, 19)),
        col = "red", border = NA)

# R^2
(SS.total <- SS.groups + SS.error)
SS.groups/SS.total

# With aov()
aov.obj <- aov(shift ~ treatment, data = JetLagKnees)
summary(aov.obj)

library(seqSAFI)

fun <- function(x) x$x1 %*% c(seq(3, 0, , 333), rep(0, 667))/1000 - 
  apply(x$x2, 1, function(a) (1/30 * mean(a[300:1000]) + 3)^3) + 0.8 * sin(x$x3)

### step 1: 
### the initial designs are set up
### the functional inputs are immediately split up in the middle 
### x is not split further (scalar variable)

design <- createSafiDesign(method = "SB", d.f = 3, variable.names = c("g1", "g2", "x"))
design <- splitSafiDesign(s.d = design, new.split.points = list(0.5, 0.5, NULL))
plot(design)

x <- accessSafiDesign(s.d = design, n.timepoints = c(1000, 1000, 1))
y <- fun(x)

model <- safiModel(s.d = design, y = y)
plot(model)

### step 2: 
### g1: as on the first half shows influence, only this intervals is split up 
### g2: both intervals show influence, both are split up

design <- splitSafiDesign(s.d = design, new.split.points = list(c(1/6, 2/6), c(0.25, 0.75), NULL))
plot(design)

x <- accessSafiDesign(s.d = design, n.timepoints = c(1000, 1000, 1))
y <- fun(x)

model <- safiModel(s.d = design, y = y)
plot(model)

### step 3: 
### g1: the first and the second intervals seem active and are split once more 
### g2: the three last intervals seems active and are split once more

design <- splitSafiDesign(s.d = design, new.split.points = list(c(1/12, 3/12), c(3/8, 5/8, 7/8), NULL))
plot(design)

x <- accessSafiDesign(s.d = design, n.timepoints = c(1000, 1000, 1))
y <- fun(x)

model <- safiModel(s.d = design, y = y)
plot(model)


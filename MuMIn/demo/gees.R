###
# Example of model selection with GEE ranked by QIC
###

library(MuMIn)
require(geepack)
require(gee)
require(yags)
options(na.action = na.pass)

data(dietox, package = 'geepack')
dietox$Cu <- as.factor(dietox$Cu)

# Compare GEE fits from alternative implementations:

fggm <- geeglm(Weight ~ Cu * (Time + I(Time^2)), id = Pig, data = dietox,
			   family = gaussian, corstr = "exchangeable")

fgee <- gee(Weight ~ Cu * (Time + I(Time^2)), id = Pig, data = dietox,
			   family = gaussian, corstr = "exchangeable")

fygs <- yags(Weight ~ Cu * (Time + I(Time^2)), id = Pig, data = dietox,
			   family = gaussian, corstr = "exchangeable",
			   alphainit = 0.01)

model.sel(fggm, fgee, fygs, rank = QIC)

QIC(fggm, fgee, fygs, typeR = TRUE)
QIC(fggm, fgee, fygs, typeR = FALSE)

system.time(dd.ggm <- dredge(fggm, rank = QIC, ct.args = list(type = "robust"),
	   fixed = ~Cu))
system.time(dd.gee <- dredge(fgee, rank = QIC, ct.args = list(type = "robust"),
	   fixed = ~Cu))
system.time(dd.ygs <- dredge(fygs, rank = QIC, ct.args = list(type = "robust"),
	   fixed = ~Cu))

# 'geeglm' seems to be the slowest, and the fitted models stand out slightly
# from the other two.

dd.ggm
dd.gee
dd.ygs


(dd.gee.n <- dredge(fgee, rank = QIC, ct.args = list(type = "naive"),
	   fixed = ~Cu, typeR = T))
(dd.gee.n <- dredge(fgee, rank = QIC, ct.args = list(type = "naive"),
	   fixed = ~Cu, typeR = F))

# model averaged parameters (with naive covariance)
# note use of ct.args argument
model.avg(dd.gee.n)
# model averaged parameters (with robust covariance)
model.avg(dd.gee)

# the same result, but re-fitting the models
models <- get.models(dd.gee, subset = NA)
summary(mavg <- model.avg(models, rank = QIC, ct.args = list(type = "naive")))
summary(mavg <- model.avg(models, rank = QIC, ct.args = list(type = "robust")))

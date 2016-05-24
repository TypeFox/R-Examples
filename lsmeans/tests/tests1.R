# Tests of lsmeans for lm and mlm objects

require(lsmeans)

# ---------- multivariate ---------------------------------

MOats.lm <- lm (yield ~ Block + Variety, data = MOats)
MOats.rg <- ref.grid (MOats.lm, 
                mult.levs = list(nitro = c(0,.2,.4,.6)))
lsmeans(MOats.rg, ~ nitro | Variety)

# Try putting missing values whenever Yield is "Marvellous"
# plus another one for good measure
mo = MOats
mo$yield[mo$Variety == "Marvellous", 3] <- NA
mo$yield[2,4] <- NA
mo.lm <- lm (yield ~ Block + Variety, data = mo)
lsmeans(mo.lm, "Variety")

# Same as above, but use na.exclude
## In R 3.0.2, this will return NAs for the SEs and test stats
## Reported as Bug 15693 - should be fixed in later versions
mo.excl.lm <- lm (yield ~ Block + Variety, data = mo, na.action = na.exclude)
lsmeans(mo.excl.lm, "Variety")


# ------------ univariate -------------
# make an unbalanced, collinear, dataset with covariates
set.seed(19841776)
warp <- warpbreaks[-c(1,2,3,5,8,13,21,34), ]
warp$x1 <- rnorm(nrow(warp), 17, 3)
warp$x2 <- warp$x1^3 / 1007
warp.lm <- lm(breaks ~ poly(x1,3) + x2 + wool*tension, data=warp)
# Note: This model is not full-rank
( warp.lsm <- lsmeans(warp.lm, "tension", by = "wool") )
# (Nothing is estimable)

# However, contrasts ARE estimable:
(warp.pairs <- pairs(warp.lsm))

#switcheroo of by variables:
(tmp = pairs(warp.lsm, by = "tension"))

# compare these contrasts
pairs(tmp, by = "contrast")

# Joint tests
test(warp.pairs, joint = TRUE)  # all 6 but reduces to 4

test(warp.pairs, joint = TRUE, rows = 1:3)  # just wool A

test(warp.pairs, joint = TRUE, rows = 2:3)  # just wool A but not lin dep
                                            # should be identical result

test(warp.pairs, joint = TRUE, rows = 4:5)  # just wool B


# Test different ways of accessing data
## ... using "with" ...
warp.lm2 <- with(warp, lm(breaks ~ x1 + x2 + wool*tension))
lsmeans(warp.lm2, ~ tension)

## ... using "attach" ...
attach(warp)
warp.lm3 <- lm(breaks ~ x1 + x2 + wool*tension)
lsmeans(warp.lm3, "tension")

detach("warp")
# won't work if detached
try(lsmeans(warp.lm3, "tension")) 

# However, we're OK again if we use 'data'
lsmeans(warp.lm3, "tension", data = warp)


# --- aovlist objects ----
# dataset borrowed from pbkrtest
beets <- data.frame (
    harvest = factor(rep(rep(c("h1","h2"), each=5) , 3)),
    block = factor(rep(rep(c("blk1","blk2","blk3"), each=5), 2)),
    sow = factor(letters[c(3,4,5,2,1,3,2,4,5,1,5,2,3,4,1,
                           2,1,5,4,3,4,1,3,2,5,1,4,3,2,5)]),
    yield = c(128,118,95,131,136.5,136.5,150,140,99.5,156,
              99,128,126,120.5,137.5,147,153.5,100,139,141,
              115.5,135,130,134,91.5,155,140.5,142,151,102) )
# Use true contrasts for coding...
old.opt <- options(contrasts = c("contr.helmert","contr.poly"))
beets.aov <- aov(yield ~ harvest*sow + Error(block/harvest), data=beets)
lsmeans(beets.aov, eff ~ sow | harvest)

# restore old 'contrasts' that confound the intercept
options(old.opt)


# --------------- Other stuff -------------------
# using cld
cld(lsmeans(warp.lm2, ~ tension | wool))

# passing to glht
require(multcomp)
# This will fail because glht can't deal with rank deficiency
# Hope this changes.
try( as.glht(pairs(warp.lsm)) )

# However, warp.lm2 isn't rank-deficient
warp.lsm2 <- lsmeans(warp.lm2, ~ tension)
warp.con <- contrast(warp.lsm2, "eff")
summary(warp.con, adjust = "sidak")
summary(as.glht(warp.con))

summary(glht(warp.lm2, lsm(eff ~ tension | wool)))

# confint
confint(contrast(warp.lsm2, "trt.vs.ctrl1"))

# lstrends
warp.lm4 <- lm(breaks ~ tension*wool*x1, data = warp)
lstrends(warp.lm4, ~tension|wool, var = "x1")

# exotic chain rule example
lstrends(warp.lm4, ~tension|wool, var = "sqrt(x1-7)")



# -------- Transformations -------------
## ... of response ...
warp.lm5 <- lm(log(breaks) ~ x1 + x2 + tension*wool, data = warp)
warp.lsm5 <- lsmeans(warp.lm5, ~tension | wool)
summary(warp.lsm5)
summary(warp.lsm5, type = "resp")

## In a GLM
# One of the glm examples...
d.AD <- data.frame(treatment = gl(3,3), outcome = gl(3,1,9), 
    counts = c(18,17,15,20,10,20,25,13,12))
glm.D93 <- glm(counts ~ outcome + treatment, family = poisson(), data = d.AD)

( lsm.D93 <- lsmeans(glm.D93, ~ outcome) )
# un-log the results to obtain rates
summary(lsm.D93, type = "resp")

# un-log some comparisons to obtain ratios
summary(contrast(lsm.D93, "trt.vs.ctrl", ref = 2), 
	type = "resp", adjust = "none")


# weighting
nutr.lm <- lm(gain ~ (age + group + race)^2, data = nutrition)
lsmeans(nutr.lm, "race", weights = "equal")

lsmeans(nutr.lm, "race", weights = "prop")

lsmeans(nutr.lm, "race", weights = "outer")

lsmeans(nutr.lm, "race", weights = "cells")


# covariate predictions
feedlot.add <- lm(swt ~ ewt + herd + diet, data = feedlot)
lsmeans(feedlot.add, "herd", cov.reduce = ewt ~ herd)


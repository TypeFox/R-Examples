library(MonoPoly)

par <- rep(1,8)
lapply(MonoPoly:::evalGradCoef(par, ptype="El", ctype="c2"), zapsmall)
lapply(MonoPoly:::evalGradCoef(par, ptype="El", ctype="cge0"), zapsmall)
lapply(MonoPoly:::evalGradCoef(par, ptype="EHH", ctype="c2"), zapsmall)
lapply(MonoPoly:::evalGradCoef(par, ptype="EHH", ctype="cge0"), zapsmall)
lapply(MonoPoly:::evalGradCoef(par, ptype="P", ctype="c2"), zapsmall)
lapply(MonoPoly:::evalGradCoef(par, ptype="P", ctype="cge0"), zapsmall)

par <- 1:8
lapply(MonoPoly:::evalGradCoef(par, ptype="El", ctype="c2"), zapsmall)
lapply(MonoPoly:::evalGradCoef(par, ptype="El", ctype="cge0"), zapsmall)
lapply(MonoPoly:::evalGradCoef(par, ptype="EHH", ctype="c2"), zapsmall)
lapply(MonoPoly:::evalGradCoef(par, ptype="EHH", ctype="cge0"), zapsmall)
lapply(MonoPoly:::evalGradCoef(par, ptype="P", ctype="c2"), zapsmall)
lapply(MonoPoly:::evalGradCoef(par, ptype="P", ctype="cge0"), zapsmall)

par <- 8:1
lapply(MonoPoly:::evalGradCoef(par, ptype="El", ctype="c2"), zapsmall)
lapply(MonoPoly:::evalGradCoef(par, ptype="El", ctype="cge0"), zapsmall)
lapply(MonoPoly:::evalGradCoef(par, ptype="EHH", ctype="c2"), zapsmall)
lapply(MonoPoly:::evalGradCoef(par, ptype="EHH", ctype="cge0"), zapsmall)
lapply(MonoPoly:::evalGradCoef(par, ptype="P", ctype="c2"), zapsmall)
lapply(MonoPoly:::evalGradCoef(par, ptype="P", ctype="cge0"), zapsmall)

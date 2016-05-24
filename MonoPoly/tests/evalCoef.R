library(MonoPoly)

par <- rep(1,8)
MonoPoly:::evalCoef(par, ptype="El", ctype="c2")
MonoPoly:::evalCoef(par, ptype="El", ctype="cge0")
MonoPoly:::evalCoef(par, ptype="EHH", ctype="c2")
MonoPoly:::evalCoef(par, ptype="EHH", ctype="cge0")
MonoPoly:::evalCoef(par, ptype="P", ctype="c2")
MonoPoly:::evalCoef(par, ptype="P", ctype="cge0")

par <- 1:8
MonoPoly:::evalCoef(par, ptype="El", ctype="c2")
MonoPoly:::evalCoef(par, ptype="El", ctype="cge0")
MonoPoly:::evalCoef(par, ptype="EHH", ctype="c2")
MonoPoly:::evalCoef(par, ptype="EHH", ctype="cge0")
MonoPoly:::evalCoef(par, ptype="P", ctype="c2")
MonoPoly:::evalCoef(par, ptype="P", ctype="cge0")

par <- 8:1
MonoPoly:::evalCoef(par, ptype="El", ctype="c2")
MonoPoly:::evalCoef(par, ptype="El", ctype="cge0")
MonoPoly:::evalCoef(par, ptype="EHH", ctype="c2")
MonoPoly:::evalCoef(par, ptype="EHH", ctype="cge0")
MonoPoly:::evalCoef(par, ptype="P", ctype="c2")
MonoPoly:::evalCoef(par, ptype="P", ctype="cge0")

### R code from vignette source 'residual-shadings.Rnw'

###################################################
### code chunk number 1: preliminaries
###################################################
library("grid")
library("vcd")
rseed <- 1071


###################################################
### code chunk number 2: Arthritis-data
###################################################
data("Arthritis", package = "vcd")
(art <- xtabs(~ Treatment + Improved, data = Arthritis, subset = Sex == "Female"))


###################################################
### code chunk number 3: Arthritis-classic (eval = FALSE)
###################################################
## mosaic(art)
## assoc(art)


###################################################
### code chunk number 4: Arthritis-classic1
###################################################
grid.newpage()
pushViewport(viewport(layout = grid.layout(1, 2)))
pushViewport(viewport(layout.pos.col=1, layout.pos.row=1))
mosaic(art, newpage = FALSE, margins = c(2.5, 4, 2.5, 3))
popViewport()
pushViewport(viewport(layout.pos.col=2, layout.pos.row=1))
assoc(art, newpage = FALSE, margins = c(5, 2, 5, 4))
popViewport(2)


###################################################
### code chunk number 5: Arthritis-max
###################################################
set.seed(rseed)
(art_max <- coindep_test(art, n = 5000))


###################################################
### code chunk number 6: Arthritis-sumsq
###################################################
ss <- function(x) sum(x^2)
set.seed(rseed)
coindep_test(art, n = 5000, indepfun = ss)


###################################################
### code chunk number 7: Arthritis-extended (eval = FALSE)
###################################################
## mosaic(art, gp = shading_Friendly(lty = 1, eps = NULL))
## mosaic(art, gp = shading_hsv, gp_args = list(
##   interpolate = art_max$qdist(c(0.9, 0.99)), p.value = art_max$p.value))
## set.seed(rseed)
## mosaic(art, gp = shading_max, gp_args = list(n = 5000))


###################################################
### code chunk number 8: pistonrings-data
###################################################
data("pistonrings", package = "HSAUR")
pistonrings


###################################################
### code chunk number 9: shadings
###################################################
mymar <- c(1.5, 0.5, 0.5, 2.5)
grid.newpage()
pushViewport(viewport(layout = grid.layout(2, 3)))
pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 1))
mosaic(art, margins = mymar, newpage = FALSE,
  gp = shading_Friendly(lty = 1, eps = NULL))
popViewport()
pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 2))
mosaic(art, gp = shading_hsv, margins = mymar, newpage = FALSE,
  gp_args = list(interpolate = art_max$qdist(c(0.9, 0.99)), p.value = art_max$p.value))
popViewport()
pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 3))
set.seed(rseed)
mosaic(art, gp = shading_max, margins = mymar, newpage = FALSE,
  gp_args = list(n = 5000))
popViewport()

pushViewport(viewport(layout.pos.row = 2, layout.pos.col = 1))
mosaic(pistonrings, margins = mymar, newpage = FALSE,
  gp = shading_Friendly(lty = 1, eps = NULL, interpolate = c(1, 1.5)))
popViewport()
pushViewport(viewport(layout.pos.row = 2, layout.pos.col = 2))
mosaic(pistonrings, gp = shading_hsv, margins = mymar, newpage = FALSE,
  gp_args = list(p.value = 0.069, interpolate = c(1, 1.5)))
popViewport()
pushViewport(viewport(layout.pos.row = 2, layout.pos.col = 3))
mosaic(pistonrings, gp = shading_hcl, margins = mymar, newpage = FALSE,
  gp_args = list(p.value = 0.069, interpolate = c(1, 1.5)))
popViewport(2)


###################################################
### code chunk number 10: pistonrings-inference
###################################################
set.seed(rseed)
coindep_test(pistonrings, n = 5000)
set.seed(rseed)
(pring_ss <- coindep_test(pistonrings, n = 5000, indepfun = ss))


###################################################
### code chunk number 11: pistonrings-plot (eval = FALSE)
###################################################
## mosaic(pistonrings, gp = shading_Friendly(lty = 1, eps = NULL, interpolate = c(1, 1.5)))
## mosaic(pistonrings, gp = shading_hsv, gp_args = list(p.value = pring_ss$p.value, interpolate = c(1, 1.5)))
## mosaic(pistonrings, gp = shading_hcl, gp_args = list(p.value = pring_ss$p.value, interpolate = c(1, 1.5)))


###################################################
### code chunk number 12: alzheimer-data
###################################################
data("alzheimer", package = "coin")
alz <- xtabs(~ smoking + disease + gender, data = alzheimer)
alz


###################################################
### code chunk number 13: alzheimer-plot1
###################################################
set.seed(rseed)
cotabplot(~ smoking + disease | gender, data = alz, panel = cotab_coindep, n = 5000)


###################################################
### code chunk number 14: alzheimer-inference
###################################################
set.seed(rseed)
coindep_test(alz, 3, n = 5000)
set.seed(rseed)
coindep_test(alz, 3, n = 5000, indepfun = ss)
set.seed(rseed)
coindep_test(alz, 3, n = 5000, indepfun = ss, aggfun = sum)


###################################################
### code chunk number 15: alzheimer-plot (eval = FALSE)
###################################################
## set.seed(rseed)
## cotabplot(~ smoking + disease | gender, data = alz, panel = cotab_coindep, n = 5000)


###################################################
### code chunk number 16: Punishment-data
###################################################
data("Punishment", package = "vcd")
pun <- xtabs(Freq ~ memory + attitude + age + education, data = Punishment)
ftable(pun, row.vars = c("age", "education", "memory"))


###################################################
### code chunk number 17: Punishment-assoc1
###################################################
set.seed(rseed)
cotabplot(~ memory + attitude | age + education, data = pun, panel = cotab_coindep,
  n = 5000, type = "assoc", test = "maxchisq", interpolate = 1:2)


###################################################
### code chunk number 18: Punishment-mosaic1
###################################################
set.seed(rseed)
cotabplot(~ memory + attitude | age + education, data = pun, panel = cotab_coindep,
  n = 5000, type = "mosaic", test = "maxchisq", interpolate = 1:2)


###################################################
### code chunk number 19: Punishment-inference
###################################################
set.seed(rseed)
coindep_test(pun, 3:4, n = 5000)
set.seed(rseed)
coindep_test(pun, 3:4, n = 5000, indepfun = ss)
set.seed(rseed)
coindep_test(pun, 3:4, n = 5000, indepfun = ss, aggfun = sum)


###################################################
### code chunk number 20: Punishment-assoc (eval = FALSE)
###################################################
## set.seed(rseed)
## cotabplot(~ memory + attitude | age + education, data = pun, panel = cotab_coindep,
##   n = 5000, type = "assoc", test = "maxchisq", interpolate = 1:2)


###################################################
### code chunk number 21: Punishment-mosaic (eval = FALSE)
###################################################
## set.seed(rseed)
## cotabplot(~ memory + attitude | age + education, data = pun, panel = cotab_coindep,
##   n = 5000, type = "mosaic", test = "maxchisq", interpolate = 1:2)



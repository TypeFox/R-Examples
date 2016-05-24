### R code from vignette source 'strucplot.Rnw'

###################################################
### code chunk number 1: preliminaries
###################################################
set.seed(1071)
library(grid)
library(vcd)
data(Titanic)
data(HairEyeColor)
data(PreSex)
data(Arthritis)
art <- xtabs(~Treatment + Improved, data = Arthritis)


###################################################
### code chunk number 2: Arthritis
###################################################
mosaic(art, gp = shading_max, split_vertical = TRUE)


###################################################
### code chunk number 3: UCBAdmissions
###################################################
cotabplot(UCBAdmissions, panel = cotab_coindep, shade = TRUE, legend = FALSE,
          type = "assoc")


###################################################
### code chunk number 4: PreSex
###################################################
presextest <- coindep_test(PreSex, margin = c(1,4),
                           indepfun = function(x) sum(x^2), n = 5000)
mosaic(PreSex, condvars = c(1, 4), shade = TRUE,
       gp_args = list(p.value = presextest$p.value))


###################################################
### code chunk number 5: Titanic
###################################################
doubledecker(Survived ~ ., data = Titanic, labeling_args =
             list(set_varnames = c(Sex = "Gender")))


###################################################
### code chunk number 6: vcdlayout
###################################################
    pushViewport(vcd:::vcdViewport(legend = T, mar =4))
    seekViewport("main")
    grid.rect(gp = gpar(lwd = 3))
    grid.text("main", gp = gpar(fontsize = 20))
    seekViewport("sub")
    grid.rect(gp = gpar(lwd = 3))
    grid.text("sub", gp = gpar(fontsize = 20))
    seekViewport("plot")
    grid.rect(gp = gpar(lwd = 3))
    grid.text("plot", gp = gpar(fontsize = 20))
    seekViewport("legend")
    grid.text("legend", rot = 90, gp = gpar(fontsize = 20))
    grid.rect(gp = gpar(lwd = 3))
    seekViewport("legend_sub")
    grid.rect(gp = gpar(lwd = 3))
    grid.text("[F]", gp = gpar(fontsize = 20))
    seekViewport("legend_top")
    grid.rect(gp = gpar(lwd = 3))
    grid.text("[E]", gp = gpar(fontsize = 20))
    seekViewport("margin_top")
    grid.rect(gp = gpar(lwd = 3))
    grid.text("margin_top", gp = gpar(fontsize = 20))
    seekViewport("margin_bottom")
    grid.rect(gp = gpar(lwd = 3))
    grid.text("margin_bottom", gp = gpar(fontsize = 20))
    seekViewport("margin_right")
    grid.rect(gp = gpar(lwd = 3))
    grid.text("margin_right", rot = 90, gp = gpar(fontsize = 20))
    seekViewport("margin_left")
    grid.rect(gp = gpar(lwd = 3))
    grid.text("margin_left", rot = 90, gp = gpar(fontsize = 20))
    seekViewport("corner_top_left")
    grid.rect(gp = gpar(lwd = 3))
    grid.text("[A]", gp = gpar(fontsize = 20))
    seekViewport("corner_top_right")
    grid.rect(gp = gpar(lwd = 3))
    grid.text("[B]", gp = gpar(fontsize = 20))
    seekViewport("corner_bottom_left")
    grid.rect(gp = gpar(lwd = 3))
    grid.text("[C]", gp = gpar(fontsize = 20))
    seekViewport("corner_bottom_right")
    grid.rect(gp = gpar(lwd = 3))
    grid.text("[D]", gp = gpar(fontsize = 20))


###################################################
### code chunk number 7: structable
###################################################
(HEC <- structable(Eye ~ Sex + Hair, data = HairEyeColor))


###################################################
### code chunk number 8: Observed
###################################################
mosaic(HEC)


###################################################
### code chunk number 9: Observed2
###################################################
mosaic(~ Sex + Eye + Hair, data = HairEyeColor)


###################################################
### code chunk number 10: Observedfig
###################################################
mosaic(HEC)


###################################################
### code chunk number 11: Expected
###################################################
mosaic(HEC, type = "expected")


###################################################
### code chunk number 12: Expectedfig
###################################################
mosaic(HEC, type = "expected")


###################################################
### code chunk number 13: sieve
###################################################
sieve(~ Sex + Eye + Hair, data = HEC, spacing = spacing_dimequal(c(2,0,0)))


###################################################
### code chunk number 14: sievefig
###################################################
sieve(~ Sex + Eye + Hair, data = HEC, spacing = spacing_dimequal(c(2,0,0)))


###################################################
### code chunk number 15: Residuals
###################################################
assoc(HEC, compress = FALSE)


###################################################
### code chunk number 16: Residualsfig
###################################################
assoc(HEC, compress = FALSE)


###################################################
### code chunk number 17: strucplot.Rnw:592-593
###################################################
options(width=60)


###################################################
### code chunk number 18: split1
###################################################
mosaic(HEC, split_vertical = c(TRUE, FALSE, TRUE),
labeling_args = list(abbreviate_labs = c(Eye = 3)))


###################################################
### code chunk number 19: strucplot.Rnw:601-602
###################################################
options(width=70)


###################################################
### code chunk number 20: splitfig
###################################################
mosaic(HEC, split_vertical = c(TRUE, FALSE, TRUE),
labeling_args = list(abbreviate_labs = c(Eye = 3)))


###################################################
### code chunk number 21: split2
###################################################
mosaic(HEC, direction = c("v","h","v"))


###################################################
### code chunk number 22: doubledecker1
###################################################
doubledecker(Titanic)


###################################################
### code chunk number 23: doubledecker2
###################################################
doubledecker(Survived ~ Class + Sex + Age, data = Titanic)


###################################################
### code chunk number 24: strucplot.Rnw:665-666
###################################################
options(width=75)


###################################################
### code chunk number 25: subsetting
###################################################
(STD <- structable(~ Sex + Class + Age, data = Titanic[,,2:1,]))
STD["Male",]
STD["Male", c("1st","2nd","3rd")]


###################################################
### code chunk number 26: strucplot.Rnw:675-676
###################################################
options(width=70)


###################################################
### code chunk number 27: conditioning
###################################################
STD[["Male",]]
STD[[c("Male", "Adult"),]]
STD[["Male","1st"]]


###################################################
### code chunk number 28: Variables1
###################################################
pushViewport(viewport(layout = grid.layout(ncol = 2)))


###################################################
### code chunk number 29: Variables2
###################################################
pushViewport(viewport(layout.pos.col = 1))
mosaic(STD[["Male"]], margins = c(left = 2.5, top = 2.5, 0), sub = "Male", newpage = FALSE)
popViewport()


###################################################
### code chunk number 30: Variables3
###################################################
pushViewport(viewport(layout.pos.col = 2))
mosaic(STD[["Female"]], margins = c(top = 2.5, 0), sub = "Female", newpage = FALSE)
popViewport(2)


###################################################
### code chunk number 31: Variablesfig
###################################################
pushViewport(viewport(layout = grid.layout(ncol = 2)))
pushViewport(viewport(layout.pos.col = 1))
mosaic(STD[["Male"]], margins = c(left = 2.5, top = 2.5, 0), sub = "Male", newpage = FALSE)
popViewport()
pushViewport(viewport(layout.pos.col = 2))
mosaic(STD[["Female"]], margins = c(top = 2.5, 0), sub = "Female", newpage = FALSE)
popViewport(2)


###################################################
### code chunk number 32: cotabplot
###################################################
cotabplot(~ Class + Age | Sex, data = STD, split_vertical = TRUE)


###################################################
### code chunk number 33: cotabplotfig
###################################################
cotabplot(~ Class + Age | Sex, data = STD, split_vertical = TRUE)


###################################################
### code chunk number 34: Conditioning1
###################################################
mosaic(STD, condvars = "Sex", split_vertical = c(TRUE, TRUE, FALSE))


###################################################
### code chunk number 35: Conditioning2
###################################################
mosaic(~ Class + Age | Sex, data = STD, split_vertical = c(TRUE, TRUE, FALSE))


###################################################
### code chunk number 36: Conditioningfig
###################################################
mosaic(~ Class + Age | Sex, data = STD, split_vertical = c(TRUE, TRUE, FALSE))


###################################################
### code chunk number 37: pairs
###################################################
pairs(STD, highlighting = 2, diag_panel = pairs_diagonal_mosaic,
      diag_panel_args = list(fill = grey.colors))


###################################################
### code chunk number 38: pairsfig
###################################################
pairs(STD, highlighting = 2, diag_panel = pairs_diagonal_mosaic,
      diag_panel_args = list(fill = grey.colors))


###################################################
### code chunk number 39: viewportnames
###################################################
mosaic(~ Hair + Eye, data = HEC, pop = FALSE)

seekViewport("cell:Hair=Blond")
grid.rect(gp = gpar(col = "red", lwd = 4))

seekViewport("cell:Hair=Blond,Eye=Blue")
grid.circle(r = 0.2, gp = gpar(fill = "cyan"))


###################################################
### code chunk number 40: viewportnamesfig
###################################################
mosaic(~ Hair + Eye, data = HEC, pop = FALSE)

seekViewport("cell:Hair=Blond")
grid.rect(gp = gpar(col = "red", lwd = 4))

seekViewport("cell:Hair=Blond,Eye=Blue")
grid.circle(r = 0.2, gp = gpar(fill = "cyan"))


###################################################
### code chunk number 41: changeplot
###################################################
assoc(Eye ~ Hair, data = HEC, pop = FALSE)
getNames()[1:6]
grid.edit("rect:Hair=Blond,Eye=Blue", gp = gpar(fill = "red"))


###################################################
### code chunk number 42: changeplotfig
###################################################
x <- tab <- margin.table(HairEyeColor, 1:2)
x[] <- "light gray"
x["Blond","Blue"] <- "Red"
assoc(tab, gp = gpar(fill = x))


###################################################
### code chunk number 43: ucb
###################################################
(ucb <- margin.table(UCBAdmissions, 1:2))

(fill_colors <- matrix(c("dark cyan","gray","gray","dark magenta"), ncol = 2))

mosaic(ucb, gp = gpar(fill = fill_colors, col = 0))


###################################################
### code chunk number 44: ucbfig
###################################################
(ucb <- margin.table(UCBAdmissions, 1:2))

(fill_colors <- matrix(c("dark cyan","gray","gray","dark magenta"), ncol = 2))

mosaic(ucb, gp = gpar(fill = fill_colors, col = 0))


###################################################
### code chunk number 45: recycling
###################################################
mosaic(Titanic, gp = gpar(fill = c("gray","dark magenta")),
                spacing = spacing_highlighting,
                labeling_args = list(abbreviate_labs = c(Age = 3), rep = c(Survived = FALSE))
      )


###################################################
### code chunk number 46: recyclingfig
###################################################
mosaic(Titanic, gp = gpar(fill = c("gray","dark magenta")),
                spacing = spacing_highlighting,
                labeling_args = list(abbreviate_labs = c(Age = 3), rep = c(Survived = FALSE))
      )


###################################################
### code chunk number 47: shading1
###################################################
expected <- independence_table(ucb)
(x <- (ucb - expected)  / sqrt(expected))

(shading1_obj <- ifelse(x > 0, "royalblue4", "mediumorchid4"))

mosaic(ucb, gp = gpar(fill = shading1_obj))


###################################################
### code chunk number 48: shading1fig
###################################################
expected <- independence_table(ucb)
(x <- (ucb - expected)  / sqrt(expected))

(shading1_obj <- ifelse(x > 0, "royalblue4", "mediumorchid4"))

mosaic(ucb, gp = gpar(fill = shading1_obj))


###################################################
### code chunk number 49: shading2
###################################################
shading2_fun <- function(x) gpar(fill = ifelse(x > 0, "royalblue4", "mediumorchid4"))


###################################################
### code chunk number 50: shading3
###################################################
mosaic(ucb, gp = shading2_fun)


###################################################
### code chunk number 51: shading3
###################################################
shading3a_fun <- function(col = c("royalblue4", "mediumorchid4")) {
  col <- rep(col, length.out = 2)
  function(x) gpar(fill = ifelse(x > 0, col[1], col[2]))
}


###################################################
### code chunk number 52: shading4
###################################################
mosaic(ucb, gp = shading3a_fun(c("royalblue4","mediumorchid4")))


###################################################
### code chunk number 53: shading4
###################################################
shading3b_fun <- function(observed = NULL, residuals = NULL, expected = NULL,
    df = NULL, col = c("royalblue4", "mediumorchid4")) {
  col <- rep(col, length.out = 2)
  function(x) gpar(fill = ifelse(x > 0, col[1], col[2]))
}
class(shading3b_fun) <- "grapcon_generator"


###################################################
### code chunk number 54: shading5
###################################################
mosaic(ucb, gp = shading3b_fun, gp_args = list(col = c("red","blue")))


###################################################
### code chunk number 55: haireye1
###################################################
haireye <- margin.table(HairEyeColor, 1:2)
mosaic(haireye, gp = shading_hsv)


###################################################
### code chunk number 56: haireye2
###################################################
mosaic(haireye, gp = shading_hcl)


###################################################
### code chunk number 57: haireye3
###################################################
mosaic(haireye, gp = shading_hcl, gp_args = list(h = c(130, 43), c = 100, l = c(90, 70)))


###################################################
### code chunk number 58: haireyefig1
###################################################
mosaic(haireye, gp = shading_hsv, margin = c(bottom = 1), keep_aspect_ratio = FALSE)


###################################################
### code chunk number 59: haireyefig2
###################################################
mosaic(haireye, gp = shading_hcl, margin = c(bottom = 1), keep_aspect_ratio = FALSE)


###################################################
### code chunk number 60: haireyefig3
###################################################
mosaic(haireye, gp = shading_hcl, margin = c(bottom = 1), gp_args = list(h = c(130, 43), c = 100, l = c(90, 70)), keep_aspect_ratio = FALSE)


###################################################
### code chunk number 61: interpolate
###################################################
mosaic(haireye, shade = TRUE, gp_args = list(interpolate = 1:4))


###################################################
### code chunk number 62: continuous1
###################################################
ipol <- function(x) pmin(x/4, 1)


###################################################
### code chunk number 63: continuous2
###################################################
mosaic(haireye, shade = TRUE, gp_args = list(interpolate = ipol),
labeling_args = list(abbreviate_labs = c(Sex = TRUE)))


###################################################
### code chunk number 64: interpolatefig
###################################################
pushViewport(viewport(layout = grid.layout(ncol = 2)))
pushViewport(viewport(layout.pos.col = 1))
mosaic(haireye, gp_args = list(interpolate = 1:4), margin = c(right = 1), keep_aspect_ratio= FALSE,newpage = FALSE,legend_width=5.5,shade = TRUE)
popViewport(1)
pushViewport(viewport(layout.pos.col = 2))
mosaic(haireye, gp_args = list(interpolate = ipol), margin = c(left=3,right = 1), keep_aspect_ratio = FALSE, newpage = FALSE, shade = TRUE)
popViewport(2)


###################################################
### code chunk number 65: bundesliga
###################################################
BL <- xtabs(~ HomeGoals + AwayGoals, data = Bundesliga, subset = Year == 1995)
mosaic(BL, shade = TRUE)


###################################################
### code chunk number 66: friendly
###################################################
mosaic(BL, gp = shading_Friendly, legend = legend_fixed, zero_size = 0)


###################################################
### code chunk number 67: bundesligafig
###################################################
pushViewport(viewport(layout = grid.layout(ncol = 2)))
pushViewport(viewport(layout.pos.col = 1))
mosaic(BL, margin = c(right = 1), keep_aspect_ratio= FALSE, newpage = FALSE,
       legend_width=5.5, shade = TRUE)
popViewport(1)
pushViewport(viewport(layout.pos.col = 2))
mosaic(BL, gp = shading_Friendly, legend = legend_fixed, zero_size = 0,
       margin = c(right = 1), keep_aspect_ratio= FALSE, newpage = FALSE,
       legend_width=5.5)
popViewport(2)


###################################################
### code chunk number 68: arthritis
###################################################
set.seed(4711)
mosaic(~ Treatment + Improved, data = Arthritis, subset = Sex == "Female",
       gp = shading_max)


###################################################
### code chunk number 69: arthritisfig
###################################################
set.seed(4711)
mosaic(~ Treatment + Improved, data = Arthritis, subset = Sex == "Female",
       gp = shading_max)


###################################################
### code chunk number 70: default
###################################################
mosaic(Titanic)


###################################################
### code chunk number 71: clipping
###################################################
mosaic(Titanic, labeling_args = list(clip = c(Survived = TRUE, Age = TRUE)))


###################################################
### code chunk number 72: abbreviating
###################################################
mosaic(Titanic, labeling_args = list(abbreviate_labs = c(Survived = TRUE, Age = 3)))


###################################################
### code chunk number 73: rotate
###################################################
mosaic(Titanic, labeling_args = list(rot_labels = c(bottom = 90, right = 0),
offset_varnames = c(right = 1), offset_labels = c(right = 0.3)),
margins = c(right = 4, bottom = 3))


###################################################
### code chunk number 74: repeat
###################################################
mosaic(Titanic, labeling_args = list(rep = c(Survived = FALSE, Age = FALSE)))


###################################################
### code chunk number 75: label1fig
###################################################
pushViewport(viewport(layout = grid.layout(ncol = 2,nrow=3)))
pushViewport(viewport(layout.pos.col = 1, layout.pos.row = 1))
mosaic(Titanic, newpage = FALSE, keep = TRUE, margin = c(right = 3),
       gp_labels = gpar(fontsize = 10))
popViewport()
pushViewport(viewport(layout.pos.col = 2, layout.pos.row = 1))
mosaic(Titanic, labeling_args = list(clip = c(Survived = TRUE, Age = TRUE)),
       newpage = FALSE, keep = TRUE, margin = c(left = 3), gp_labels = gpar(fontsize = 10))
popViewport()
pushViewport(viewport(layout.pos.col = 1, layout.pos.row = 2))
mosaic(Titanic, labeling_args = list(abbreviate_labs = c(Survived = TRUE, Age = 2)),
       newpage = FALSE, keep = TRUE, margin = c(right = 3), gp_labels = gpar(fontsize = 10))
popViewport()
pushViewport(viewport(layout.pos.col = 2, layout.pos.row = 2))
mosaic(Titanic, labeling_args = list(rep = c(Survived = FALSE, Age = FALSE)),
       newpage = FALSE, keep = TRUE, margin = c(left = 3), gp_labels = gpar(fontsize = 10))
popViewport()
pushViewport(viewport(layout.pos.col = 1:2, layout.pos.row = 3))
pushViewport(viewport(width = 0.55))
mosaic(Titanic, labeling_args = list(rot_labels = c(bottom = 90, right = 0),
                  offset_varnames = c(right = 1), offset_labels = c(right = 0.3)),
       margins = c(right = 4, bottom = 3), newpage = FALSE, keep = FALSE, gp_labels = gpar(fontsize = 10))
popViewport(3)


###################################################
### code chunk number 76: left
###################################################
mosaic(Titanic, labeling_args = list(pos_varnames = "left", pos_labels = "left",
       just_labels = "left", rep = FALSE))


###################################################
### code chunk number 77: left2
###################################################
mosaic(Titanic, labeling = labeling_left)


###################################################
### code chunk number 78: margins
###################################################
mosaic(Titanic, labeling_args = list(tl_labels = FALSE, tl_varnames = TRUE, abbreviate_labs = c(Survived = 1, Age = 3)))


###################################################
### code chunk number 79: boxes
###################################################
mosaic(Titanic, labeling_args = list(tl_labels = FALSE, tl_varnames = TRUE,
       boxes = TRUE, clip = TRUE))


###################################################
### code chunk number 80: boxes2
###################################################
mosaic(Titanic, labeling = labeling_cboxed)


###################################################
### code chunk number 81: labbl
###################################################
mosaic(Titanic, labeling_args = list(tl_labels = TRUE, boxes = TRUE,
       clip = c(Survived = FALSE, Age = FALSE, TRUE), abbreviate_labs = c(Age = 4),
       labbl_varnames = TRUE), margins = c(left = 4, right = 1, 3))


###################################################
### code chunk number 82: labbl2
###################################################
mosaic(Titanic, labeling = labeling_lboxed, margins = c(right = 4, left = 1, 3))


###################################################
### code chunk number 83: label2fig
###################################################
pushViewport(viewport(layout = grid.layout(ncol = 2, nrow = 2)))
pushViewport(viewport(layout.pos.col = 1, layout.pos.row = 1))
mosaic(Titanic, labeling_args = list(pos_varnames = "left", pos_labels = "left",
       just_labels = "left", rep = FALSE), newpage = FALSE, keep = TRUE,
       gp_labels = gpar(fontsize = 12))
popViewport()
pushViewport(viewport(layout.pos.col = 2, layout.pos.row = 1))
mosaic(Titanic, labeling_args = list(tl_labels = FALSE, tl_varnames = TRUE,
                  abbreviate_labs = c(Survived = 1, Age = 3)), newpage = FALSE, keep = TRUE,
       margins = c(left = 4, right = 1, 3), gp_labels = gpar(fontsize = 12))
popViewport()
pushViewport(viewport(layout.pos.col = 1, layout.pos.row = 2))
mosaic(Titanic, labeling_args = list(tl_labels = FALSE, tl_varnames = TRUE,
       boxes = TRUE, clip = TRUE), newpage = FALSE, keep = TRUE,
       gp_labels = gpar(fontsize = 12))
popViewport()
pushViewport(viewport(layout.pos.col = 2, layout.pos.row = 2))
mosaic(Titanic, labeling_args = list(tl_labels = TRUE, boxes = TRUE,
       clip = c(Survived = FALSE, Age = FALSE, TRUE),
       labbl_varnames = TRUE, abbreviate_labs = c(Age = 4)),
       margins = c(left = 4, right = 1, 3),
       newpage = FALSE, keep = TRUE, gp_labels = gpar(fontsize = 12))
popViewport(2)


###################################################
### code chunk number 84: cell
###################################################
mosaic(~ MaritalStatus + Gender, data = PreSex, labeling = labeling_cells)


###################################################
### code chunk number 85: cell2
###################################################
mosaic(~ PremaritalSex + ExtramaritalSex, data = PreSex,
labeling = labeling_cells(abbreviate_labels = TRUE, abbreviate_varnames = TRUE, clip = FALSE))


###################################################
### code chunk number 86: conditional
###################################################
mosaic(~ PremaritalSex + ExtramaritalSex | MaritalStatus + Gender, data = PreSex, labeling = labeling_conditional(abbreviate_varnames = TRUE, abbreviate_labels = TRUE, clip = FALSE, gp_text = gpar(col = "red")))


###################################################
### code chunk number 87: text
###################################################
mosaic(Titanic, labeling_args = list(abbreviate_labs = c(Survived = 1, Age = 4)), pop = FALSE)

tab <- ifelse(Titanic < 6, NA, Titanic)
labeling_cells(text = tab, clip = FALSE)(Titanic)


###################################################
### code chunk number 88: label3fig
###################################################
grid.newpage()
pushViewport(viewport(layout = grid.layout(ncol = 2, nrow = 2)))
pushViewport(viewport(layout.pos.col = 1, layout.pos.row = 1))
mosaic(~ MaritalStatus + Gender, data = PreSex, labeling = labeling_cells, newpage = FALSE, keep = TRUE, gp_labels = gpar(fontsize = 10))
popViewport()
pushViewport(viewport(layout.pos.col = 2, layout.pos.row = 1))
mosaic(~ PremaritalSex + ExtramaritalSex, data = PreSex,
       labeling = labeling_cells(abbreviate_labels = TRUE, abbreviate_varnames = TRUE,
         clip = FALSE), newpage = FALSE, keep = TRUE, gp_labels = gpar(fontsize = 10))
popViewport()
pushViewport(viewport(layout.pos.col = 1, layout.pos.row = 2))
mosaic(~ PremaritalSex + ExtramaritalSex | MaritalStatus + Gender, data = PreSex,
       labeling = labeling_conditional(abbreviate_varnames = TRUE,
         abbreviate_labels = TRUE, clip = FALSE, gp_text = gpar(col = "red")),
       newpage = FALSE, keep = TRUE, gp_labels = gpar(fontsize = 10))
popViewport()
pushViewport(viewport(layout.pos.col = 2, layout.pos.row = 2))
mosaic(Titanic, labeling_args = list(abbreviate_labs = c(Survived = 1, Age = 3)),
       pop = FALSE, newpage = FALSE, keep = TRUE, gp_labels = gpar(fontsize = 10))
tab <- ifelse(Titanic < 6, NA, Titanic)
labeling_cells(text = tab, clip = FALSE)(Titanic)


###################################################
### code chunk number 89: list
###################################################
mosaic(Titanic, labeling = labeling_list, margins = c(bottom = 5))


###################################################
### code chunk number 90: listfig
###################################################
mosaic(Titanic, labeling = labeling_list, margins = c(bottom = 5), keep = TRUE)


###################################################
### code chunk number 91: artspine
###################################################
(art <- structable(~Treatment + Improved, data = Arthritis, split_vertical = TRUE))
(my_spacing <- list(unit(0.5, "lines"), unit(c(0, 0), "lines")))
my_colors <- c("lightgray", "lightgray", "black")
mosaic(art, spacing = my_spacing, gp = gpar(fill = my_colors, col = my_colors))


###################################################
### code chunk number 92: artspinefig
###################################################
(art <- structable(~Treatment + Improved, data = Arthritis, split_vertical = TRUE))
(my_spacing <- list(unit(0.5, "lines"), unit(c(0, 0), "lines")))
my_colors <- c("lightgray", "lightgray", "black")
mosaic(art, spacing = my_spacing, gp = gpar(fill = my_colors, col = my_colors))


###################################################
### code chunk number 93: artspine
###################################################
mosaic(Improved ~ Treatment, data = Arthritis, split_vertical = TRUE)


###################################################
### code chunk number 94: space1
###################################################
mosaic(art, spacing = spacing_equal(unit(2, "lines")))


###################################################
### code chunk number 95: space2
###################################################
mosaic(art, spacing = spacing_dimequal(unit(1:2, "lines")))


###################################################
### code chunk number 96: space3
###################################################
mosaic(art, spacing = spacing_increase(start = unit(0.5, "lines"), rate = 1.5))


###################################################
### code chunk number 97: spine4
###################################################
mosaic(art, spacing = spacing_highlighting, gp = my_colors)


###################################################
### code chunk number 98: spacingfig
###################################################
pushViewport(viewport(layout = grid.layout(ncol = 2, nrow = 2)))
pushViewport(viewport(layout.pos.col = 1, layout.pos.row = 1))
mosaic(art, spacing = spacing_equal(unit(2, "lines")), keep = TRUE, newpage = FALSE)
popViewport()
pushViewport(viewport(layout.pos.col = 2, layout.pos.row = 1))
mosaic(art, spacing = spacing_dimequal(unit(c(0.5, 2), "lines")), keep = TRUE, newpage = FALSE)
popViewport()
pushViewport(viewport(layout.pos.col = 1, layout.pos.row = 2))
mosaic(art, spacing = spacing_increase(start = unit(0.3, "lines"), rate = 2.5), keep = TRUE, newpage = FALSE)
popViewport()
pushViewport(viewport(layout.pos.col = 2, layout.pos.row = 2))
mosaic(art, spacing = spacing_highlighting, keep = TRUE, newpage = FALSE)
popViewport(2)


###################################################
### code chunk number 99: oc1
###################################################
tab <- xtabs(Freq ~ stage + operation + xray + survival, data = OvaryCancer)


###################################################
### code chunk number 100: oc2
###################################################
structable(survival ~ ., data = tab)


###################################################
### code chunk number 101: oc3
###################################################
dpa <- list(var_offset = 1.2, rot = -30, just_leveltext= "left")
pairs(tab, diag_panel = pairs_barplot, diag_panel_args = dpa)


###################################################
### code chunk number 102: ocpairs
###################################################
dpa <- list(var_offset = 1.2, rot = -30, just_leveltext= "left")
pairs(tab, diag_panel = pairs_barplot, diag_panel_args = dpa)


###################################################
### code chunk number 103: oc4
###################################################
doubledecker(survival ~ stage + operation + xray, data = tab)


###################################################
### code chunk number 104: ocdoubledecker
###################################################
doubledecker(survival ~ stage + operation + xray, data = tab)


###################################################
### code chunk number 105: oc6
###################################################
split <- c(TRUE, TRUE, TRUE, FALSE)
mosaic(tab, expected = ~ survival + operation * xray * stage, split_vertical = split)


###################################################
### code chunk number 106: ocmosaicnull
###################################################
split <- c(TRUE, TRUE, TRUE, FALSE)
mosaic(tab, expected = ~ survival + operation * xray * stage, split_vertical = split)


###################################################
### code chunk number 107: oc7
###################################################
mosaic(tab, expected = ~ (survival + operation * xray) * stage, split_vertical = split)


###################################################
### code chunk number 108: ocmosaicstage
###################################################
mosaic(tab, expected = ~ (survival + operation * xray) * stage, split_vertical = split)




## @knitr echo=FALSE, message=FALSE
library(pathdiagram)


## @knitr waffles_pancakes, echo=FALSE
# manifest variables
ingredients = list(
  eggs = manifest("Eggs", x = 0.25, y = 0.8, width = 0.1, height = 0.08),
  milk = manifest("Milk", x = 0.25, y = 0.65, width = 0.1, height = 0.08),
  flour = manifest("Flour", x = 0.25, y = 0.5, width = 0.1, height = 0.08),
  sugar = manifest("Sugar", x = 0.25, y = 0.35, width = 0.1, height = 0.08),
  butter = manifest("Butter", x = 0.25, y = 0.2, width = 0.1, height = 0.08)
)

# latent variables
pancakes = latent("Pancakes", x = 0.8, y = 0.65, rx = 0.08, ry = 0.06)
waffles = latent("Waffles", x = 0.8, y = 0.35, rx = 0.08, ry = 0.06)


## @knitr waf_pan_diag, fig.cap='example of a path diagram', fig.width=6, fig.height=4, out.width='.8\\linewidth', out.height='.5\\linewidth', fig.align='center', fig.pos='h', echo=FALSE
# open a new wall
op = par(mar = rep(0.5, 4))
wall(xlim=c(0.1, 0.9), ylim=c(.15,.85))

# draw latent variables
draw(pancakes)
draw(waffles)

# draw ingredients
for (i in 1:5) {
  draw(ingredients[[i]])
}

for (i in 1:5) {
arrow(from = ingredients[[i]], to = pancakes, start = "east", end = "west")
arrow(from = ingredients[[i]], to = waffles, start = "east", end = "west")
}
par(op)


## @knitr eval=FALSE
## # installation
## install.packages("pathdiagram")


## @knitr eval=FALSE
## # load pathdiagram
## library(pathdiagram)


## @knitr manifest_variables, tidy=FALSE
# manifest variables
ingredients = list(
  eggs = manifest("Eggs", x = 0.25, y = 0.8, width = 0.1, height = 0.08),
  milk = manifest("Milk", x = 0.25, y = 0.65, width = 0.1, height = 0.08),
  flour = manifest("Flour", x = 0.25, y = 0.5, width = 0.1, height = 0.08),
  sugar = manifest("Sugar", x = 0.25, y = 0.35, width = 0.1, height = 0.08),
  butter = manifest("Butter", x = 0.25, y = 0.2, width = 0.1, height = 0.08)
)


## @knitr latent_variables
# latent variables
pancakes = latent("Pancakes", x = 0.8, y = 0.65, rx = 0.08, ry = 0.06)
waffles = latent("Waffles", x = 0.8, y = 0.35, rx = 0.08, ry = 0.06)


## @knitr diagram1, fig.width=6, fig.height=4, out.width='.8\\linewidth', out.height='.5\\linewidth', fig.align='center', echo=-c(1,11)
op = par(mar = rep(0.5, 4))
# open a new wall
wall(xlim = c(0.1, 0.9), ylim = c(0.1, 0.9))

# draw latent variables
draw(pancakes)
draw(waffles)

# draw ingredients
for (i in 1:5) {
  draw(ingredients[[i]])
}
par(op)


## @knitr add_arrows, eval=FALSE
## # arrows
## for (i in 1:5) {
## arrow(from = ingredients[[i]], to = pancakes, start = "east", end = "west")
## arrow(from = ingredients[[i]], to = waffles, start = "east", end = "west")
## }


## @knitr diagram2, echo=FALSE, fig.width=6, fig.height=4, out.width='.8\\linewidth', out.height='.5\\linewidth', fig.align='center', echo=FALSE
op = par(mar = rep(0.5, 4))
# open a new wall
wall(xlim=c(0.1, 0.9), ylim=c(.1,.9))

# draw latent variables
draw(pancakes)
draw(waffles)

# draw ingredients
for (i in 1:5) {
  draw(ingredients[[i]])
}

for (i in 1:5) {
arrow(from = ingredients[[i]], to = pancakes, start = "east", end = "west")
arrow(from = ingredients[[i]], to = waffles, start = "east", end = "west")
}
par(op)


## @knitr barcelona_blocks, tidy=FALSE
# define Attack block
attack = list(
  att1 = manifest("Messi", x=0.15, y=0.9, width=0.09, height=0.08, fill="#d199a4"),
  att2 = manifest("Xavi", x=0.15, y=0.75, width=0.09, height=0.08, fill="#d199a4"),
  att3 = manifest("Iniesta", x=0.15, y=0.6, width=0.09, height=0.08, fill="#d199a4"))
ATTACK = latent("ATTACK", x=0.35, y=0.75, rx=0.08, ry=0.07, fill="#a12b43", font=1)

# define Defense block
defense = list(
  def1 = manifest("Puyol", x=0.15, y=0.4, width=0.09, height=0.08, fill="#a0bee1"),
  def2 = manifest("Pique", x=0.15, y=0.25, width=0.09, height=0.08, fill="#a0bee1"),
  def3 = manifest("Abidal", x=0.15, y=0.1, width=0.09, height=0.08, fill="#a0bee1"))
DEFENSE = latent("DEFENSE", x=0.35, y=0.25, rx=0.08, ry=0.07, fill="#1e67ba", font=1)

# define Success block
success = list(
  suc1 = manifest("2008-2009", x=0.85, y=0.65, width=0.14, height=0.08, fill="gold2"),
  suc2 = manifest("2009-2010", x=0.85, y=0.5, width=0.14, height=0.08, fill="gold2"),
  suc3 = manifest("2010-2011", x=0.85, y=0.35, width=0.14, height=0.08, fill="gold2"))
SUCCESS = latent("SUCCESS", x=0.65, y=0.5, rx=0.08, ry=0.07, fill="gold2", font=1)


## @knitr barcelona_pathdiagram, fig.cap='FC Barcelona path diagram', fig.width=8, fig.height=4, out.width='.9\\linewidth', out.height='.55\\linewidth', fig.align='center', fig.pos='h', echo=-c(1,16)
op = par(mar = rep(0.5, 4))
# open plot window
wall(ylim=c(0.1, 0.9))

# draw latent variables
draw(ATTACK)
draw(DEFENSE)
draw(SUCCESS)

# draw manifest variables
for (i in 1:3) {
  draw(attack[[i]])
  arrow(from=attack[[i]], to=ATTACK, start="east", end="west", col="#d199a4")
  draw(defense[[i]])
  arrow(from=defense[[i]], to=DEFENSE, start="east", end="west", col="#a0bee1")
  draw(success[[i]])
  arrow(from=SUCCESS, to=success[[i]], start="east", end="west", col="gold1")
}

# arrows of inner model
arrow(from=ATTACK, to=SUCCESS, start="east", end="west", col="#d199a4")
arrow(from=DEFENSE, to=SUCCESS, start="east", end="west", col="#a0bee1")
par(op)


## @knitr Gryffindor_blocks, tidy=FALSE
# Gryffindor block
gryff = list(
  harry = manifest("Harry \nPotter", x=0.15, y=0.8, width=0.12, height=0.08, 
    cex=0.8, fill="#f2d22e", col="#7c4f87", family="serif"),
  
  ron = manifest("Ron\nWeasley", x=0.15, y=0.7, width=0.12, height=0.08, 
    cex=0.8, fill="#f2d22e", col="#7c4f87", family="serif"),
  
  hermione = manifest("Hermione\nGranger", x=0.15, y=0.6, width=0.12, height=0.08, 
    cex=0.8, fill="#f2d22e", col="#7c4f87", family="serif"),
  
  albus = manifest("Albus\nDumbledore", x=0.15, y=0.5, width=0.12, height=0.08, 
    cex=0.8, fill="#f2d22e", col="#7c4f87", family="serif"),
  
  neville = manifest("Neville\nLongbottom", x=0.15, y=0.4, width=0.12, height=0.08, 
    cex=0.8, fill="#f2d22e", col="#7c4f87", family="serif"),
  
  sirius = manifest("Sirius\nBlack", x=0.15, y=0.3, width=0.12, height=0.08, 
    cex=0.8, fill="#f2d22e", col="#7c4f87", family="serif"),
  
  rubeus = manifest("Rubeus\nHagrid", x=0.15, y=0.2, width=0.12, height=0.08, 
    cex=0.8, fill="#f2d22e", col="#7c4f87", family="serif"))


## @knitr Slytherin_blocks, tidy=FALSE
# Slytherin block
slyth = list(
  tom = manifest("Tom\nRiddle", x=0.85, y=0.8, width=0.12, height=0.08, 
    cex=0.8, fill="gray70", col="#467d70", family="serif"),
  
  severus = manifest("Severus\nSnape", x=0.85, y=0.7, width=0.12, height=0.08, 
    cex=0.8, fill="gray70", col="#467d70", family="serif"),
  
  bella = manifest("Bellatrix\nLestrange", x=0.85, y=0.6, width=0.12, height=0.08, 
    cex=0.8, fill="gray70", col="#467d70", family="serif"),
  
  regulus = manifest("Regulus\nBlack", x=0.85, y=0.5, width=0.12, height=0.08, 
    cex=0.8, fill="gray70", col="#467d70", family="serif"),
  
  phineas = manifest("Phineas\nBlack", x=0.85, y=0.4, width=0.12, height=0.08, 
    cex=0.8, fill="gray70", col="#467d70", family="serif"),
  
  draco = manifest("Draco\nMalfoy", x=0.85, y=0.3, width=0.12, height=0.08, 
    cex=0.8, fill="gray70", col="#467d70", family="serif"),
  
  horace = manifest("Horace\nSlughorn", x=0.85, y=0.2, width=0.12, height=0.08, 
    cex=0.8, fill="gray70", col="#467d70", family="serif"))


## @knitr Gry_Sly, tidy=FALSE
# latent variables
gry = latent("Gryffindor", x=0.375, y=0.5, rx=0.07, ry=0.06, cex=0.85, 
    fill="#7c4f87", family="serif")

sly = latent("Slytherin", x=0.625, y=0.5, rx=0.07, ry=0.06, cex=0.85, 
    fill="#467d70", family="serif")


## @knitr harry_potter_diagram, fig.cap='Gryffindor vs Slytherin path diagram', fig.width=6.5, fig.height=4.5, out.width='.9\\linewidth', out.height='.63\\linewidth', fig.align='center', fig.pos='h', echo=c(-1,-7)
op = par(mar = rep(0.5, 4))
# open plot window
wall(xlim=c(0.1, 0.9), ylim=c(0.15, 0.85))

# draw variables
for (i in 1:7)
{
  # arrows between each block and its latent
  arrow(from=gryff[[i]], to=gry, start="east", end="west", 
    col="#b095b7", angle=5, lwd=1)
  arrow(from=slyth[[i]], to=sly, start="west", end="east", 
    col="#90b1a9", angle=5, lwd=1)
  
  # variables
  draw(gryff[[i]])
  draw(slyth[[i]])
  draw(gry)
  draw(sly)

  # arrows between latent variables
  arrow(from=gry, to=sly, start="east", end="west", col="#dddddd", angle=20)
  arrow(from=sly, to=gry, start="west", end="east", col="#dddddd", angle=20)
}
par(op)


## @knitr my_model_blocks, tidy=FALSE
# latent variables
optimism = latent("Optimism", x=0.35, y=0.75, rx=0.08, ry=0.06, 
                  fill="gray90", col="#1B9E77", font=1)

dedication = latent("Dedication", x=0.2, y=0.6, rx=0.08, ry=0.06, 
                    fill="gray90", col="#D95F02", font=1)

patience = latent("Patience", x=0.2, y=0.4, rx=0.08, ry=0.06, 
                  fill="gray90", col="#7570B3", font=1)

sacrifice = latent("Sacrifice", x=0.35, y=0.25, rx=0.08, ry=0.06, 
                   fill="gray90", col="#E7298A", font=1)

work = latent("Work", x=0.5, y=0.5, rx=0.08, ry=0.06, 
              fill="gray90", col="#1F78B4", font=1)

achievement = latent("Achievement", x=0.8, y=0.5, rx=0.10, ry=0.075, 
                     fill="gray90", col="tomato", font=1)

luck = latent("Luck", x=0.85, y=0.7, rx=0.065, ry=0.06, 
              fill="gray90", col="#E6AB02", font=1)


## @knitr my_model_diagram, fig.cap='Personal Achievement Model path diagram', fig.width=8, fig.height=4, out.width='.9\\linewidth', out.height='.55\\linewidth', fig.align='center', fig.pos='h', echo=-c(1,21)
op = par(mar = rep(0.5, 4))
# open wall to plot
wall(ylim = c(0.15, 0.85))

# draw latent variables
draw(optimism)
draw(dedication)
draw(patience)
draw(sacrifice)
draw(work)
draw(achievement)
draw(luck)

# add arrows
arrow(from=optimism, to=work, start="east", end="north", col="gray90")
arrow(from=dedication, to=work, start="east", end="west", col="gray90")
arrow(from=patience, to=work, start="east", end="west", col="gray90")
arrow(from=sacrifice, to=work, start="east", end="south", col="gray90")
arrow(from=work, to=achievement, start="east", end="west", col="gray90")
arrow(from=luck, to=achievement, start="south", end="north", col="gray90")
par(op)


## @knitr block1, tidy=FALSE
# Block 1
X1 = list(
  x11 = manifest(expression(x[11]), x=0.4, y=0.9, width=0.06, height=0.06, 
                 cex=0.8, fill="gray90", col="gray50"),
  x12 = manifest(expression(x[12]), x=0.47, y=0.9, width=0.06, height=0.06, 
                 cex=0.8, fill="gray90", col="gray50"),
  x13 = manifest(expression(x[13]), x=0.53, y=0.9, width=0.06, height=0.06, 
                 cex=0.8, fill="gray90", col="gray50"),
  x14 = manifest(expression(x[14]), x=0.6, y=0.9, width=0.06, height=0.06, 
                 cex=0.8, fill="gray90", col="gray50")
)

LV1 = latent(expression(xi[1]), x=0.5, y=0.75, rx=0.06, ry=0.06, 
             cex=1.2, fill="gray50")


# Block 2
X2 = list(
  x21 = manifest(expression(x[21]), x=0.1, y=0.63, width=0.06, height=0.06, 
                 cex=0.8, fill="gray90", col="gray50"),
  x22 = manifest(expression(x[22]), x=0.1, y=0.55, width=0.06, height=0.06, 
                 cex=0.8, fill="gray90", col="gray50"),
  x23 = manifest(expression(x[23]), x=0.1, y=0.47, width=0.06, height=0.06, 
                 cex=0.8, fill="gray90", col="gray50")
)

LV2 = latent(expression(xi[2]), x=0.25, y=0.55, rx=0.06, ry=0.06, 
             cex=1.2, fill="gray50")


# Block 3
X3 = list(
  x31 = manifest(expression(x[31]), x=0.1, y=0.33, width=0.06, height=0.06, 
                 cex=0.8, fill="gray90", col="gray50"),
  x32 = manifest(expression(x[32]), x=0.1, y=0.25, width=0.06, height=0.06, 
                 cex=0.8, fill="gray90", col="gray50"),
  x33 = manifest(expression(x[33]), x=0.1, y=0.17, width=0.06, height=0.06, 
                 cex=0.8, fill="gray90", col="gray50")
)

LV3 = latent(expression(xi[3]), x=0.25, y=0.25, rx=0.06, ry=0.06, 
             cex=1.2, fill="gray50")


# Block 4
X4 = list(
  x41 = manifest(expression(x[41]), x=0.4, y=0.45, width=0.06, height=0.06, 
                 cex=0.8, fill="gray90", col="gray50"),
  x42 = manifest(expression(x[42]), x=0.4, y=0.35, width=0.06, height=0.06, 
                 cex=0.8, fill="gray90", col="gray50")
)

LV4 = latent(expression(xi[4]), x=0.55, y=0.4, rx=0.06, ry=0.06, 
             cex=1.2, fill="gray50")


# Block 5
X5 = list(
  x51 = manifest(expression(x[51]), x=0.95, y=0.47, width=0.06, height=0.06, 
                 cex=0.8, fill="gray90", col="gray50"),
  x52 = manifest(expression(x[52]), x=0.95, y=0.4, width=0.06, height=0.06, 
                 cex=0.8, fill="gray90", col="gray50"),
  x53 = manifest(expression(x[53]), x=0.95, y=0.33, width=0.06, height=0.06, 
                 cex=0.8, fill="gray90", col="gray50")
)

LV5 = latent(expression(xi[5]), x=0.8, y=0.4, rx=0.06, ry=0.06, 
             cex=1.2, fill="gray50")


## @knitr formal_diagram, fig.cap='Formal Notation Model path diagram', fig.width=6.5, fig.height=4, out.width='.9\\linewidth', out.height='.6\\linewidth', fig.align='center', fig.pos='h', echo=-c(1,32)
op = par(mar = rep(0.5, 4))
# open plot window
wall(xlim=c(0, 1), ylim=c(0.05, 0.95))

# block 1
draw(LV1)
for (i in 1:length(X1)) {
  draw(X1[[i]])
  arrow(from=X1[[i]], to=LV1, start="south", end="north", lwd=1, col="gray80")
}

# block 2
draw(LV2)
for (i in 1:length(X2)) {
  draw(X2[[i]])
  arrow(from=X2[[i]], to=LV2, start="east", end="west", lwd=1, col="gray80")
}

# block 3
draw(LV3)
for (i in 1:length(X3)) {
  draw(X3[[i]])
  arrow(from=X3[[i]], to=LV3, start="east", end="west", lwd=1, col="gray80")
}

# block 4
draw(LV4)
for (i in 1:length(X4)) {
  draw(X4[[i]])
  arrow(from=X4[[i]], to=LV4, start="east", end="west", lwd=1, col="gray80")
}

# block 5
draw(LV5)
for (i in 1:length(X5)) {
  draw(X5[[i]])
  arrow(from=X5[[i]], to=LV5, start="west", end="east", lwd=1, col="gray80")
}

# arrows between latent variables
arrow(from=LV1, to=LV2, start="west", end="north", col="gray80")
arrow(from=LV1, to=LV5, start="east", end="north", col="gray80")
arrow(from=LV2, to=LV3, start="south", end="north", col="gray80")
arrow(from=LV2, to=LV4, start="east", end="north", col="gray80")
arrow(from=LV3, to=LV4, start="east", end="south", col="gray80")
arrow(from=LV4, to=LV5, start="east", end="west", col="gray80")
par(op)



### R code from vignette source 'extensibility.Rnw'

###################################################
### code chunk number 1: extensibility.Rnw:32-34
###################################################
library(grid)
library(gridSVG)


###################################################
### code chunk number 2: extensibility.Rnw:60-61
###################################################
tg <- grob(name="tg", cl="timegrob")


###################################################
### code chunk number 3: extensibility.Rnw:67-71
###################################################
timegrob <- function(x) {
    textGrob(paste("text generated at", Sys.time(), sep="\n"),
             gp=x$gp, name=x$name)
}


###################################################
### code chunk number 4: extensibility.Rnw:77-80
###################################################
drawDetails.timegrob <- function(x, ...) {
    grid.draw(timegrob(x))
}


###################################################
### code chunk number 5: simplegrob
###################################################
grid.draw(tg)


###################################################
### code chunk number 6: extensibility.Rnw:110-112
###################################################
bt <- grob(x=unit(.5, "npc"), y=unit(.5, "npc"), 
                  label="hi", name="bt", cl="boxedtext")


###################################################
### code chunk number 7: extensibility.Rnw:117-126
###################################################
boxedtext <- function(x) {
    tg <- textGrob(x$label, x$x, x$y,
                   name=paste(x$name, "text", sep=".")) 
    rg <- rectGrob(x$x, x$y, 
                   width=grobWidth(tg) + unit(2, "mm"),
                   height=grobHeight(tg) + unit(2, "mm"),
                   name=paste(x$name, "rect", sep="."))
    gTree(children=gList(tg, rg), gp=x$gp, name=x$name)
}


###################################################
### code chunk number 8: extensibility.Rnw:130-133
###################################################
drawDetails.boxedtext <- function(x, ...) {
    grid.draw(boxedtext(x))
}


###################################################
### code chunk number 9: notsimplegrob
###################################################
grid.draw(bt)


###################################################
### code chunk number 10: extensibility.Rnw:178-181
###################################################
primToDev.timegrob <- function(x, dev) {
    primToDev(timegrob(x), dev)
}


###################################################
### code chunk number 11: extensibility.Rnw:187-190
###################################################
grid.newpage()
grid.draw(tg)
gridToSVG("simpleclass.svg")


###################################################
### code chunk number 12: extensibility.Rnw:199-200
###################################################
library(XML)


###################################################
### code chunk number 13: extensibility.Rnw:203-205
###################################################
simpleclasssvg <- xmlParse("simpleclass.svg")
cat(saveXML(simpleclasssvg))


###################################################
### code chunk number 14: extensibility.Rnw:227-230
###################################################
primToDev.boxedtext <- function(x, dev) {
    primToDev(boxedtext(x), dev)
}


###################################################
### code chunk number 15: extensibility.Rnw:236-239
###################################################
grid.newpage()
grid.draw(bt)
gridToSVG("notsimpleclass.svg")


###################################################
### code chunk number 16: extensibility.Rnw:244-246
###################################################
notsimpleclasssvg <- xmlParse("notsimpleclass.svg")
cat(saveXML(notsimpleclasssvg))


###################################################
### code chunk number 17: extensibility.Rnw:280-286
###################################################
animate.timegrob <- function(x, ...) {
    tg <- timegrob(x)
    tg$animationSets <- x$animationSets
    tg$groupAnimationSets <- x$groupAnimationSets
    animate(tg, ...)
}


###################################################
### code chunk number 18: extensibility.Rnw:293-297
###################################################
grid.newpage()
grid.draw(tg)
grid.animate("tg", x=c(.3, .7))
gridToSVG("animsimpleclass.svg")


###################################################
### code chunk number 19: extensibility.Rnw:305-307
###################################################
animsimpleclasssvg <- xmlParse("animsimpleclass.svg")
cat(saveXML(animsimpleclasssvg))


###################################################
### code chunk number 20: extensibility.Rnw:325-337
###################################################
animate.boxedtext <- function(x, ...) {
    bt <- boxedtext(x)
    bt$groupAnimationSets <- x$groupAnimationSets
    animate(bt, ...)
    # Animate the children of bt 
    btrect <- getGrob(bt, "bt.rect")
    btrect$animationSets <- x$animationSets
    animate(btrect, ...)
    bttext <- getGrob(bt, "bt.text")
    bttext$animationSets <- x$animationSets
    animate(bttext, ...)
}


###################################################
### code chunk number 21: extensibility.Rnw:347-353
###################################################
grid.newpage()
grid.draw(bt)
grid.animate("bt", x=c(.3, .7))
grid.animate("bt", visibility=c("visible", "hidden"),
             begin=1, duration=0.1, group=TRUE)
gridToSVG("animnotsimpleclass.svg")


###################################################
### code chunk number 22: extensibility.Rnw:358-360
###################################################
animnotsimpleclasssvg <- xmlParse("animnotsimpleclass.svg")
cat(saveXML(animnotsimpleclasssvg))


###################################################
### code chunk number 23: extensibility.Rnw:391-397
###################################################
garnish.timegrob <- function(x, ...) {
    tg <- timegrob(x)
    tg$attributes <- x$attributes
    tg$groupAttributes <- x$groupAttributes 
    garnish(tg, ...)
}


###################################################
### code chunk number 24: extensibility.Rnw:405-409
###################################################
grid.newpage()
grid.draw(tg)
grid.garnish("tg", onmousedown="alert('ouch')")
gridToSVG("garnishsimpleclass.svg")


###################################################
### code chunk number 25: extensibility.Rnw:414-416
###################################################
garnishsimpleclasssvg <- xmlParse("garnishsimpleclass.svg")
cat(saveXML(garnishsimpleclasssvg))


###################################################
### code chunk number 26: extensibility.Rnw:423-429
###################################################
garnish.boxedtext <- function(x, ...) {
    bt <- boxedtext(x)
    bt$attributes <- x$attributes
    bt$groupAttributes <- x$groupAttributes 
    garnish(bt, ...)
}


###################################################
### code chunk number 27: extensibility.Rnw:443-449
###################################################
grid.newpage()
grid.draw(bt)
grid.garnish("bt", onmousedown="alert('ouch')")
grid.garnish("bt", onmouseover=c(bt.text="alert('watch it!')"),
             group=FALSE)
gridToSVG("garnishnotsimpleclass.svg")


###################################################
### code chunk number 28: extensibility.Rnw:454-456
###################################################
garnishnotsimpleclasssvg <- xmlParse("garnishnotsimpleclass.svg")
cat(saveXML(garnishnotsimpleclasssvg))



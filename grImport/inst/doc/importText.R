### R code from vignette source 'importText.Rnw'

###################################################
### code chunk number 1: importText.Rnw:22-30
###################################################
options(prompt="R> ")
options(continue = "+  ")
options(width = 60)
options(useFancyQuotes = FALSE)
strOptions(strict.width = TRUE)
library(grid)
library(lattice)



###################################################
### code chunk number 2: importText.Rnw:54-62
###################################################
addBBox <- function(existingPic) {
    PostScriptTrace("helloBBox.ps", "helloBBox.xml")
    helloBBox <- readPicture("helloBBox.xml")
    grid.picture(helloBBox[seq(1, 9, 2)], gp=gpar(lex=.1),
                 xscale=existingPic@summary@xscale,
                 yscale=existingPic@summary@yscale)
}



###################################################
### code chunk number 3: importText.Rnw:73-74
###################################################
library(grImport)


###################################################
### code chunk number 4: simpleimport
###################################################
PostScriptTrace("hello.ps", "hello.xml")
hello <- readPicture("hello.xml")
grid.picture(hello)


###################################################
### code chunk number 5: simpleimportsmooth
###################################################
PostScriptTrace("hello.ps", "hello-smooth.xml", setflat=0.5)
helloSmooth <- readPicture("hello-smooth.xml")
grid.picture(helloSmooth)


###################################################
### code chunk number 6: simpleimporttextsrc
###################################################
PostScriptTrace("hello.ps", "helloText.xml", charpath=FALSE)
helloText <- readPicture("helloText.xml")
grid.picture(helloText)


###################################################
### code chunk number 7: simpleimporttext
###################################################
PostScriptTrace("hello.ps", "helloText.xml", charpath=FALSE)
helloText <- readPicture("helloText.xml")
grid.picture(helloText)
addBBox(helloText)



###################################################
### code chunk number 8: simpleimporttextfontsrc
###################################################
grid.picture(helloText, gp=gpar(fontfamily="serif"))


###################################################
### code chunk number 9: simpleimporttextfont
###################################################
grid.picture(helloText, gp=gpar(fontfamily="serif"))
addBBox(helloText)



###################################################
### code chunk number 10: simpleimporttextheightsrc
###################################################
grid.picture(helloText, sizeByWidth=FALSE)


###################################################
### code chunk number 11: simpleimporttextheight
###################################################
grid.picture(helloText, sizeByWidth=FALSE)
addBBox(helloText)



###################################################
### code chunk number 12: simpleimporttextcharpossrc
###################################################
PostScriptTrace("hello.ps", "helloChar.xml", 
                charpath=FALSE, charpos=TRUE)
helloChar <- readPicture("helloChar.xml")
grid.picture(helloChar, sizeByWidth=FALSE)


###################################################
### code chunk number 13: simpleimporttextcharpos
###################################################
PostScriptTrace("hello.ps", "helloChar.xml", 
                charpath=FALSE, charpos=TRUE)
helloChar <- readPicture("helloChar.xml")
grid.picture(helloChar, sizeByWidth=FALSE)
addBBox(helloText)



###################################################
### code chunk number 14: simpleimporttextcharposwidthsrc
###################################################
grid.picture(helloChar)


###################################################
### code chunk number 15: simpleimporttextcharposwidth
###################################################
grid.picture(helloChar)
addBBox(helloText)




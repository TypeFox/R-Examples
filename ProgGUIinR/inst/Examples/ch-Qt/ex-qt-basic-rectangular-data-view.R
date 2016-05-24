
###################################################
### code chunk number 268: Widgets-MVC.Rnw:2-4
###################################################
require(qtbase)
require(MASS)


###################################################
### code chunk number 269: qt-mvc-qdfm
###################################################
model <- qdataFrameModel(mtcars)
view <- Qt$QTableView()
view$setModel(model)


###################################################
### code chunk number 270: qt-mvc-qdfm-access
###################################################
DF <- qdataFrame(model)
DF[1:3, 1:10]


###################################################
### code chunk number 271: Widgets-MVC.Rnw:72-73
###################################################
qdataFrame(model)$hpToMpg <- with(qdataFrame(model), hp / mpg)


###################################################
### code chunk number 272: qt-mvc-table-header-align
###################################################
header <- view$horizontalHeader()
header$defaultAlignment <- Qt$Qt$AlignLeft


###################################################
### code chunk number 273: qt-mvc-table-alternating
###################################################
view$alternatingRowColors <- TRUE

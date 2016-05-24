.rescalefct <- RobAStBase:::.rescalefct
.plotRescaledAxis <- RobAStBase:::.plotRescaledAxis
.makedotsP <- RobAStBase:::.makedotsP
.makedotsLowLevel <- RobAStBase:::.makedotsLowLevel
.SelectOrderData <- RobAStBase:::.SelectOrderData
### helper function to recursively evaluate list
.evalListRec <- RobAStBase:::.evalListRec


if(packageVersion("distrMod")<"2.5"){
.isUnitMatrix <- function(m){
### checks whether m is unit matrix
              m.row <- nrow(m)
              isTRUE(all.equal(m, diag(m.row), check.attributes = FALSE))
              }

.deleteDim <- function(x){
     attribs <- attributes(x)
     attribs$dim <- NULL
     attribs$dimnames <- NULL
     attributes(x) <- attribs
     x
     }

}

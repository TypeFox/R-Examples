## $Id: 00virtual.R 5596 2014-07-16 22:17:15Z bpikouni $ 

## Virtual Classes for cg package
## For Polymorphic slots

setOldClass("lm")
setOldClass("rlm")
setOldClass("survreg")
setOldClass("gls")

setClass("olsfit")
setClass("rrfit")
setClass("aftfit")
setClass("uvfit")

setClassUnion("olsfit", c("character", "lm"))
setClassUnion("rrfit", c("character", "rlm", "lm"))
setClassUnion("aftfit", c("character", "survreg"))
setClassUnion("uvfit", c("character", "gls"))

setClass("numericOrNULL")
setClassUnion("numericOrNULL", c("numeric", "NULL"))

setClass("dataframeMatrixOrNULL")
setClassUnion("dataframeMatrixOrNULL", c("data.frame", "matrix", "NULL"))

setClass("dataframeOrNULL")
setClassUnion("dataframeOrNULL", c("data.frame", "NULL"))

setClass("characterOrNULL")
setClassUnion("characterOrNULL", c("character", "NULL"))

setClass("characterOrExpression")
setClassUnion("characterOrExpression", c("character", "expression"))



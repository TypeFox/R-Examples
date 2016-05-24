tdlnre <- function (model, x, ...)  UseMethod("tdlnre")
tplnre <- function (model, q, lower.tail=FALSE, ...)  UseMethod("tplnre")
tqlnre <- function (model, p, lower.tail=FALSE, ...)  UseMethod("tqlnre")

dlnre <- function (model, x, ...)  UseMethod("dlnre")
plnre <- function (model, q, lower.tail=TRUE, ...)  UseMethod("plnre")
qlnre <- function (model, p, lower.tail=TRUE, ...)  UseMethod("qlnre")

ltdlnre <- function (model, x, base=10, log.x=FALSE, ...)  UseMethod("ltdlnre")
ldlnre <- function (model, x, base=10, log.x=FALSE, ...)  UseMethod("ldlnre")

rlnre <- function (model, n, ...)  UseMethod("rlnre")

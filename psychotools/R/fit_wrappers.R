## *.fit wrappers for new *model() functions.
RaschModel.fit <- function(y, ...) raschmodel(y = y, ...)
RSModel.fit <- function (y, ...) rsmodel(y = y, ...)
PCModel.fit <- function (y, ...) pcmodel(y = y, ...)
btReg.fit <- function(y, ...) btmodel(y = y, ...)

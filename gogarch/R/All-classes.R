##
## Class definition of orthogonal matrices
##
setClass(Class = "Orthom", representation(M = "matrix"))
##
## Class definition of initial GO-GARCH objects
##
setClass(Class = "Goinit", representation(X = "matrix", V = "matrix", P = "matrix", Dsqr = "matrix", garchf = "formula", name = "character"))
##
## Class definition of GO-GARCH objects
##
setClass(Class = "GoGARCH", representation(Z = "matrix", U = "matrix", Y = "matrix", H = "list", models = "list", estby = "character", CALL = "call"), contains = "Goinit")
##
## Class definition of GO-GARCH objects, estimated by ICA
##
setClass(Class = "Goestica", representation(ica = "list"), contains = "GoGARCH")
##
## Class definition of GO-GARCH objects, estimated by Methods of Moments
##
setClass(Class = "Goestmm", representation(weights = "numeric", Umatched = "list"), contains = "GoGARCH")
##
## Class definition of GO-GARCH objects, estimated by Maximum-Likelihood
##
setClass(Class = "Goestml", representation(opt = "list"), contains = "GoGARCH")
##
## Class definition of GO-GARCH objects, estimated by Non-linear Least-Squares
##
setClass(Class = "Goestnls", representation(nls = "list"), contains = "GoGARCH")
##
## Class definition for summary objects from GoGARCH
##
setClass(Class = "Gosum", representation(name = "character", method = "character", model = "formula", garchc = "list", Zinv = "matrix"))
##
## Class definition for predict objects from GoGARCH
##
setClass(Class = "Gopredict", representation(Hf = "list", Xf = "matrix", CGARCHF = "list"))



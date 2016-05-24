
set.seed(290875)
library("party")
if (!require("TH.data"))
    stop("cannot load package TH.data")

### get rid of the NAMESPACE
attach(asNamespace("party"))

gtctrl <- new("GlobalTestControl")
tlev <- levels(gtctrl@testtype)

data(GlaucomaM, package = "TH.data")
inp <- initVariableFrame(GlaucomaM[,-63,drop = FALSE], trafo = NULL) #, fun = rank)
resp <- initVariableFrame(GlaucomaM[,"Class",drop = FALSE], trafo = NULL, response = TRUE)
ls <- new("LearningSample", inputs = inp, responses = resp,
          weights = rep(1, inp@nobs), nobs = nrow(GlaucomaM), 
          ninputs = inp@ninputs)
tm <- ctree_memory(ls, TRUE)
ctrl <- ctree_control()
node <- .Call("R_Node", ls, ls@weights, tm, ctrl, PACKAGE = "party")
stopifnot(isequal(node[[5]][[3]], 0.059))

### and now with ranked inputs -> Wilcoxon-Mann-Whitney tests
inp <- initVariableFrame(GlaucomaM[,-63,drop = FALSE], trafo = function(data)
ptrafo(data, numeric_trafo = rank))
resp <- initVariableFrame(GlaucomaM[,"Class",drop = FALSE], trafo = NULL, response = TRUE)
ls <- new("LearningSample", inputs = inp, responses = resp,
          weights = rep(1, inp@nobs), nobs = nrow(GlaucomaM), 
          ninputs = inp@ninputs)
tm <- ctree_memory(ls, TRUE)
ctrl <- ctree_control()
node <- .Call("R_Node", ls, ls@weights, tm, ctrl, PACKAGE = "party")
stopifnot(isequal(node[[5]][[3]], 0.059))

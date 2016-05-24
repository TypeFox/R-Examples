setClass("multimodel", representation(data = "list", datasetind =
"vector", modelspec = "list", modellist = "list", modeldiffs = "list",
fit = "fit", parorder = "list", parorderdiff = "list", parorderchange
= "list", finished = "logical", groups = "list", stderrclp =
"logical", getXsuper = "logical", nnls = "logical", nnlscrit =
"list", optlist = "list", trilinear = "logical", nclp = "numeric",
                                      algorithm = "character"),

         prototype = list(data=list(), modelspec = list(),
modellist = list(), modeldiffs = list(), datasetind = vector(), fit =
fit(), parorder = list(),parorderdiff = list(),parorderchange =
list(), finished=FALSE, groups = list(), nnls = FALSE, stderrclp = FALSE,
getXsuper = FALSE, nnlscrit = list(), optlist = list(), trilinear=FALSE,
           nclp = 0, algorithm="nls"))



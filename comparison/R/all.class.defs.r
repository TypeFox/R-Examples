setClass("ece", representation(prior="numeric", ece.null="numeric", ece="numeric", ece.cal="numeric"))
setClass("compcovar", representation(v.within="matrix", v.between="matrix", n.observations="numeric", n.items="numeric", item.n="matrix", item.means="matrix", n.vars="numeric", overall.means="numeric", multivariate="logical", balanced="logical", s.within="matrix", s.between="matrix", warn.type="character"))
setClass("compitem", representation(item.means="numeric", n.replicates="numeric", n.vars="numeric", multivariate="logical", observed="matrix", warn.type="character"))



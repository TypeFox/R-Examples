enumerateCycles <- function(model,rxnList,solver = SYBIL_SETTINGS("SOLVER") ){
fva=lrFVA(model,rxnList,solver)
ldet=fva[[2]]
lp=ldet[ldet[,"llen"]>0,]
ulp=unique(lp[,c("lp","llen")])
return (ulp)
}
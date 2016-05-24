attrib.dens <-
function(optim.param)
{
	if(identical(optim.param,optim.const.ordi)|identical(optim.param,optim.noconst.ordi)) attr(optim.param,"type") <- "ordi"
	if(identical(optim.param,optim.diff.norm)|identical(optim.param,optim.equal.norm)|identical(optim.param,
					optim.gene.norm)|identical(optim.param,optim.indep.norm)) attr(optim.param,"type") <- "norm"
	optim.param
}


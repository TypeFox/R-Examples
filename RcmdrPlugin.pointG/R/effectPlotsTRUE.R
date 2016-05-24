effectPlotsTRUE <- function(){
	Library("effects")
	.activeModel <- ActiveModel()
	if (is.null(.activeModel) || !checkMethod("Effect", .activeModel)) return()
	doItAndPrint('trellis.device(theme="col.whitebg")')
	command <- paste("plot(allEffects(", .activeModel, "), ask=TRUE)", sep="")
	justDoIt(command)
	logger(command)
	activateMenus()
	NULL
}

plot.FactorNetworks <-
function (x, interactiveMode=T, ...)
{
	interactivePrompt<-"Press return for next plot..."
	for (Level in names(x))
	{
		plot(x[[Level]]$Network,Level)
		if (interactiveMode & Level!=names(x)[length(x)]){line <- readline(interactivePrompt)}
	}
	if(interactiveMode){line <- readline("Press return to close the plot...");dev.off()}
}

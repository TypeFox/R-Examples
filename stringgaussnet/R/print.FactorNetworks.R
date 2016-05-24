print.FactorNetworks <-
function (x, nlimit=2, ...)
{
	cat("Object of class FactorNetworks (package stringgaussnet)\n\n")
	Levels<-names(x)
	LevelsTable<-vector()
	for (Level in Levels){LevelsTable[Level]<-nrow(x[[Level]]$DEGeneExpr$DataExpression)}
	cat("Levels distribution:\n")
	print(LevelsTable)
	for (Level in Levels)
	{
		cat("\n",Level,":\n\n",sep="")
		print(x[[Level]]$DEGeneExpr,nlimit)
		cat("\n")
		print(x[[Level]]$Network,nlimit)
	}
}

print.summary.sla <-
function(x, ...)
{
	title.string.Pretty.Fit.Table <- ' Description of Fits for 4 ANCOVA Models'
	title.string.Pretty.Test.Table <- ' ANCOVA Tests: Two Groups/Straight Line Fits'	
	Fit.Table.Pretty <- x$Fit.Table.Pretty
	Test.Table.Pretty <- x$Test.Table.Pretty
	cat('\nCall:  ')
	print(x$Call)
	cat('\nSummary of ANCOVA Tests. . .\n')
	cat('\n', title.string.Pretty.Fit.Table, '\n\n')
	print(Fit.Table.Pretty)
	cat('\n', title.string.Pretty.Test.Table, '\n\n')
	print(Test.Table.Pretty)
	cat('\n')
}

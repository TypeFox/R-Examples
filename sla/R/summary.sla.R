summary.sla <-
function(object, ...)
{
	Call <- object$Call
	Fit.Table.Pretty <- object$Fit.Table.Pretty
	Test.Table.Pretty <- object$Test.Table.Pretty
	summary.list <- list(Call = Call, Fit.Table.Pretty = Fit.Table.Pretty, Test.Table.Pretty = Test.Table.Pretty)
	class(summary.list) <- 'summary.sla'
	summary.list
}

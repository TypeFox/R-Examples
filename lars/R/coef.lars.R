coef.lars <-
function(object, ...)
{
	predict(object, type = "coefficient", ...)$coef
}


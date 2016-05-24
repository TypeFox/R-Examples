as.AtkAttributeSet <-
function(x)
{
	x <- as.GSList(x)
	class(x) <- "AtkAttributeSet"
	x
}
as.AtkAttribute <-
function(x)
{
	attr <- as.character(x)
	names(attr) <- names(x)
	attr
}

as.AtkTextRectangle <-
function(x)
{
	x <- as.GdkRectangle(x)
	class(x) <- "AtkTextRectangle"
	x
}

as.AtkRectangle <-
function(x)
{
	x <- as.GdkRectangle(x)
	class(x) <- "AtkRectangle"
	x
}

DEGeneExpr.default <- function (x,y,Identifier=0,...)
{
	if (is.data.frame(x))
	{
		x=as.matrix(x)
	}
	if (is.matrix(y))
	{
		y=as.data.frame(y)
	}
	if (!is.matrix(x))
	{
		stop("x is not a matrix or data frame")
	}
	if (!is.data.frame(y))
	{
		stop("y is not a matrix or data frame")
	}
	if (ncol(x)!=nrow(y))
	{
		stop("The number of columns in x differs of the number of row in y")
	}
	if (Identifier==0)
	{
		Identifiers <- rownames(y)
	}
	if (Identifier!=0)
	{
		Identifiers <-y[,Identifier]
	}
	if (is.null(Identifiers) | length(which(!(colnames(x) %in% Identifiers)))>0)
	{
		stop("Problem with genes identifiers")
	}
	rownames(y)=Identifiers
	x <- x[,Identifiers]
	DEGeneExpr <- list(DataExpression=x,DEGenesResults=y)
	
	class(DEGeneExpr) <- "DEGeneExpr"
	return(DEGeneExpr)
}

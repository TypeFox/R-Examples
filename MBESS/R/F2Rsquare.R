"F2Rsquare" <-
function(F.value=NULL, df.1=NULL, df.2=NULL)
{
if(is.null(df.1) | is.null(df.2)) stop("You have not specified \'df.1\' and/or \'df.2\'.")
return(F.value*df.1/(F.value*df.1 + df.2))
}

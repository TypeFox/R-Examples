"Rsquare2F" <-
function(R2=NULL, df.1=NULL, df.2=NULL, p=NULL, N=NULL)
{
if(is.null(df.1) & is.null(df.2) & !is.null(N) & !is.null(p))
{
df.1 <- p
df.2 <- N-p-1
}
if(is.null(df.1) | is.null(df.2)) stop("You have not specified \'df.1\', \'df.2\', \'N\', and/or \'p\' correctly.")
return((R2/df.1)/((1-R2)/df.2))
}

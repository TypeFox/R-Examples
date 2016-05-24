rchi <-
function (n, df = stop("no 'df' arg")) 
{
    sqrt(rchisq(n, df))
}

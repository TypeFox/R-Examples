block_cross <-
function(A,B, sortcol,sortcol2)
{
combine <- NULL
for (i in 1:length(sortcol))
{
result <- crossprod(B[,sortcol[i]],A[,sortcol2[i]])
combine <- c(combine,result)
}
return(combine)
}

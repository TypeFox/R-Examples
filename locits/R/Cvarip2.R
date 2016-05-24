Cvarip2 <-
function (i, p, ll, S, Pmat, PsiJL) 
{
    ans <- 0
    ans <- .C("Cvarip2", i = as.integer(i), ll = as.integer(ll), 
        S = as.double(S), lS = as.integer(length(S)), Pmat = as.double(Pmat), 
        ncP = as.integer(ncol(Pmat)), nrP = as.integer(nrow(Pmat)), 
        PsiJL = as.integer(PsiJL), ans = as.double(ans),
	PACKAGE="locits")$ans
    ans
}

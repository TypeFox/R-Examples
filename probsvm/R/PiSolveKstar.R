PiSolveKstar <-
function (Kstar, vec, con, ridge) 
{
    al <- c(con, vec)
    if (ridge > 0) 
        Kstar <-  Kstar + diag(length(al)) * ridge
    solve(Kstar, al)

}

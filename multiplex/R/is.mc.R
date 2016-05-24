is.mc <-
function (B, C, A, ord = NULL) 
{
    return(!(all(dichot(transf(B, "listmat", ord = ord) + transf(C, 
        "listmat", ord = ord), c = 2) == transf(A, "listmat", 
        ord = ord))))
}

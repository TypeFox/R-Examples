transformer.bord.j <-
function(bord)
{
    m <- length(bord)
    
    bord <- bord[c(2:m, 1)]
    
    return(bord)
}

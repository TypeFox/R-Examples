transformer.bord.i <-
function(bord, part)
{
    bord <- unlist(by(bord, part, transformer.bord.j))
    
    return(bord)
}

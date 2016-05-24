transformer.bord <-
function(bord)
{
    b <- apply(bord[,c("X","Y")], MARGIN=2, transformer.bord.i, bord$id)
    bord <- cbind(b,bord)
    colnames(bord) <- c("X1","Y1","X2","Y2","id")

    return(bord)
}

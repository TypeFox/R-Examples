inter.sgmts <-
function(bord, coeff)
{
    m <- nrow(bord)
    bordX <- bord$X
    bordY <- bord$Y
    coup <- rep(0, m)
    coup.X <- rep(0, 2)
    coup.Y <- rep(0, 2)
    
    res <- .C("interSegments",
              as.integer(m), as.double(bordX), as.double(bordY), as.double(coeff),
              as.integer(coup), as.double(coup.X), as.double(coup.Y)
              
             ) #ajout ,PACKAGE="SPODT"
             
    resultats <- list(part=res[[5]], X=res[[6]], Y=res[[7]])
    
    return(resultats)
}

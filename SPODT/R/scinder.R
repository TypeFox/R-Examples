scinder <-
function(coup, bord, coeff)
{
    m <- nrow(bord)

    bordX1 <- c(bord$X[which(coup$part==2)], rev(coup$X))
    bordY1 <- c(bord$Y[which(coup$part==2)], rev(coup$Y))
    bordX2 <- c(bord$X[which(coup$part==1)], coup$X, bord$X[which(coup$part==3)])
    bordY2 <- c(bord$Y[which(coup$part==1)], coup$Y, bord$Y[which(coup$part==3)])

    if (bordY1[1] <= coeff[1] * bordX1[1] + coeff[2])
    {
        bord.g <- data.frame(X=bordX1, Y=bordY1)
        bord.d <- data.frame(X=bordX2, Y=bordY2)
    }
    else
    {
        bord.d <- data.frame(X=bordX1, Y=bordY1)
        bord.g <- data.frame(X=bordX2, Y=bordY2)
    }

    return(list(g=bord.g, d=bord.d ))
}

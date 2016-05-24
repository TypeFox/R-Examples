segments.greffe <-
function(cple, bord)
{
    ind <- which(bord$id == cple[1]  |  bord$id == cple[2])
    
    X1 <- bord$X1[ind]
    Y1 <- bord$Y1[ind]
    X2 <- bord$X2[ind]
    Y2 <- bord$Y2[ind]
    
    ind.cl <- c(length(which(bord$id == cple[1])), length(ind))

    grf <- rep(0, 4*length(which(bord$id == cple[1]))*length(which(bord$id == cple[2])))

    res <- .C("sgmtsGrf",
              as.integer(ind.cl), as.double(X1), as.double(Y1), as.double(X2), as.double(Y2),
              as.double(grf)
              
             ) #ajouté ,PACKAGE="SPODT"

    grf <- t(matrix(res[[6]], nrow=4))
    grf <- grf[which(grf[,1]!=grf[,3]  &  grf[,2]!=grf[,4]), ]


    bord$id[ind] <- rep(max(bord$id)+2, length(ind))

    bord <- rbind.data.frame(bord)[c(which(!(1:nrow(bord)%in%ind)), ind),]

    return(list(grf=grf, bord=bord))
}

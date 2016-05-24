part.vqt <-
function(data, vqt, ordre.vqt, min.fils, ponderer)
{
    num.vqt <- 0
    seuil <- 0
    vic <- 0
    ordre <- ind.ordre.vqt(ordre.vqt, rownames(data))
    
    res <- .C("partVqt",
              as.integer(nrow(data)), as.double(data$x), as.double(data$y), as.double(data$z),
              as.double(unlist(data[,vqt],use.names=FALSE)),
              as.integer(length(vqt)), as.integer(as.vector(ordre)), as.double(sum(data$z)),
              as.integer(ponderer), as.integer(min.fils),
              as.integer(num.vqt), as.double(seuil), as.double(vic)
              
             )#ajout ,PACKAGE="SPODT"

    if(res[[13]] != 0)
    {
        vqt <- vqt[res[[11]]]
        partition <- rep(0, nrow(data))
        partition[which(data[,vqt] <= res[[12]])] <- -1
        partition[which(data[,vqt] > res[[12]])] <- 1
    }
    else
    {
        vqt <- 0
        partition <- 0
    }

    return(list(vrbl=vqt, seuil=res[[12]], vic=res[[13]], part=partition))
}

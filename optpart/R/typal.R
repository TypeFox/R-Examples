typal <- function (clustering,dist,k=1)
{
    clustering <- as.integer(clustify(clustering))
    if (class(dist) != "dist" ) {
        stop("typal is not defined for classes other than dist")
    }
    classes <- 1:length(table(clustering))

    part <- partana(clustering,dist)
    sil <- silhouette(clustering,dist)

    part.out <- matrix(NA,nrow=max(classes),ncol=k)
    for (i in classes) {
        tmp <- clustering==i
        names <- part$names[tmp]
        vals <- part$ptc[tmp,i]
        part.out[i,] <- names[rev(order(vals))][1:k]
    }
    part.out <- data.frame(part.out)
    names(part.out) <- as.character(1:k)

    sil.out <- matrix(NA,nrow=max(classes),ncol=k)
    for (i in classes) {
        tmp <- clustering==i
        names <- attr(dist,'Labels')[tmp]
        vals <- sil[tmp,3]
        sil.out[i,] <- names[rev(order(vals))][1:k]
    }
    sil.out <- data.frame(sil.out)
    names(sil.out) <- as.character(1:k)

    out <- list(partana=part.out,silhouette=sil.out)
    out
}

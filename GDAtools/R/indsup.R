indsup <- function(resmca,supdata) {
    z <- dichotom(supdata,out='numeric')
    z <- as.matrix(z)
    classe <- class(resmca)[1] # new
    if(classe %in% c('speMCA','csMCA')) z <- z[,-resmca$call$excl] # new
    Q <- ncol(supdata)
    delta <- 1/sqrt(resmca$eig$eigen[1:resmca$call$ncp])
    #pk <- resmca$call$marge.col
    vcoord <- resmca$var$coord
    #pkvcoord <- sweep(vcoord,1,pk,'*')
    coord <- (1/Q)*z%*%vcoord
    coord <- sweep(coord,2,delta,'*')
    #coord <- z%*%vcoord*(1/Q)-z%*%pkvcoord
    GM2 <- rowSums(coord^2)
    cos2 <- sweep(coord^2,1,GM2,'/')
    cos2 <- round(cos2,6)
    rownames(coord) <- rownames(supdata)
    rownames(cos2) <- rownames(supdata)
    return(list(coord=coord,cos2=cos2))
    }
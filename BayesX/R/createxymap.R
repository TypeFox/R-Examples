createxymap <- function(x, y, districts=NULL, p=2, max.dist){

    if(length(x)!=length(y))
        stop("Lengths of x and y differ!")

    if(p < 1)
        stop("p have to be p >= 1!")

    S <- length(x)
    if(is.null(districts))
        districts <- paste(seq(1,S))

    ## gra format
    gra <- matrix(0,S,S)

    if(p==Inf){
        for(i in 1:(S-1)){
            for(j in (i+1):S){
                xdiff <- x[j] - x[i]
                ydiff <- y[j] - y[i]

                if(max(abs(xdiff),abs(ydiff)) <= max.dist){
                    gra[i,j] <- gra[j,i] <- -1
                }
            }
        }
    }

    else{
        for(i in 1:(S-1)){
            for(j in (i+1):S){
                xdiff <- x[j] - x[i]
                ydiff <- y[j] - y[i]

                if(((abs(xdiff))^p+(abs(ydiff))^p)^(1/p) <= max.dist){
                    gra[i,j] <- gra[j,i] <- -1
                }
            }
        }
    }

    diag(gra) <- -apply(gra, 1, sum)

    rownames(gra) <- districts
    colnames(gra) <- districts

    class(gra) <- "gra"

    ## bnd format
    x.range <- max(x) - min(x)
    y.range <- max(y) - min(y)
    step <- max(x.range,y.range)/100
    x.range <- x.range + 2*step
    y.range <- y.range + 2*step
    height2width <- round(y.range/x.range, digits=2)

    bnd <- list()
    for(i in 1:S){
        bnd[[i]] <- matrix(data=c(x[i]-step,y[i]-step,
                           x[i]-step,y[i]+step,
                           x[i]+step,y[i]+step,
                           x[i]+step,y[i]-step,
                           x[i]-step,y[i]-step),byrow=T,nrow=5,ncol=2)
    }

    names(bnd) <- districts

    attr(bnd, "surrounding") <- replicate(n=length(bnd),
                                          expr=character())
    attr(bnd, "height2width") <- height2width
    class(bnd) <- "bnd"

    maps <- list(gra,bnd)
    names(maps) <- c("gra","bnd")
    return(maps)
}

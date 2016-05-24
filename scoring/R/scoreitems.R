scoreitems <- function(param, data, fam, ordered, decomp=FALSE, group=NULL){
    if(ordered){
        item.res <- ordwrap(data, param, fam)
    } else {
        scorefun <- paste(fam,"fam",sep="")

        nalts <- ncol(data)-1

        item.res <- apply(data, 1, function(x){
            do.call(scorefun, list(x[1:nalts], x[length(x)], param, fam))
        })
    }

    if(decomp){
        if(nalts>2) stop("Brier score decompositions are only valid for 2-alternative forecasts.\n")
        d.res <- bdecomp(data[,1:nalts], data[,ncol(data)], group)

        item.res <- list(rawscores=item.res, decomp=d.res)
    }
    
    item.res
}

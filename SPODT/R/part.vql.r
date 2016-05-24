part.vql <- function(data, vql, min.fils, ponderer)
{
    if (!is.null(dim(data[,vql])))
    {
        tab.vql <- apply(data[,vql], MARGIN=2, table)
        n.mod <- unlist(lapply(tab.vql, length))
        vql.ok <- names(n.mod[which(n.mod > 1)])
        vql.sortie <- names(n.mod[which(n.mod <= 1)])
    }
    else
    {
        tab.vql <- table(as.vector(data[,vql]))
        n.mod <- length(tab.vql)
        if (n.mod < 1)
        {
            return(list(vic=0, vql.s=vql))
        }
        vql.ok <- vql
        vql.sortie <- character(0)
    }

    if (length(vql.ok) < 1)
    {
        return(list(vic=0, vql.s=vql.sortie))
    }
    
    if(!is.vector(vql.ok)){
        tab.vql <- tab.vql[vql.ok]
        n.mod <- n.mod[vql.ok]
    }

    n.vql <- length(vql.ok)
    somz <- sum(data$z)

    if (!is.null(dim(data[,vql.ok])))
    {
        z.mod <- unlist(apply(data[,vql.ok], MARGIN=2, som.vql, data.som=data$z), use.names=FALSE)
        x.mod <- unlist(apply(data[,vql.ok], MARGIN=2, som.vql, data.som=data$x), use.names=FALSE)
        y.mod <- unlist(apply(data[,vql.ok], MARGIN=2, som.vql, data.som=data$y), use.names=FALSE)
        x2.mod <- unlist(apply(data[,vql.ok], MARGIN=2, som2.vql, data.som2=data$x), use.names=FALSE)
        y2.mod <- unlist(apply(data[,vql.ok], MARGIN=2, som2.vql, data.som2=data$y), use.names=FALSE)
        xy.mod <- unlist(apply(data[,vql.ok], MARGIN=2, som.vql, data.s=data$x*data$y), use.names=FALSE)
        eff.mod <- unlist(apply(data[,vql.ok], MARGIN=2, table), use.names=FALSE)
    }
    else
    {
        fact <- as.vector(data[,vql.ok])
        z.mod <- unlist(som.vql(fact, data.som=data$z), use.names=FALSE)
        x.mod <- unlist(som.vql(fact, data.som=data$x), use.names=FALSE)
        y.mod <- unlist(som.vql(fact, data.som=data$y), use.names=FALSE)
        x2.mod <- unlist(som2.vql(fact, data.som2=data$x), use.names=FALSE)
        y2.mod <- unlist(som2.vql(fact, data.som2=data$y), use.names=FALSE)
        xy.mod <- unlist(som.vql(fact, data.som=data$x*data$y), use.names=FALSE)
        eff.mod <- unlist(table(fact))        
    }

    num.vql <- 0
    mod <- 0
    vic <- 0
          
    res <- .C("partVql",
              as.integer(nrow(data)), as.double(x.mod), as.double(y.mod),
              as.double(x2.mod), as.double(y2.mod), as.double(xy.mod),
              as.integer(n.vql), as.integer(n.mod),
              as.double(z.mod), as.double(somz), as.integer(eff.mod),
              as.integer(ponderer), as.integer(min.fils),
              as.integer(num.vql), as.integer(mod), as.double(vic)
             )
    
    if (res[[16]] != 0)
    {
        vql.opt <- vql.ok[res[[14]]]
        if(!is.vector(vql.ok)){
            mod.opt <- unlist(dimnames(tab.vql[[ res[[14]] ]]))[ res[[15]] ]
        } else{
            mod.opt <- dimnames(tab.vql)[[ res[[14]] ]][ res[[15]] ]
        }
        partition <- rep(0, nrow(data))
        partition[which(data[,vql.opt] == mod.opt)] <- -1
        partition[which(data[,vql.opt] != mod.opt)] <- 1
    }
    else
    {
        vql.opt <- 0
        mod.opt <- 0
        partition <- 0
    }

    return(list(vrbl=vql.opt, mod=mod.opt, vic=res[[16]], part=partition, vql.s=vql.sortie))
}

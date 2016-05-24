test.spodt <- function(formula, data, R2.obs,  rdist, par.rdist, nb.sim, weight=FALSE, graft=0, level.max=5, min.parent=10, min.child=5, rtwo.min=0.001)
{
    if (class(data)!="SpatialPointsDataFrame") stop("use a SpatialPointsDataFrame")
    if (is.na((is.projected(data)))|(! is.projected(data)) ) warning("the coordinates are not projected. Please, provide projected coordinates or be sure to use euclidian coordinates!")

    coord.x    <- coordinates(data)[,1]
    coord.y    <- coordinates(data)[,2]
    loc.data   <- as.numeric(row.names(data@data))
    data.temp  <- data@data

    Call <- match.call()
    indx <- match(c("formula", "data"), names(Call), nomatch = 0L)
    if (indx[1] == 0L) stop("a 'formula' with the cofactors is required\n for single spatial analysis (with no cofactor) the right hand side should be z~1")

    dataset.prep <- model.frame(Call$formula, data=data.temp)

    dataset <- cbind(loc.data, coord.x, coord.y, dataset.prep)
    colnames(dataset)[1:4] <- c("loc", "x", "y", "z")
    n <- nrow(dataset)
    rownames(dataset) <- 1:n

    mv<-data.frame(dataset[,-c(1:4)])
    colnames(mv)<-colnames(dataset)[-c(1:4)]
    L <- dim(mv)[2]
    ind_num <- 0
    ind_fact <- 0

    if(L ==0){
        vqt<-NULL
        vql<-NULL
    }
    if(L != 0){
        for (i in 1:L){
            if (class(mv[,i])=="numeric"){
                ind_num <- c(ind_num, i)
            }
            else if(class(mv[,i])=="factor"){
                ind_fact <- c(ind_fact, i)
            }
            else{
                print(paste("class of ", colnames(mv)[i], " must be numeric or factor", sep=""))
            }
        }
    ind_num  <- ind_num[-1]
    ind_fact <- ind_fact[-1]
    vqt <- colnames(mv)[ind_num]
    vql <- colnames(mv)[ind_fact]
    }

    tloi<-get(rdist)
    if (rdist=="rbinom" | rdist=="runif" | rdist=="rnorm")
    {
        data.sim <- data.frame(replicate(nb.sim, tloi(par.rdist[1],par.rdist[2],par.rdist[3])))
    }
    else if  (rdist=="rpois")
    {
        data.sim <- data.frame(replicate(nb.sim, tloi(par.rdist[1],par.rdist[2])))
    }
    else if  (rdist=="rnbinom")
    {
        data.sim <- data.frame(replicate(nb.sim, tloi(par.rdist[1],par.rdist[2],par.rdist[3],par.rdist[4])))
    }
    else
    {
        warning("wrong distribution")
    }
    R2.sim <- apply(data.sim, MARGIN=2, simuler.spodt, data=dataset,
                    vqt=vqt, vql=vql, weight=weight, graft=graft,
                    level.max=level.max, min.parent=min.parent, min.child=min.child, rtwo.min=rtwo.min)

    hist(R2.sim, xlim=c(0,1))
    abline(v=R2.obs, col="red")
   
    return(list(R2.sim=R2.sim, p=length(R2.sim[which(R2.sim > R2.obs)])/(length(R2.sim))))
}

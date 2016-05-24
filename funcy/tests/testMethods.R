
##regular data set**********************************************

if(1){
    ##detach("package:funcy",
    ##unload=TRUE)
    start <- proc.time()
    parallel=TRUE
    library(funcy)

    ##regular data sets******************************************
    k <- 5

    set.seed(3804)
    ds <- sampleFuncy(obsNr=100, timeNr=20, reg=TRUE, k=k, sd=.4)
    plotFuncy(ds@data)
    
    ##try to cluster
    methods <- "ALL"

    res1 <- funcit(methods=methods, data=ds@data, k=k,
                                   clusters=ds@clusters, seed=2405,
                   save.data=TRUE, parallel=parallel)

    plot(res1)

    ##irregular data set which is made regular*******************
    set.seed(3805)
    ds <- sampleFuncy(reg=FALSE, timeNrMin=5, timeNrMax=10, k=k, sd=0.3)
    plotFuncy(ds)
    
    ##make it reg
    data <- regFuncy(ds@data, baseType="splines", timeNr=10, nbasis=5,
                    method="interpolate")$data
    data <- regFuncy(ds@data, baseType="splines", timeNr=10, nbasis=10, method="project")$data
    ##data <- regFuncy(ds@data, method="pace", timeNr=10)$data
    
    res2 <- funcit(methods=methods, seed=2506, data=data, k=4)


    ##irregular data set*****************************************
    k <- 4
    set.seed(3806)
    ds <- sampleFuncy(timeNrMin=5, timeNrMax=10, reg=FALSE, k=k, sd=0.3)

    a <- list(coeffsCalc="estimate", average=TRUE)
    fpcCtrl <- as(a, "fpcCtrl")

    b <- list(maxit=5, baseType="eigenbasis", flexDim=TRUE)
    funcyCtrlMbc <- as(b, "funcyCtrlMbc")
    ##funcyCtrlMbc@baseType <- "splines"

    res3 <- funcit(methods=c(1,2,3),data=ds@data,k=k,
                   clusters=ds@clusters,
                   fpcCtrl=fpcCtrl,
                   funcyCtrl=funcyCtrlMbc,
                   save.data=TRUE,
                   parallel=parallel)

    plot(res1, select=1, type="all")
    plot(res1, select=c(1,2), type="centers")
    plot(res1, select=1, type="dist2centers")

    ##special plot functions for fitclust-object
    plot(res1,select="fitfclust", type="discrim")
    plot(res1, select="fitfclust", type="conf")

    ##special plot functions for FSCM-object
    plot(res1,select="fscm", type="overview")
    plot(res1,select="fscm", type="deviations")
    plot(res1, select="fscm", type="locations")


    plot(res3, type="all")
    plot(res3, type="accordance")
    plot(res3, select=3, type="fpc")
    b <- proc.time()-start
}

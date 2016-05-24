#
#  Copyright (C) 2011-2015 Christina Yassouridis
#  
#

setGeneric("plotFuncy", function(data, regTime=NULL, col=NULL,
                                 ctr=NULL, ctrOnly=FALSE,
                                 ctrCols=NULL, showLegend=TRUE,
                                 legendPlace="bottomleft", lty=3,
                                 lwd=NULL, xlim=NULL, ylim=NULL, ...)
    standardGeneric("plotFuncy"))

setMethod("plotFuncy", signature(data="matrix"),
          function(data, ...){
               format <- checkFormat(data, reformat=FALSE)$format
               if(format!="Format1")
                   data <- formatFuncy(data=data, regTime=regTime,
                                       format="Format1")
               time <- data[,3]
               curve <- data[,1]
               evals <- data[,2]
               t_unique <- sort(unique(time))
               
               if(!is.null(col)) 
                   col = col
               else
                   col = rep(1, max(curve))
               
               if(is.numeric(col) & is.null(ctrCols))
                   ctrCols <- unique(sort(col))
               
               if(length(lwd)==1 | is.null(lwd))
                   lwd <- rep(lwd, max(curve))
               
               add=FALSE
               if(!ctrOnly){
                   if(is.null(xlim))
                       xlim <- c(min(time), max(time))
                   if(is.null(ylim))
                      ylim <- c(min(evals,ctr),max(ctr,evals))
                   plot(time[curve==1], evals[curve==1],
                        xlim=xlim,
                        ylim=ylim,
                        type='l',
                        col=col[1], xlab="time", ylab="", lty=lty, lwd=lwd[1], ...)
                   for(i in 1:max(curve))
                       lines(time[curve==i],evals[curve==i], col=col[i],
                             lwd=lwd[i],lty=lty,  ...)
                   add=TRUE
               }
               
               if(!is.null(ctr))
                   matplot(t_unique, ctr, add=add, col=ctrCols, lty=1,
                           lwd=2, type='l', ...)
               
               if(showLegend)
                   legend(legendPlace, paste("cl",match(ctrCols,ctrCols)),text.col=ctrCols, horiz=FALSE)
               
           }
           )

setMethod("plotFuncy", signature(data="sampleFuncy"),
          function(data, showLegend=TRUE, legendPlace="bottomleft", lty=3, lwd=NULL, ...){
              dat <- data@data
              clusters <- data@clusters
              plotFuncy(data=dat, col=clusters, showLegend=showLegend, legendPlace=legendPlace, lty=lty, lwd=lwd,...)
          }
          )


setMethod("plot", signature(x="funcyOut", y="missing"),
          function(x, y, data, type="all", showLegend=TRUE,
                   legendPlace="bottomleft", main, ...){
              methodName <- x@methodName
              cls <- col <- x@cluster
              time <- x@time
              ctr <- x@centers
              nc <- length(cls)
              ncl <- max(cls)
              if(missing(main))
                  main <- methodName
              
              if(type=="all" ){
                  plotFuncy(data=data, regTime=time, col=cls, ctr=ctr,
                            showLegend=showLegend,
                            legendPlace=legendPlace, main=main, ...)
                  
              }else if(type=="centers"){
                  plotFuncy(data=data, regTime=time, col=cls, ctr=ctr,
                            showLegend=showLegend,
                            ctrOnly=TRUE, legendPlace=legendPlace,
                            main=main,...)
                  
              }else if(type=="shadow"){
                  sh <- new("shadow")
                  sh@cluster <- as.integer(cls)
                  sh@values <- (2*x@cldist[,1])/rowSums(x@cldist[,1:2])
                  sh@size <- as.integer(table(cls))
                  sh@k <- length(sh@size)
                  sh@similarity <- TRUE
                  plot(sh,...)
                  
              }else if(type=="dist2centers"){
                  format <- checkFormat(data, reformat=FALSE)$format
                  if(format!="Format1")
                      data <- formatFuncy(data=data, format="Format1")
                  if(missing(main))
                      main <- paste("cluster",1:ncl)
                  par(mar=par("mar"))
                  col=rgb(1,0,0,alpha=0.4)
                  ctr <- x@centers
                  ctrDist <- x@dist2centers
                  ctrDist <- (max(ctrDist)-ctrDist)/max(ctrDist)+0.01
                  par(mfrow=squareGrid(x=ncl, round=TRUE))
                  for(i in 1:ncl){
                      indx <- sum(cls==i)
                      subdat <- data[data[,1]%in%which(cls==i),]
                      subN <- table(subdat[,1])
                      nrSubcv <- length(subN)
                      subdat[,1] <- rep(1:nrSubcv, subN)
                      plotFuncy(data=subdat,
                                lwd=ctrDist[,i],
                                col=rep(col,indx),
                                lty=1,
                                main=main[i], showLegend=FALSE,...)
                      lines(ctr[,i],col=1, type="l",lwd=2, lty=2)
                  }
                  par(mfrow=c(1,1))
                  
              }else if(type=="fpc"){
                  if(is.null(x@plotParams))
                      stop("You must specify eigenbasis as baseType to use this plot.")
                  else
                      plotFPCs(x@plotParams, legendPlace=legendPlace)
              }else
                  stop(paste("The plot type",type,
                             "does not exist for the chosen algorithm."))
          }
          )


setMethod("plot", signature(x="funcyOutMbc-fitfclust", y="missing"),
          function(x, y, data, newdata=NULL, type="all", showLegend=TRUE,
                   legendPlace="bottomleft", main, ...){
              fit <- x@fit
              chf <- checkFormat(data, reformat=TRUE)
              reg <- chf$reg
              data <- chf$data
              
              if(reg)
                  plotCurvesFct <- fitfclust.plotcurves
              else
                  plotCurvesFct <- fitfclust.plotcurvesIrreg
              discrimFct <- fitfclust.discrim
              
              if(type== "discrim")
                  plotFct <- discrimFct
              
              else if(type=="conf")
                  plotFct <- plotCurvesFct
              else{
                  x <- as(x, "funcyOut")
                  plot(x=x, data=data, type=type,
                       showLegend=showLegend, legendPlace=legendPlace,
                       main=main,...)
                  return()
              }
              plotFct(fit=fit, ...)
          }
          )

setMethod("plot", signature(x="funcyOutMbc-fscm", y="missing"),
          function(x, y, data, type="all", showLegend=TRUE,
                   legendPlace="bottomleft", main, ...){
              
              if(type== "overview")
                  plotFct <- plotOverview
              else if(type=="deviations")
                  plotFct <- plotDeviations
              else if(type=="locations")
                  plotFct <- plotLoc
              else{
                  x <- as(x, "funcyOut")
                  plot(x, data=data, type=type, showLegend=showLegend,
                       legendPlace=legendPlace, main=main, ...)
                  return()
              }
              plotFct(object=x, showLegend=showLegend, ...)
          }
          )

setMethod("plot", signature(x="funcyOutList", y="missing"),
          function(x, y, data=NULL, select=NULL, type="all", showLegend=TRUE,
                   legendPlace="bottomleft", main, ...){
              
              if((type=="shadow"| type=="dist2centers") &
                 length(select)>1){
                  warning("Only the first model is used.")
                  select <- select[1]
              }
                  
              data <- x@data
              reg <- x@reg
              if(is.null(data))
                  stop("You must choose option save.data=TRUE to use plot.")
              
              if(type=="accordance"){
                  accordance <- x@accordance
                  votedCluster <- x@votedCluster
                  res <- colorFct(cluster=votedCluster,
                                  accordance=accordance,
                                  method="accordance")
                  if(missing(main))
                      main <- "Accordance plot"
                  plotFuncy(data=data,
                            col=res$cols,
                            ctrCols=res$ctrCols,
                            main=main,
                            showLegend=showLegend,
                            legendPlace=legendPlace,
                            lty=1,
                            lwd=accordance, ...)
                  
              }else{
                  if(!is.null(select))
                      models <- select
                  else
                      models <- seq.int(length(x@models))
                  n <- length(models)
                  if(missing(main))
                      main <- x@methodName
                  par(mfrow=squareGrid(x=n, round=TRUE))
                  par(mar=c(2,2,3,2))
                  for(i in models)
                      plot(x[[i]], data=data, type=type, legendPlace=legendPlace,
                           showLegend=showLegend, main=main[i],  ...)
                  par(mfrow=c(1,1))
                  
              }
              
          }
          )


colorFct <- function(cluster, dist2centers, accordance=NULL,
                    method="dist"){
    nc <- length(cluster)
    nrCl <- max(cluster)
    
    if(method=="dist"){
        intensity <- sapply(1:nc, function(x) dist2centers[x,cluster[x]])
        intensity <- (max(intensity)-intensity)/(max(intensity))
    }else if(method=="accordance")
         intensity <- accordance/max(accordance)
    
    cols <- rep(0,nc)
    hue <- seq(1,360,length.out=nrCl+1)[1:nrCl]-1
    ctrCols <- hcl(hue,l=50)
    
    for(i in 1:nrCl){
        indx <- which(cluster==i)
        tempIntensity <- intensity[indx]
        tempN <- length(indx)
        cols[indx] <- hcl(hue[i], c=tempIntensity*100, alpha=0.5)
    }
    
    return(list(cols=cols, ctrCols=ctrCols, intensity=intensity))
}




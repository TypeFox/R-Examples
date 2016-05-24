####################################################################
## clValid Methods
####################################################################

####################################################################
## Accessor Functions
####################################################################

## cluster methods accessor
setGeneric("clusterMethods", function(object, ...) standardGeneric("clusterMethods"))
setMethod("clusterMethods",signature(object="clValid"),
          function(object) return(object@clMethods))

## number of clusters accessor
setGeneric("nClusters", function(object, ...) standardGeneric("nClusters"))
setMethod("nClusters",signature(object="clValid"),
          function(object) return(object@nClust))

## measure names accessor
setGeneric("measNames", function(object, ...) standardGeneric("measNames"))
setMethod("measNames",signature(object="clValid"),
          function(object) return(object@measNames))

## clusters accessor
setGeneric("clusters", function(object, ...) standardGeneric("clusters"))
setMethod("clusters",signature(object="clValid"),
          function(object,method=clusterMethods(object)) {
            method <- match.arg(method,clusterMethods(object)) ##, several.ok=TRUE)
            return(object@clusterObjs[[method]])})

## measures accessor
setGeneric("measures", function(object, ...) standardGeneric("measures"))
setMethod("measures",signature(object="clValid"),
          function(object,measures=measNames(object)) {
            measures <- match.arg(measures,measNames(object),several.ok=TRUE)
            return(object@measures[measures,,,drop=FALSE])})


####################################################################
## Print and Show Methods
####################################################################

setMethod("print","clValid",
          function(x) {
            cat("\nCall:\n")
            print(x@call); cat("\n")
            cat("Clustering Methods:\n",clusterMethods(x),"\n\n")
            cat("Cluster sizes:\n",nClusters(x),"\n\n")
            cat("Validation measures:\n",measNames(x),"\n\n")
          })

setMethod("show","clValid",
          function(object) {
            cat("\nCall:\n")
            print(object@call); cat("\n")
            cat("Clustering Methods:\n",clusterMethods(object),"\n\n")
            cat("Cluster sizes:\n",nClusters(object),"\n\n")
            cat("Validation measures:\n",measNames(object),"\n\n")
          })

####################################################################
## Summary Method
####################################################################

setMethod("summary","clValid",
          function(object, digits = max(3,getOption("digits")-3)) {
            cat("\nClustering Methods:\n",clusterMethods(object),"\n\n")
            cat("Cluster sizes:\n",nClusters(object),"\n\n")
            cat("Validation Measures:\n")
            print(ftable(round(measures(object),digits),row.vars=c(3,1)))
            cat("\n")
            ## Find best scores
            ## APN, AD, ADM, Connectivity, FOM minimized
            ## BHI, BSI, Dunn, Silhouette maximized
            measNames <- measNames(object)
            best <- numeric(length(measNames))
            bestMeth <- character(length(measNames))
            bestNc <- character(length(measNames))
            names(best) <- names(bestMeth) <- names(bestNc) <- measNames
            minmeas <- c("APN", "AD", "ADM", "FOM", "Connectivity")
            maxmeas <- c("BHI","BSI","Dunn","Silhouette")
            ## Measures to minimize
            if (any(a <- minmeas%in%measNames)) {
              best[minmeas[a]] <- apply(measures(object)[minmeas[a],,,drop=FALSE],1,min,na.rm=TRUE)
              bestInd <- apply(measures(object)[minmeas[a],,,drop=FALSE],1,function(x) which(x==min(x,na.rm=TRUE),arr.ind=TRUE)[1,])
              bestNc[minmeas[a]] <- nClusters(object)[bestInd[1,]]
              bestMeth[minmeas[a]] <- clusterMethods(object)[bestInd[2,]]
            }
            ## Measures to maximize
            if (any(a <- maxmeas%in%measNames)) {
              best[maxmeas[a]] <- apply(measures(object)[maxmeas[a],,,drop=FALSE],1,max,na.rm=TRUE)
              bestInd <- apply(measures(object)[maxmeas[a],,,drop=FALSE],1,function(x) which(x==max(x,na.rm=TRUE),arr.ind=TRUE)[1,])
              bestNc[maxmeas[a]] <- nClusters(object)[bestInd[1,]]
              bestMeth[maxmeas[a]] <- clusterMethods(object)[bestInd[2,]]
            }
            
            cat("Optimal Scores:\n\n") 
            print(data.frame("Score"=round(best,digits),"Method"=bestMeth,"Clusters"=bestNc), right=FALSE)
            cat("\n")
          })



####################################################################
## Plot Method
####################################################################

setMethod("plot",c("clValid","missing"),
          function(x,y,measures=measNames(x), legend=TRUE, legendLoc="topright", main=NULL,
                   pch=NULL, type="b", ask=prod(par("mfcol")) < length(measures) && dev.interactive(), ...) {
            measures <- match.arg(measures,measNames(x),several.ok=TRUE)
            methods <- clusterMethods(x)
            nclust <- nClusters(x)
            k <- length(methods)
            if (ask) {
              op <- par(ask = TRUE)
              on.exit(par(op))
            }
            ##        if(is.null(main))
            ##          main <- paste("Validation Measures for ", deparse(substitute(x, sys.frame(-1))))
            if (is.null(pch)) 
              pch <- c(paste(c(1:9, 0)), letters)[1:k]            
            for(i in 1:length(measures)) {
              if (is.null(main)) {
                main <- switch(measures[i],
                               APN="Stability validation",
                               AD="Stability validation",
                               ADM="Stability validation",
                               FOM="Stability validation",
                               Connectivity="Internal validation",
                               Dunn="Internal validation",
                               Silhouette="Internal validation",
                               BHI="Biological validation",
                               BSI="Biological validation")
              }
              matplot(measures(x)[measures[i],,],type=type,ylab=measures[i],
                      xlab="Number of Clusters",col=1:k,
                      lty=1:k,main=main,xaxt="n", pch=pch, ...)
              axis(1,at=1:length(nclust),labels=nclust)
              if(legend) legend(x=legendLoc,methods,lty=1:k,col=1:k, pch=pch,...)
            }
          })


######################################################################
## Some other utility functions
######################################################################

## optimalScores method
setGeneric("optimalScores",function(object, ...) standardGeneric("optimalScores"))
setMethod("optimalScores", signature(object="clValid"),
          function(object,measures=measNames(object)) {
            measNames <- match.arg(measures, measNames(object), several.ok=TRUE)
            best <- numeric(length(measNames))
            bestMeth <- character(length(measNames))
            bestNc <- character(length(measNames))
            names(best) <- names(bestMeth) <- names(bestNc) <- measNames
            minmeas <- c("APN", "AD", "ADM", "FOM", "Connectivity")
            maxmeas <- c("BHI","BSI","Dunn","Silhouette")
            ## Measures to minimize
            if (any(a <- minmeas%in%measNames)) {
              best[minmeas[a]] <- apply(measures(object)[minmeas[a],,,drop=FALSE],1,min)
              bestInd <- apply(measures(object)[minmeas[a],,,drop=FALSE],1,function(x) which(x==min(x),arr.ind=TRUE)[1,])
              bestNc[minmeas[a]] <- nClusters(object)[bestInd[1,]]
              bestMeth[minmeas[a]] <- clusterMethods(object)[bestInd[2,]]
            }
            ## Measures to maximize
            if (any(a <- maxmeas%in%measNames)) {
              best[maxmeas[a]] <- apply(measures(object)[maxmeas[a],,,drop=FALSE],1,max)
              bestInd <- apply(measures(object)[maxmeas[a],,,drop=FALSE],1,function(x) which(x==max(x),arr.ind=TRUE)[1,])
              bestNc[maxmeas[a]] <- nClusters(object)[bestInd[1,]]
              bestMeth[maxmeas[a]] <- clusterMethods(object)[bestInd[2,]]
            }
            return(data.frame("Score"=best,"Method"=bestMeth,"Clusters"=bestNc))
          })

####################################################################


#
#  Copyright (C) 2011-2015 Christina Yassouridis
#  
#

funcit <- function(data,  k, 
                   methods=c("fitfclust", "distclust", "iterSubspace", "funclust",
                       "funHDDC", "fscm", "waveclust"),
                   seed=NULL,
                   regTime=NULL, clusters=NULL,
                   funcyCtrl=NULL, fpcCtrl=NULL,
                   parallel=FALSE, save.data=TRUE, ...){
    
    ##Check for missing arguments
    if(missing(methods))
        stop("Please select one method or methods='ALL'.")
    else if(length(methods)==1)
        if(methods=="ALL")
            methods <- 1:7
    ##Method names
    allMethods <- c("fitfclust", "distclust", "iterSubspace", "funclust",
                    "funHDDC", "fscm", "waveclust")
    
    if(is.numeric(methods))
        usedMethods <- allMethods[methods]
    else
        usedMethods <- match.arg(methods, allMethods, several.ok=TRUE)
    nrMethods <- length(usedMethods)
    
    if(missing(k))
        stop(paste(k, "is missing"))
    if(!(class(data) %in% c("matrix", "data.frame")))
        stop(paste(data,
                   "must be given in matrix or data.frame format."))
    
    ##Check if data is in the right format
    chf <- checkFormat(data)
    data <- chf$data
    reg <- chf$reg

    ##Check if funcyCtrl class is given
    if(is.null(funcyCtrl))
        funcyCtrl <- new("funcyCtrl")
    funcyCtrl@seed <- seed

    ##Convert funcyCtrl automatically to funcyCtrl if model based cluster
    ##algorithm but fpcCtrl was chosen 
    if(sum(usedMethods%in%allMethods[c(1,3,4,5,6,7)])>0 & class(funcyCtrl)=="funcyCtrl")
        funcyCtrl <- as(funcyCtrl, "funcyCtrlMbc")
    
    ##Check if fpcCtrl object if defined for eigenbasis
    if(funcyCtrl@baseType!="eigenbasis" & !is.null(fpcCtrl))
        warning("fpcCtrl is ignored since it controls only eigenbasis.")
    else if(funcyCtrl@baseType=="eigenbasis")
        fpcCtrl <- fpcCtrlCheck(fpcCtrl=fpcCtrl, data=data, reg=reg)
    
    ##Check if correct method for the given dataset was chosen.
    if(reg==0 & sum(usedMethods%in%c("fscm", "funclust", "funHDDC"))>0){
        notWork <- usedMethods[which(usedMethods%in%allMethods[-c(1:3)])]
        stop(paste("Algorithm", notWork,
                   "works only on regular data!\n Please choose one of fitfclust, distclust or iterSubspace."))}
    
    
    if(.Platform$OS.type!="unix"){
        parallel <- FALSE
        warning("Parallel computing is only supported on Unix platforms.")
        mcparallel <- identity
    }else if(.Platform$OS.type=="unix" & parallel==FALSE){
        mcparallel <- identity
    }else if(.Platform$OS.type=="unix" & parallel){
        coresNr <- detectCores()-1
        options("cores"=coresNr)
    }
      
    parallelFct <- mcparallel

    
    RES <-  list()
    ##Method1--------------------
    if("fitfclust" %in% usedMethods){
        indx <- match("fitfclust",usedMethods)
        RES[[indx]] <-
            parallelFct(fitfclustWrapper(data=data, k=k, 
                                         reg=reg, regTime=regTime, fpcCtrl=fpcCtrl,
                                         funcyCtrlMbc=funcyCtrl,
                                         ...))
    }
    ##Method2----------------------
    if("distclust" %in% usedMethods){
        indx <- match("distclust", usedMethods)
        RES[[indx]] <- 
            parallelFct(distclustWrapper(data=data, k=k,
                                         reg=reg, regTime=regTime, fpcCtrl=fpcCtrl,
                                         funcyCtrl=funcyCtrl, ...))           
    }
    ##Method 3----------------------
    if("iterSubspace" %in% usedMethods){
        indx <- match("iterSubspace",usedMethods)
        RES[[indx]] <- 
            parallelFct(iterSubspaceWrapper(data=data, k=k, reg=reg, regTime=regTime,
                                            fpcCtrl=fpcCtrl,
                                            funcyCtrlMbc=funcyCtrl,  ...))
    }
    ##Method 4-----------
    if("funclust" %in% usedMethods){
        indx <- match("funclust",usedMethods)
        RES[[indx]] <-
            parallelFct(funclustWrapper(data=data, k=k, 
                                        reg=reg, regTime=regTime,
                                        funcyCtrlMbc=funcyCtrl,
                                        ...))
    }
    ##Method 5-----------
    if("funHDDC" %in% usedMethods){
        indx <- match("funHDDC", usedMethods)
        RES[[indx]] <-
            parallelFct(funHDDCWrapper(data=data, k=k,
                                       reg=reg, regTime=regTime,
                                       funcyCtrlMbc=funcyCtrl, ...))
    }
    ##Method 6-----------
    if("fscm" %in% usedMethods){
        indx <- match("fscm", usedMethods)
        RES[[indx]] <-
            parallelFct(fscmWrapper(data=data, k=k, reg=reg,
                                    regTime=regTime,
                                    funcyCtrlMbc=funcyCtrl, ...))
    }
    ##Method 7-----------
    if("waveclust" %in% usedMethods){
        indx <- match("waveclust", usedMethods)
        RES[[indx]] <-
            parallelFct(waveclustWrapper(data=data, k=k, reg=reg,
                                         regTime=regTime,
                                         funcyCtrlMbc=funcyCtrl, ...))
    }
    
    FRES <- new("funcyOutList")
    FRES@call <- match.call()
    
    if(parallel)
        FRES@models <- parallel::mccollect(RES)
    else
        FRES@models <- RES
    names(FRES@models) <- usedMethods

    ##Check if error appeard (only for parallel computing)----
    error <- which(sapply(FRES@models, class) == "try-error")
    if(sum(error)!=0)
        stop(paste("Method", usedMethods[error[1]], ":",
                   attributes(FRES@models[[error[1]]])$condition$message))
    
    
    allClusters <- sapply(FRES@models, function(x) x@cluster)
    allCenters <- lapply(FRES@models, function(x) x@centers)
    names(allCenters) <- colnames(allClusters) <- usedMethods
    rI <- rIMethods(methodNames=usedMethods, cls=allClusters, trueCluster=clusters)

    ##Relabel cluster output for better comparability in plots
    if(nrMethods>1){
        rel <- relabelMethods(methodNames=usedMethods, cls=allClusters,
                              ctrs=allCenters)
        allClusters <- rel$allClusters
        allCenters <- rel$allCenters
        for(i in 1:nrMethods){
            FRES@models[[i]]@cluster <- allClusters[,i]
            FRES@models[[i]]@centers <- allCenters[[i]]
            FRES@models[[i]]@correctCl <- rI[i,i]
        }
    }

    ##Warning if cluster size is smaller than 3
    smallCl <-  which(apply(allClusters, 2, function(x)
        min(table(x)))<2)

    if(length(smallCl)!=0){
        warning(paste("Method", usedMethods[smallCl],
                      "has clusters with less than 3 obervations!\n"), immediate.=TRUE)
    }
    
    accord <- accordance(cls=allClusters, relabel=FALSE)

    if(save.data)
        FRES@data <- data
    else
        FRES@data <- as.matrix(NULL)
    
    FRES@timeNr <- calcTimeNr(data, reg)
    FRES@reg <- reg
    FRES@k <- k
    FRES@methodName <- usedMethods
    FRES@allClusters <- allClusters
    FRES@randIndex <-  rI
    FRES@votedCluster <- accord$votedCluster
    FRES@accordance <- accord$accordance

    return(FRES)
}


fpcCtrlCheck <- function(fpcCtrl=NULL, data, reg){
    
    if(is.null(fpcCtrl))
        fpcCtrl <- new("fpcCtrl")

    fct.exist1 <- try(match.fun(fpcCtrl@sm1Dim),silent=TRUE)
    fct.exist2 <- try(match.fun(fpcCtrl@sm2Dim),silent=TRUE)
    
    if(class(fct.exist1)=="try-error" |
       class(fct.exist2)=="try-error")
        stop("sm1Dim and/or sm2Dim are no valid function names.")
    
    if(fpcCtrl@select=="automatic")
        {
            res <- selBw(data=data, reg=reg)
            fpcCtrl@h1Dim <- res$h1Dim
            fpcCtrl@h2Dim <- res$h2Dim
        }
    return(fpcCtrl)
}


setMethod("[[", signature(x="funcyOutList", i="ANY", j="missing"),
          function(x, i, j) x@models[[i]])


setGeneric("calcTime",
           function(object) standardGeneric("calcTime"))

setMethod("calcTime", "funcyOutList",
function(object){
    cat("\nSummary of the Calculation Time:\n")
    calcTime <- t(sapply(object@models, function(x) x@calcTime))
    rownames(calcTime) <- object@methodName
    print(calcTime)
})

setMethod("Cluster", "funcyOutList",
          function(object){
              n <- length(object@models)
              allClusters <- sapply(object@models,function(x)
                  x@cluster)
              colnames(allClusters) <- object@methodName
              if(n==1)
                  allClusters <- as.numeric(allClusters)
              allClusters
          }
          )


setGeneric("Center",
function(object) standardGeneric("Center"))

setMethod("Center", "funcyOutList",
          function(object){
              n <- length(object@models)
              allCenters <- lapply(object@models,function(x)
                  x@centers)
              names(allCenters) <- object@methodName
              if(n==1)
                  allCenters <- allCenters[[1]]
              allCenters
          }
          )


setGeneric("props",
function(object) standardGeneric("props"))

setMethod("props", "funcyOutList",
          function(object){
              a <- lapply(object@models, function(x) as.data.frame(t(x@props)))
              props <- data.frame(matrix(nrow=length(object@models),ncol=object@k))
              colnames(props) <- paste("cl", 1:object@k)
              props <- rbind.fill(a, props)
              rownames(props) <- object@methodName
              cat("\nSummary of the Cluster Proportions:\n")
              print(props)
})


setGeneric("randIndex")

setMethod("randIndex", signature(x="funcyOutList"),
function(x){
    cat("\nSummary of the Rand Indices:\n")
    print(x@randIndex)
})

setMethod("summary", "funcyOutList",
          function(object)
              {
                  outlines <- paste0(sQuote(class(object)),
                                    "\nobject with called algorithm(s):\n\n",
                                    paste(object@methodName,
                                          collapse=" "))
                  cat(writeLines(strwrap(outlines,
                                         width=0.75*getOption("width"))))
                  cat("\n\n")
                  cat("call:", deparse(object@call,0.75*getOption("width")),
                      sep="\n")
                  props(object)
                  randIndex(object)
                  calcTime(object)
              }
          )

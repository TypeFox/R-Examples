#
#  Copyright (C) 2011-2015 Christina Yassouridis
#  
#


##3 Different Formats:
##Format1: c("curveIndex", "yin", "tin")
##Format2: matrix of dimension nr_time x nr_curves - only if data is regular curves can be stored in a matrix
##Fomat3: res <- structure(list(t(data), Tin, N, isobs, t_all), names=c("Yin","Tin","N", "isobs", "t_all"))

##Default method, used if data=NULL
setGeneric("formatFuncy", function(data, format, ...) standardGeneric("formatFuncy"))


setMethod("formatFuncy", signature(data="list", format="character"),
          function(data,
                   format="Format1",
                   regTime=NULL){
              Yin <- data$Yin
              Tin <- data$Tin
              isobs <- data$isobs
               if(sum(!isobs)==0)
                   reg <- TRUE
               else
                   reg <- FALSE
               
               ##dimension
               nt <- dim(Yin)[1]
               nc <- dim(Yin)[2]
               
               switch(format,
                      "Format1"={
                          dataVec <- t(Yin)[t(isobs)==1]
                          timeVec <- t(Tin)[t(isobs)==1]
                          curveIndx <- rep(seq.int(nt),rowSums(isobs))
                          res <-
                              cbind(curveIndx, dataVec,timeVec)
                          colnames(res) <- c("curveIndx", "yin", "tin")
                      },
                      
                      "Format3"={  
                          N <- table(data[,1])
                          Yin <- Tin  <- matrix(NA, length(N), max(N))
                          isobs <- matrix(0, length(N), max(N))
                          for(i in seq.int(length(N))){
                              Yin[i,1:N[i]] <- data[,2][data[,1]==i]
                              Tin[i,1:N[i]] <- data[,3][data[,1]==i]
                              isobs[i,1:N[i]] <- 1
                          }
                          t_all <- sort(unique(data[,3]))
                          res <- structure(list(Yin, Tin, isobs, N, t_all), names
                                           = c("Yin", "Tin","isobs","N","t_all"))
                      }         
                      )
               return(res)
           }
           )

##Method used if data!=NULL and Yin and Tin and isobs are missing:
setMethod("formatFuncy", signature(data="matrix", format="character"),
          function(data,
                   format="Format1",
                   regTime=NULL){
              nt <- dim(data)[1]
              nc <- dim(data)[2]
              chf <- checkFormat(data, reformat=FALSE)
              reg <- chf$reg
              origFormat <- chf$format
              
              if(is.null(regTime))
                  regTime <-seq.int(nt)

              switch(format,
                     "Format1"={
                         if(origFormat=="Format1")
                             stop("Your data is already in Format1!")
                         dataVec <- as.vector(data)
                         timeVec <- rep(regTime, nc)
                         curveIndx <- rep(1:nc, each=nt)
                         res <-cbind(curveIndx, dataVec, timeVec)
                         colnames(res) <-  c("curveIndx", "yin",
                                             "tin")
                     },
                     
                     "Format2"={
                         if(origFormat=="Format2")
                             stop("Your data is already in Format2!")
                         mat <- NULL
                             for(i in seq.int(length(table(data[,1]))))
                                 mat <- rbind(mat,data[,2][data[,1]==i])
                             res <- t(mat)
                     },
                     
                     "Format3"={
                         if(origFormat=="Format2"){ 
                             if(is.null(regTime))
                                 Tin <- matrix(rep(regTime,nc), nrow=nc, byrow=nt)
                             else
                                 Tin <-matrix(regTime, nrow=nc, ncol=nt, byrow=TRUE)
                             N <- rep(nt, nc)
                             Yin <- t(data)
                             isobs <- matrix(1, nrow=nc, ncol=nt)
                             t_all <- sort(unique(as.vector(Tin)))
                             res <- structure(list(Yin, Tin, isobs, N, t_all),
                                              names=c("Yin","Tin","isobs","N", "t_all"))
                         }else{
                             N <- table(data[,1])
                             Yin <- Tin  <- matrix(NA, length(N), max(N))
                             isobs <- matrix(0, length(N), max(N))
                             for(i in seq.int(length(N))){
                                 Yin[i,1:N[i]] <- data[,2][data[,1]==i]
                                 Tin[i,1:N[i]] <- data[,3][data[,1]==i]
                                 isobs[i,1:N[i]] <- 1
                             }
                             t_all <- sort(unique(data[,3]))
                             res <- structure(list(Yin, Tin, isobs, N, t_all), names
                                              = c("Yin", "Tin", "isobs", "N", "t_all"))

                         }
                         
                     }
                     
                     )
              return(res)
          }
        )


checkFormat <- function(data, reformat=TRUE){
    if(dim(data)[2]!=3){
        format <- "Format2"
        reg <- TRUE
    }else if(dim(data)[2]==3 & length(unique(table(data[,1])))==1){
        format <- "Format1"
        if(reformat)
            data <- formatFuncy(data=data, format="Format2")
        reg <- TRUE
    }else if(dim(data)[2]==3 & length(unique(table(data[,1])))>1){
        format <- "Format1"
        reg <- FALSE
    }
    return(check=list(data=data, reg=reg, format=format))
}


regFuncy <- function(data,
                     timeNr=10,
                     method="project",
                     baseType=NULL,
                     nbasis=4,
                     plot=TRUE){
    chf <- checkFormat(data, reformat=FALSE)
    if(chf$reg)
        warning("Your data is already in regular format!")
    if(is.null(baseType) & method=="project")
        stop("Please choose the projection basis for method=project!")
    if(!is.null(baseType) & method%in%c("interpolate","pace"))
        warning(paste("Base type is ignored for method",method,"!"))
    
    res <- formatFuncy(data, format="Format3")
    t_all <- res$t_all
    if(plot)
        par(mfrow=c(1,2))
    time <-  makeCommonTime(data, timeNr=timeNr, plot)
    timeIndx <- unlist(lapply(time, function(x)
        which.min(abs(x-t_all))))
    
    switch(method,
           "interpolate"={
               Yin <- res$Yin;
               Tin <- res$Tin;
               isobs <- res$isobs;
               nc <- dim(Yin)[1]
               int <- lapply(seq.int(nc), function(x){
                   ys <- approxfun(Tin[x,][isobs[x,]==1],
                                   Yin[x,][isobs[x,]==1],
                                   rule=2)
                   approx <- ys(time)
               }
                             )
               result <- do.call(cbind, int)
               
           },
           "project"={
               if(baseType!="eigenbasis"){
                   if(is.null(nbasis))
                       stop("Please determine number of basis functions!")
                   res <- makeCoeffs(data=data, reg=FALSE, dimBase=nbasis,
                                     grid=t_all, pert=0.01, baseType=baseType)
                   temp <- res$base%*%t(res$coeffs)
                   result <- temp
                   if(!is.null(timeNr)){
                       result <- temp[timeIndx,]
                       time <- t_all[timeIndx]
                   }
               }else if(baseType=="eigenbasis"){
                   res <-fpca(data=data, dimBase=nbasis)
                   result <- res$yreg
                   time <- res$time
               }
           },
           "pace"={
               isAvailable <- search()
               if(! ("package:funcyOctave" %in% isAvailable))
                   stop("Please install and load package funcyOctave to use this method.")
               
               sparseFPCA <- get("sparseFPCA", mode="function")
               res <- sparseFPCA(data, time)
               time <- res$time
               result <- res$data
           }
           )
           if(plot){
               matplot(time, result, type='l', main="Regularized data")
               par(mfrow=c(1,1))
           }
           ret <- list(data=result, time=time)
           return(ret)
       }


makeCommonTime <- function(data, timeNr, plot=TRUE){
    res <- formatFuncy(data, format="Format3")
    Tin <- res$Tin
    isobs <- res$isobs
    nc <- dim(Tin)[1]
    m <- dim(Tin)[2]
    scale <- matrix(rep(seq.int(nc), each=m), nrow=nc, ncol=m, byrow=TRUE)
    
    if(plot)
        matplot(Tin, scale, type="p", col=1, pch=".", xlab="x",
                ylab="y", main="Regular time points")
    time <- sort(Tin[isobs!=0])
    nt <- length(unique(time))
    if(timeNr >= nt)
        stop(paste("Please choose a number of time points smaller than", nt))
    res <- kmeans(time, iter.max=50, timeNr)
    newTime <- sort(res$centers)
    if(plot){
        axis(3, newTime, rep("",timeNr), col.axis = "dark violet")
        abline(v=newTime, col=2, lwd=1.2)
    }
    return(newTime)
}


makeSparse <- function(Yin, Tin, isobs, time){
    n <- dim(Tin)[1]
    m <- dim(Tin)[2]
    nt <- length(time)
    indx_tin <- makeIndex(Tin, time, isobs)
    i_ind<- row(indx_tin)[as.matrix(isobs)>0]
    j_ind <- indx_tin[as.matrix(isobs)>0]

    longYin <- sparseMatrix(i_ind, j_ind,dims=c(n,nt), x=Yin[isobs>0])
    longIsobs <- sparseMatrix(i_ind, j_ind, dims=c(n,nt),
                              x=isobs[(isobs>0)])
    longTin <- sparseMatrix(i_ind, j_ind, dims=c(n,nt),
                            x=Tin[(isobs>0)])
    btin <- matrix(0,n,m)
    btin[,2:m] <- Tin[,1:(m-1)]
    Delta <- Tin - btin
    meanDelta <- apply(Delta, 1, mean)
    Delta[,1] <- rep(mean(Delta[isobs==1]),n); Delta[isobs == 0] <- 0
    longDelta <- sparseMatrix(i_ind, j_ind, dims=c(n,nt),
                              x=Delta[(isobs>0)])
    return(list(longYin=longYin, longTin=longTin, longIsobs=longIsobs, longDelta=longDelta))
}


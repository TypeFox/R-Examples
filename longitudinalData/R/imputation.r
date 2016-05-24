cat("\n###################################################################
########################## Class LongData #########################
############################ Imputation ###########################
###################################################################\n")

matrix_imputation <- function(traj,method="copyMean",lowerBound="globalMin",upperBound="globalMax"){

    ## Imputation according to the method
    trajImp <- switch(EXPR=method,
                  "crossMean"={imput_crossMean(traj)},
                  "crossMedian"={imput_crossMedian(traj)},
                  "crossHotDeck"={imput_crossHotDeck(traj)},
                  "locf"={lowerBound <- upperBound <- NA ; imput_locf(traj)},
                  "nocb"={lowerBound <- upperBound <- NA ; imput_nocb(traj)},
                  "trajMean"={lowerBound <- upperBound <- NA ; imput_trajMean(traj)},
                  "trajMedian"={lowerBound <- upperBound <- NA ; imput_trajMedian(traj)},
                  "trajHotDeck"={lowerBound <- upperBound <- NA ; imput_trajHotDeck(traj)},
                  "spline"=imput_spline(traj),
                  "linearInterpol"=,"linearInterpol.locf"=imput_linearInterpol_locf(traj),
                  "linearInterpol.global"=imput_linearInterpol_global(traj),
                  "linearInterpol.center"=imput_linearInterpol_center(traj),
                  "linearInterpol.local"=imput_linearInterpol_local(traj),
                  "linearInterpol.bisector"=imput_linearInterpol_bisector(traj),
                  "copyMean"=,"copyMean.locf"=imput_copyMean_locf(traj),
                  "copyMean.center"=imput_copyMean_center(traj), ### Pas testé !
                  "copyMean.global"=imput_copyMean_global(traj),
                  "copyMean.local"=imput_copyMean_local(traj),
                  "copyMean.bisector"=imput_copyMean_bisector(traj),
#                  "regressionInt"=imput_regressionInt(traj),
#                  "regressionExt"=imput_regressionExt(traj,predictor),
                  stop("[imputation] Method ",method," does not exists")
    )

    ## Low bounding
    if(!identical(lowerBound,NA)){
        if(length(lowerBound)==1){
            lowerBound <- rep(lowerBound,ncol(traj))
        }else{}
        lowerBound[lowerBound=="min"] <- apply(traj,2,min,na.rm=TRUE)[lowerBound=="min"]
        lowerBound[lowerBound=="globalMin"] <- min(traj,na.rm=TRUE)
        lowerBound[lowerBound=="Inf"] <- min(traj,na.rm=TRUE)

        lowerBound <- as.numeric(lowerBound)
        for(i in 1:ncol(traj)){
            trajImp[trajImp[,i]<lowerBound[i],i] <- lowerBound[i]
        }
    }else{}

    ## Upper bounding
    if(!identical(upperBound,NA)){
        if(length(upperBound)==1){
            upperBound <- rep(upperBound,ncol(traj))
        }else{}
        upperBound[upperBound=="max"] <- apply(traj,2,max,na.rm=TRUE)[upperBound=="max"]
        upperBound[upperBound=="globalMax"] <- max(traj,na.rm=TRUE)
        upperBound[upperBound=="-Inf"] <- max(traj,na.rm=TRUE)
        upperBound <- as.numeric(upperBound)
        for(i in 1:ncol(traj)){
            trajImp[trajImp[,i]>upperBound[i],i] <- upperBound[i]
        }
    }else{}

    return(trajImp)
}

setMethod(f="imputation",
    signature=c("matrix","ANY","ANY","ANY"),
    definition=matrix_imputation
)



imputation_array <- function(traj,method="copyMean",lowerBound="globalMin",upperBound="globalMax"){
    for(i in 1:dim(traj)[3]){
        traj[,,i] <- imputation(traj[,,i],method=method,lowerBound=lowerBound,upperBound=upperBound)
    }
    return(traj)
}

setMethod(f="imputation",
    signature=c("array","ANY","ANY","ANY"),
    definition=imputation_array
)


LongData_imputation <- function(traj,method="copyMean",lowerBound="globalMin",upperBound="globalMax"){
    traj@traj <- imputation(traj["traj"],method=method,lowerBound=lowerBound,upperBound=upperBound)
    return(traj)
}

setMethod(f="imputation",
    signature=c("LongData","ANY","ANY","ANY"),
    definition=LongData_imputation
)

LongData3d_imputation <- function(traj,method="copyMean",lowerBound="globalMin",upperBound="globalMax"){
    traj@traj <- imputation(traj["traj"],method=method,lowerBound=lowerBound,upperBound=upperBound)
    return(traj)
}

setMethod(f="imputation",
    signature=c("LongData3d","ANY","ANY","ANY"),
    definition=LongData3d_imputation
)

## ### Autre methode: interpolation quadratique
## x <- 1:10
## y1 <- 10+1:10*2
## y2 <- -25+1:10*3
## plot(x[1:2],y1[1:2],ylim=c(-10,30),xlim=c(0,10),type="b")
## lines(x[9:10],y2[9:10],type="b")
## y3 <- rep(NA,10)
## for (t in 2:9){
##   y3[t] <- (y2[t]*((t-2)/7)^2 +y1[t]*((9-t)/7)^2  )/ (((t-2)/7)^2+((9-t)/7)^2)
## }
## lines(x,y3,col=2,type="p")


cat("\n+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
++++++++++++++++++++++++++ Class LongData +++++++++++++++++++++++++
++++++++++++++++++++++++++ Fin Imputation +++++++++++++++++++++++++
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n")

